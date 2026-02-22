# AI Explorer Authentication Plan

## Goal

Protect the AI Explorer (chat) functionality behind authentication so only pre-approved users can use it. The rest of the Perturbation Catalogue (search, targets, datasets, downloads) remains public.

---

## Current State

- **Zero authentication** anywhere in the stack — no middleware, no tokens, no cookies, no login
- `POST /v1/chat/stream` is completely open — anyone with the URL can call Gemini (costs money)
- Chat sessions are keyed by random UUID with no user binding
- CORS is `allow_origins=["*"]`
- Frontend chat.js sends `fetch()` with only `Content-Type` header — no auth
- Frontend and backend are separate Cloud Run services; in production they share `www.ebi.ac.uk/perturbation-catalogue` via a reverse proxy

---

## Architecture Decision: JWT + Authorization Header

**Why JWT over sessions/cookies:**
- Frontend (Dash/Gunicorn) and backend (FastAPI/Uvicorn) are separate services on separate Cloud Run instances
- chat.js makes direct browser→backend fetch calls (not proxied through frontend)
- JWT in `Authorization: Bearer <token>` header works regardless of same-origin/cross-origin setup
- No need for CSRF protection (no cookies)
- Stateless — no server-side session store needed (backend already has in-memory `_sessions` for chat history, not auth)

**Why sessionStorage for token persistence:**
- Cleared automatically when tab/browser closes
- Acceptable for a research tool with pre-approved users
- Simple to implement — chat.js reads directly via `sessionStorage.getItem()`

**Trade-off:** sessionStorage is accessible to JavaScript (XSS risk). Acceptable here because:
1. The app has no user-generated content (no stored XSS vector)
2. Users are pre-approved researchers, not the general public
3. Can upgrade to HTTP-only cookies later if needed

---

## What Needs Protection

| Endpoint | Auth? | Reason |
|----------|-------|--------|
| `POST /v1/chat/stream` | **Yes** | Calls Gemini, costs money |
| `GET /v1/auth/login`, `POST /v1/auth/login` | No | Login endpoints |
| `GET /v1/auth/me` | Yes | Returns current user info |
| `GET /`, `/health` | No | Infrastructure |
| `GET /summary`, `/search`, `/dataset/*` | No | Public catalogue data |
| `GET /v1/{modality}/search`, `/download` | No | Public catalogue data |
| Frontend `/chat` page | **Yes** | Gate the UI itself |
| All other frontend pages | No | Public |

---

## GCP Actions (Manual, Before Code Changes)

### 1. Create the `users` Table in Cloud SQL

Connect to the existing Cloud SQL PostgreSQL instance and run:

```sql
CREATE TABLE IF NOT EXISTS users (
    id SERIAL PRIMARY KEY,
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL,
    display_name VARCHAR(100) NOT NULL,
    is_active BOOLEAN DEFAULT TRUE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    last_login TIMESTAMP WITH TIME ZONE
);

-- Index for login lookups
CREATE INDEX idx_users_email ON users (email);
```

### 2. Create a JWT Signing Secret in Secret Manager

```bash
# Generate a random 256-bit key
openssl rand -base64 32

# Store it in Secret Manager
gcloud secrets create jwt-signing-key \
    --project=prj-ext-dev-pertcat-437314 \
    --replication-policy=automatic

echo -n "<the-generated-key>" | gcloud secrets versions add jwt-signing-key \
    --project=prj-ext-dev-pertcat-437314 \
    --data-file=-
```

### 3. Update Backend Cloud Run Environment Variables

Add to the backend Cloud Run service:

| Variable | Value | Source |
|----------|-------|--------|
| `JWT_SECRET` | (from Secret Manager) | Secret Manager reference: `projects/prj-ext-dev-pertcat-437314/secrets/jwt-signing-key/versions/latest` |
| `JWT_EXPIRY_HOURS` | `168` (7 days) | Plain text — adjust as needed |

```bash
# Using gcloud CLI (or via Console UI)
gcloud run services update perturbation-catalogue-be \
    --region=europe-west2 \
    --update-secrets=JWT_SECRET=jwt-signing-key:latest \
    --update-env-vars=JWT_EXPIRY_HOURS=168
```

### 4. Create Initial User Accounts

After deploying the code changes, use the `manage_users.py` CLI script (described below) to create accounts:

```bash
# Connect to Cloud SQL (via Cloud SQL Auth Proxy or from a VM with access)
python manage_users.py create --email researcher@ebi.ac.uk --name "Jane Doe"
# Prompts for password interactively, prints confirmation
```

---

## Backend Changes

### New File: `be/auth.py`

Core auth module with:

```
auth.py
├── Settings (JWT_SECRET, JWT_EXPIRY_HOURS from env)
├── Password hashing (passlib + bcrypt)
│   ├── hash_password(plain) -> str
│   └── verify_password(plain, hash) -> bool
├── JWT handling (PyJWT)
│   ├── create_token(user_id, email) -> str  (includes exp claim)
│   └── decode_token(token) -> dict          (raises on invalid/expired)
├── FastAPI dependencies
│   ├── get_current_user(token) -> User      (raises 401 if invalid)
│   └── get_optional_user(token) -> User|None (for mixed endpoints)
├── Router: /v1/auth
│   ├── POST /login   — validates email+password against DB, returns {token, user}
│   ├── GET  /me      — returns current user info (requires valid token)
│   └── POST /logout  — optional, for audit logging
└── Pydantic models
    ├── LoginRequest(email, password)
    ├── LoginResponse(token, user)
    └── UserInfo(id, email, display_name)
```

**Dependencies:** Add to `be/requirements.txt`:
```
PyJWT>=2.8.0
passlib[bcrypt]>=1.7.4
```

**Key implementation notes:**
- Use `OAuth2PasswordBearer(tokenUrl="/v1/auth/login")` as the token extractor — this also enables the "Authorize" button in FastAPI's Swagger UI for testing
- `get_current_user` reads the `Authorization: Bearer <token>` header, decodes JWT, looks up user in `db_pools["pg"]`, checks `is_active`, returns user dict
- `create_token` sets `exp` to `now + JWT_EXPIRY_HOURS`, includes `sub` (user ID) and `email` claims
- Password hashing uses `bcrypt` with automatic salting via passlib

### Modify: `be/main.py`

1. **Add settings fields:**
```python
class Settings(BaseSettings):
    # ... existing fields ...
    jwt_secret: str = ""          # Required for auth
    jwt_expiry_hours: int = 168   # 7 days default
```

2. **Register auth router:**
```python
from auth import router as auth_router, configure as configure_auth

app.include_router(auth_router)  # adds /v1/auth/* routes
```

3. **Pass settings to auth module in lifespan:**
```python
configure_auth(
    jwt_secret=settings.jwt_secret,
    jwt_expiry_hours=settings.jwt_expiry_hours,
)
```

4. **Tighten CORS** (replace wildcard with actual origins):
```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://www.ebi.ac.uk",
        "https://perturbation-catalogue-fe-*.europe-west2.run.app",  # dev
        "http://localhost:8050",   # local dev
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### Modify: `be/ai_chat.py`

Add auth dependency to the chat endpoint:

```python
from auth import get_current_user

@router.post("/stream")
async def chat_stream(request: ChatRequest, user: dict = Depends(get_current_user)):
    # user is now available — bind session to user
    session_id = request.session_id or str(uuid.uuid4())
    # ... rest unchanged ...
```

**Optionally:** Tie `_sessions` to user IDs instead of random UUIDs, so chat history persists across page refreshes for the same user.

### New File: `be/manage_users.py`

Standalone CLI script for user administration (not part of the FastAPI app):

```
Usage:
  python manage_users.py create --email user@ebi.ac.uk --name "Display Name"
  python manage_users.py list
  python manage_users.py disable --email user@ebi.ac.uk
  python manage_users.py enable --email user@ebi.ac.uk
  python manage_users.py reset-password --email user@ebi.ac.uk
```

Reads `PG_*` env vars (same as the backend), connects directly to PostgreSQL. Uses passlib for password hashing (same as auth.py). Prompts for passwords interactively via `getpass`.

---

## Frontend Changes

### New File: `fe/pages/login.py`

A login page registered at `/login`:

```python
dash.register_page(__name__, path="/login", name="Login")
```

Layout:
- Centered card with email + password fields + "Sign in" button
- Error message area for failed logins
- Styled consistently with the existing app (same green `#007B53` primary, Bootstrap)

Login flow uses a **clientside_callback** (runs in browser JavaScript, not server-side Python):
1. User fills form, clicks "Sign in"
2. Clientside callback `fetch`es `POST /v1/auth/login` with `{email, password}`
3. On success: stores `{token, user}` in sessionStorage, redirects to `/chat` via `window.location`
4. On failure: shows error message ("Invalid credentials" / "Account disabled")

```python
app.clientside_callback(
    """
    function(n_clicks, email, password) {
        if (!n_clicks) return [window.dash_clientside.no_update, ""];
        var baseUrl = ...; // from hidden div or env
        return fetch(baseUrl + "/v1/auth/login", {
            method: "POST",
            headers: {"Content-Type": "application/json"},
            body: JSON.stringify({email: email, password: password})
        })
        .then(r => r.ok ? r.json() : Promise.reject(r))
        .then(data => {
            sessionStorage.setItem("auth_token", data.token);
            sessionStorage.setItem("auth_user", JSON.stringify(data.user));
            window.location.href = "/perturbation-catalogue/chat";
            return [window.dash_clientside.no_update, ""];
        })
        .catch(() => [window.dash_clientside.no_update, "Invalid email or password"]);
    }
    """,
    [Output("login-redirect", "data"), Output("login-error", "children")],
    Input("login-btn", "n_clicks"),
    [State("login-email", "value"), State("login-password", "value")],
    prevent_initial_call=True,
)
```

### Modify: `fe/assets/chat.js`

Three changes:

**1. Auth check on load** — at the top of the initialization, check for token:
```javascript
// Before anything else
var token = sessionStorage.getItem("auth_token");
if (!token) {
    window.location.href = "/perturbation-catalogue/login";
}
```

**2. Send token with fetch** — add Authorization header:
```javascript
fetch(url, {
    method: "POST",
    headers: {
        "Content-Type": "application/json",
        "Authorization": "Bearer " + sessionStorage.getItem("auth_token")
    },
    body: body,
    signal: currentController.signal
})
```

**3. Handle 401 response** — redirect to login on expired/invalid token:
```javascript
.then(function (response) {
    if (response.status === 401) {
        sessionStorage.removeItem("auth_token");
        sessionStorage.removeItem("auth_user");
        window.location.href = "/perturbation-catalogue/login";
        return;
    }
    if (!response.ok) {
        throw new Error("HTTP " + response.status);
    }
    // ... rest of SSE handling ...
})
```

### Modify: `fe/pages/chat.py`

Add a hidden div for auth user display (optional — show username in the sidebar):
```python
html.Div(id="chat-auth-user", style={"display": "none"}),
```

### Modify: `fe/app.py`

Add a logout link in the header nav (visible only on chat page, handled client-side):
```python
html.A(
    "Logout",
    id="logout-link",
    href="#",
    className="header-link",
    style={"display": "none"},  # chat.js shows this when authenticated
),
```

---

## Implementation Order

### Step 1: Backend Auth Module
1. Create `be/auth.py` with JWT/password utilities, FastAPI dependencies, and login/me/logout endpoints
2. Add `PyJWT` and `passlib[bcrypt]` to `be/requirements.txt`
3. Add `jwt_secret` and `jwt_expiry_hours` to `Settings` in `main.py`
4. Register the auth router in `main.py`
5. Test login flow via Swagger UI (`/docs`)

### Step 2: Protect Chat Endpoint
1. Add `Depends(get_current_user)` to `chat_stream` in `ai_chat.py`
2. Verify: unauthenticated requests get 401, authenticated requests work
3. Tighten CORS origins in `main.py`

### Step 3: User Management CLI
1. Create `be/manage_users.py`
2. Create test accounts locally
3. Verify login→chat flow end-to-end via curl

### Step 4: Frontend Login Page
1. Create `fe/pages/login.py` with form + clientside_callback
2. Style login page to match the existing app design
3. Register callbacks in `fe/app.py`

### Step 5: Frontend Auth Integration
1. Update `chat.js`: auth check on load, Authorization header on fetch, 401 handling
2. Add logout link in the header
3. Test full flow: login → chat → logout → redirect to login

### Step 6: GCP Deployment
1. Create `users` table in Cloud SQL (SQL above)
2. Store JWT secret in Secret Manager
3. Update Cloud Run env vars
4. Create initial user accounts via `manage_users.py`
5. Deploy backend, then frontend
6. Verify end-to-end on dev environment

---

## Optional Enhancements (Post-Launch)

**Usage tracking:** Add a `chat_usage` table to log per-user request counts and token consumption. Useful for monitoring costs.

```sql
CREATE TABLE chat_usage (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id),
    session_id VARCHAR(36),
    message_length INTEGER,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);
```

**Rate limiting:** Add per-user rate limits (e.g., 100 requests/day) using an in-memory counter or Redis. Prevents runaway Gemini costs.

**Token refresh:** Current plan uses 7-day tokens. If shorter expiry is needed, add a `POST /v1/auth/refresh` endpoint that issues a new token given a valid (non-expired) existing token.

**Upgrade to HTTP-only cookies:** If same-origin routing is guaranteed (both services behind `www.ebi.ac.uk`), switch from sessionStorage to HTTP-only secure cookies for stronger XSS protection. Requires adding CSRF protection.

**Admin UI:** Replace CLI `manage_users.py` with a protected admin page in the frontend for user management. Lower priority since the user base is small and pre-approved.
