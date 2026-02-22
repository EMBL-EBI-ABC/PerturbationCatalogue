"""Authentication module for AI Explorer access control.

Provides JWT-based auth with pre-approved user accounts stored in PostgreSQL.
Only the AI chat endpoint requires authentication; public catalogue endpoints remain open.
"""

import logging
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, Optional

import bcrypt
import jwt
from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer
from pydantic import BaseModel

from data_query import db_pools

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/auth", tags=["Authentication"])

# --- Configuration (set via configure()) ---

_config: Dict[str, Any] = {
    "jwt_secret": "",
    "jwt_expiry_hours": 168,
}

# Token extractor — also enables the "Authorize" button in Swagger UI
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/v1/auth/login", auto_error=False)


def configure(jwt_secret: str, jwt_expiry_hours: int = 168):
    if not jwt_secret:
        logger.warning("JWT_SECRET is not set — authentication will not work")
    _config["jwt_secret"] = jwt_secret
    _config["jwt_expiry_hours"] = jwt_expiry_hours


# --- Password hashing ---


def hash_password(plain: str) -> str:
    return bcrypt.hashpw(plain.encode(), bcrypt.gensalt()).decode()


def verify_password(plain: str, hashed: str) -> bool:
    return bcrypt.checkpw(plain.encode(), hashed.encode())


# --- JWT handling ---


def create_token(user_id: int, email: str) -> str:
    payload = {
        "sub": str(user_id),
        "email": email,
        "exp": datetime.now(timezone.utc)
        + timedelta(hours=_config["jwt_expiry_hours"]),
    }
    return jwt.encode(payload, _config["jwt_secret"], algorithm="HS256")


def decode_token(token: str) -> dict:
    return jwt.decode(token, _config["jwt_secret"], algorithms=["HS256"])


# --- Pydantic models ---


class LoginRequest(BaseModel):
    email: str
    password: str


class UserInfo(BaseModel):
    id: int
    email: str
    display_name: str


class LoginResponse(BaseModel):
    token: str
    user: UserInfo


# --- FastAPI dependencies ---


async def get_current_user(token: Optional[str] = Depends(oauth2_scheme)) -> dict:
    """Require a valid JWT. Returns the user dict or raises 401."""
    if not token:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Not authenticated",
            headers={"WWW-Authenticate": "Bearer"},
        )
    try:
        payload = decode_token(token)
    except jwt.ExpiredSignatureError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Token expired",
            headers={"WWW-Authenticate": "Bearer"},
        )
    except jwt.InvalidTokenError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid token",
            headers={"WWW-Authenticate": "Bearer"},
        )

    pool = db_pools.get("pg")
    if not pool:
        raise HTTPException(status_code=503, detail="Database unavailable")

    async with pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT id, email, display_name, is_active FROM users WHERE id = $1",
            int(payload["sub"]),
        )
    if not row:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail="User not found"
        )
    if not row["is_active"]:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail="Account disabled"
        )
    return dict(row)


async def get_optional_user(
    token: Optional[str] = Depends(oauth2_scheme),
) -> Optional[dict]:
    """Like get_current_user but returns None instead of raising 401."""
    if not token:
        return None
    try:
        return await get_current_user(token)
    except HTTPException:
        return None


# --- Endpoints ---


@router.post("/login", response_model=LoginResponse)
async def login(request: LoginRequest):
    pool = db_pools.get("pg")
    if not pool:
        raise HTTPException(status_code=503, detail="Database unavailable")

    async with pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT id, email, display_name, password_hash, is_active FROM users WHERE email = $1",
            request.email,
        )

        if not row or not verify_password(request.password, row["password_hash"]):
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Invalid email or password",
            )

        if not row["is_active"]:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED, detail="Account disabled"
            )

        # Update last_login
        await conn.execute(
            "UPDATE users SET last_login = NOW() WHERE id = $1", row["id"]
        )

    token = create_token(row["id"], row["email"])
    return LoginResponse(
        token=token,
        user=UserInfo(
            id=row["id"], email=row["email"], display_name=row["display_name"]
        ),
    )


@router.get("/me", response_model=UserInfo)
async def me(user: dict = Depends(get_current_user)):
    return UserInfo(
        id=user["id"], email=user["email"], display_name=user["display_name"]
    )
