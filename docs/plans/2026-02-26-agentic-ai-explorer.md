# Agentic AI Explorer: Multi-Step Planning & Execution

**Date:** 2026-02-26
**Author:** Alexey + Claude
**Status:** Draft / For Discussion

---

## 1. Executive Summary

Evolve the AI Explorer from a **1:1 request-to-tool mapper** into a **multi-step agentic system** that can plan, execute multiple tools, and synthesize answers from combined results — while showing intermediate results in the dashboard as they arrive.

Currently, Gemini Flash 2.5 receives a user question, decides which single tool (or small batch of tools) to call, gets the result, and responds. This works for direct questions ("show me a volcano plot for TP53") but fails for sophisticated queries like:

> "Compare the druggability of TP53 vs BRCA1, show their protein structures side by side, and tell me which one has more significant perturbation effects across all datasets"

This requires calling 6+ tools across multiple categories, collecting all results, and writing a comparative synthesis — something the current single-pass architecture cannot do well.

---

## 2. Current Architecture (As-Is)

```
User message
  → Gemini 2.5 Flash (single model, temperature=0.3)
    → Tool-calling loop (max 12 iterations)
      → Each iteration: LLM decides tool(s) → execute → feed results back
    → Final text response
  → SSE stream to frontend (thinking → tool_call → visualization → text → done)
```

**Key characteristics:**
- Model: `gemini-2.5-flash` (hardcoded default, server-configured via env var)
- Loop: Up to 12 iterations, but Flash model rarely uses more than 2-3
- No explicit planning phase — the model implicitly decides what to do
- All tool results fed back into context window (can get large)
- Non-streaming LLM calls (full response per iteration, not token-by-token)
- No user visibility into the "plan" — just spinner status updates
- 23 internal tools + Open Targets MCP tools

**What works well:**
- SSE streaming with intermediate visualization rendering
- Server-side pagination for tables/heatmaps (bypasses LLM)
- Visualization persistence in PostgreSQL
- Session management with full history

**Limitations for complex queries:**
- Flash model doesn't plan well — it's optimized for speed, not deep reasoning
- No explicit plan shown to user — they can't see or guide the approach
- 12-iteration cap can be hit on complex multi-tool queries
- Large tool results consume context window quickly
- No way to prioritize or parallelize tool calls
- Single model — can't use a smarter model for planning + cheaper model for execution

---

## 3. Proposed Architecture (To-Be)

### 3.1 Two-Phase Agent: Plan → Execute

```
User message
  ┌─────────────────────────────────────────────┐
  │ PHASE 1: PLANNING (Gemini 3.1 Pro)          │
  │                                             │
  │ Input: user question + available tools list │
  │ Output: structured execution plan (JSON)    │
  │   - ordered list of steps                   │
  │   - each step: tool name, args, purpose     │
  │   - dependency graph (which steps need      │
  │     results from previous steps)            │
  │   - final synthesis instructions            │
  │                                             │
  │ → SSE: emit plan to frontend for display    │
  └─────────────────┬───────────────────────────┘
                    │
  ┌─────────────────▼───────────────────────────┐
  │ PHASE 2: EXECUTION (Gemini 2.5 Flash or     │
  │          same model based on user choice)   │
  │                                             │
  │ For each step in plan:                      │
  │   1. Execute tool call                      │
  │   2. Emit visualization to dashboard (SSE)  │
  │   3. Collect result summary for synthesis   │
  │                                             │
  │ After all steps complete:                   │
  │   → Feed all summaries to LLM               │
  │   → Generate synthesis response             │
  │   → SSE: emit final text                    │
  └─────────────────────────────────────────────┘
```

### 3.2 Model Selection

Add a model selector dropdown to the frontend, allowing users to choose:

| Model ID | Display Name | Best For | Pricing |
|---|---|---|---|
| `gemini-2.5-flash` | Gemini 2.5 Flash | Quick questions, simple lookups | Cheapest |
| `gemini-3-flash-preview` | Gemini 3 Flash | Balanced speed + reasoning | Mid |
| `gemini-3.1-pro-preview` | Gemini 3.1 Pro | Complex multi-step analysis | ~$2/$12 per 1M tokens |

**Behavior by model:**
- **Flash models**: Use the current single-pass architecture (no planning phase). Fast, cheap, good for "show me X" type questions.
- **3.1 Pro**: Activates the agentic planning mode. Shows the plan to the user, executes step-by-step, synthesizes at the end.

This means the planning overhead only kicks in when a capable model is selected — Flash users get the same fast experience they have today.

### 3.3 "Simple" vs "Agentic" Mode Decision

Rather than always planning, the system should decide based on query complexity. Two approaches (recommend Option A):

**Option A — Model-driven:** If user selects a Pro model, always use the planning flow. If Flash, always use direct flow. Simple, predictable.

**Option B — Auto-detect:** A lightweight classifier (or the LLM itself) decides whether planning is needed based on the query. More magical, but harder to debug and explain.

**Recommendation: Option A.** It's transparent — users who want deeper analysis explicitly choose the Pro model. We can always add auto-detect later.

---

## 4. Detailed Design

### 4.1 Backend Changes (`be/ai_chat.py`)

#### 4.1.1 New Request Model

```python
class ChatRequest(BaseModel):
    message: str
    session_id: Optional[str] = None
    model: Optional[str] = None  # NEW: user-selected model override
```

The `model` field lets the frontend pass the user's model choice. If `None`, fall back to the server default (`GEMINI_MODEL` env var).

#### 4.1.2 Available Models Endpoint

```python
@router.get("/models")
async def list_models():
    """Return available models for the frontend dropdown."""
    return [
        {
            "id": "gemini-2.5-flash",
            "name": "Gemini 2.5 Flash",
            "description": "Fast responses for simple questions",
            "supports_planning": False,
        },
        {
            "id": "gemini-3-flash-preview",
            "name": "Gemini 3 Flash",
            "description": "Balanced speed and reasoning",
            "supports_planning": False,
        },
        {
            "id": "gemini-3.1-pro-preview",
            "name": "Gemini 3.1 Pro",
            "description": "Deep analysis with multi-step planning",
            "supports_planning": True,
        },
    ]
```

#### 4.1.3 Planning Phase

When a `supports_planning` model is selected, the `generate()` function adds a planning phase before the execution loop:

```python
# Phase 1: Planning
if model_supports_planning:
    plan = await _generate_plan(client, model_name, message, history)
    yield _sse_event("plan", plan)  # NEW SSE event type

    # Phase 2: Execute plan steps
    for step in plan["steps"]:
        yield _sse_event("plan_step", {"step_id": step["id"], "status": "running"})
        result = await _execute_plan_step(client, model_name, step, collected_results)
        collected_results.append(result)
        yield _sse_event("plan_step", {"step_id": step["id"], "status": "done"})

    # Phase 3: Synthesis
    yield _sse_event("thinking", {"status": "Synthesizing results..."})
    synthesis = await _synthesize(client, model_name, message, collected_results)
    yield _sse_event("text", {"content": synthesis})
else:
    # Current single-pass flow (unchanged)
    ...
```

#### 4.1.4 Plan Structure

The planning LLM call uses a dedicated system prompt and returns structured JSON:

```json
{
  "plan_summary": "Compare TP53 and BRCA1 across druggability, structure, and perturbation data",
  "steps": [
    {
      "id": 1,
      "tool": "find_datasets_for_target",
      "args": {"gene_name": "TP53"},
      "purpose": "Find all datasets containing TP53 perturbations",
      "depends_on": []
    },
    {
      "id": 2,
      "tool": "find_datasets_for_target",
      "args": {"gene_name": "BRCA1"},
      "purpose": "Find all datasets containing BRCA1 perturbations",
      "depends_on": []
    },
    {
      "id": 3,
      "tool": "get_druggability",
      "args": {"gene_name": "TP53"},
      "purpose": "Assess TP53 druggability from Pharos",
      "depends_on": []
    },
    {
      "id": 4,
      "tool": "get_druggability",
      "args": {"gene_name": "BRCA1"},
      "purpose": "Assess BRCA1 druggability from Pharos",
      "depends_on": []
    },
    {
      "id": 5,
      "tool": "get_protein_structure",
      "args": {"gene_name": "TP53"},
      "purpose": "Show TP53 protein structure",
      "depends_on": []
    },
    {
      "id": 6,
      "tool": "get_protein_structure",
      "args": {"gene_name": "BRCA1"},
      "purpose": "Show BRCA1 protein structure",
      "depends_on": []
    },
    {
      "id": 7,
      "tool": "create_volcano_plot",
      "args": {"gene_name": "TP53", "dataset_id": "$FROM_STEP_1"},
      "purpose": "Visualize TP53 perturbation effects",
      "depends_on": [1]
    }
  ],
  "synthesis_instructions": "Compare druggability profiles, highlight structural differences, and summarize which gene shows stronger perturbation effects with statistical evidence."
}
```

**Key design decisions:**
- `depends_on` enables parallel execution of independent steps (1, 2, 3, 4, 5, 6 can all run concurrently)
- `$FROM_STEP_N` placeholders allow dynamic arg resolution from prior step results
- `synthesis_instructions` guide the final synthesis LLM call

#### 4.1.5 Plan Step Execution

```python
async def _execute_plan_step(client, model_name, step, prior_results):
    """Execute a single plan step, resolving dependencies."""
    # 1. Resolve $FROM_STEP_N placeholders in args
    resolved_args = _resolve_step_args(step["args"], prior_results)

    # 2. Execute the tool (reuse existing tool dispatch logic)
    result = await _execute_tool(step["tool"], resolved_args)

    # 3. If tool produces a visualization, emit it via SSE
    #    (reuse existing viz emission logic)

    # 4. Return summary for synthesis
    return {
        "step_id": step["id"],
        "tool": step["tool"],
        "purpose": step["purpose"],
        "result_summary": _summarize_result(result),  # truncate large results
    }
```

#### 4.1.6 Parallel Step Execution

Steps with no unresolved dependencies can run concurrently:

```python
async def _execute_plan(steps, ...):
    completed = {}
    pending = list(steps)

    while pending:
        # Find steps whose dependencies are all completed
        ready = [s for s in pending if all(d in completed for d in s["depends_on"])]

        # Execute ready steps in parallel
        results = await asyncio.gather(
            *[_execute_plan_step(s, completed) for s in ready]
        )

        for step, result in zip(ready, results):
            completed[step["id"]] = result
            pending.remove(step)
```

#### 4.1.7 Thinking Budget Configuration

For Gemini 3.1 Pro, enable the thinking budget at the "HIGH" level during the planning phase:

```python
config = types.GenerateContentConfig(
    system_instruction=PLANNING_SYSTEM_INSTRUCTION,
    tools=tools,
    temperature=0.3,
    thinking_config=types.ThinkingConfig(thinking_budget=10000),  # HIGH
)
```

### 4.2 New SSE Event Types

| Event | Payload | Frontend Behavior |
|---|---|---|
| `plan` | `{plan_summary, steps: [{id, tool, purpose, depends_on}]}` | Show plan card in dashboard with step list |
| `plan_step` | `{step_id, status: "running"/"done"/"error"}` | Update step status indicators in plan card |
| `model_info` | `{model_id, model_name, supports_planning}` | Show active model badge in UI |

### 4.3 Frontend Changes (`fe/assets/chat.js`)

#### 4.3.1 Model Selector Dropdown

Add a dropdown next to the chat input area:

```
┌──────────────────────────────────────────────────────┐
│ [Model: Gemini 2.5 Flash ▼]                          │
│ ┌──────────────────────────────────────────────────┐ │
│ │ Type your question...                        [→] │ │
│ └──────────────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────┘
```

- Fetch available models from `GET /v1/chat/models` on page load
- Store selected model in `sessionStorage`
- Pass `model` field in the `ChatRequest` POST body
- Show a subtle badge when Pro model is selected: "Deep analysis mode"

#### 4.3.2 Plan Card Visualization

When a `plan` SSE event arrives, render a plan card in the dashboard:

```
┌─────────────────────────────────────────────┐
│ 🔬 Execution Plan                            │
│                                              │
│ "Compare TP53 and BRCA1 across druggability, │
│  structure, and perturbation data"           │
│                                              │
│ ✅ 1. Find TP53 datasets                     │
│ ✅ 2. Find BRCA1 datasets                    │
│ ⏳ 3. Get TP53 druggability                  │
│ ⏳ 4. Get BRCA1 druggability                 │
│ ○  5. Show TP53 protein structure            │
│ ○  6. Show BRCA1 protein structure           │
│ ○  7. Visualize TP53 perturbation effects    │
│                                              │
│ Progress: 2/7 steps complete                 │
└─────────────────────────────────────────────┘
```

- Steps update in real-time as `plan_step` events arrive
- The plan card stays pinned at the top of the dashboard
- Other visualization tiles appear below as tools produce them

#### 4.3.3 Updated `handleSSEEvent()`

```javascript
case "plan":
    renderPlanCard(data);
    break;
case "plan_step":
    updatePlanStep(data.step_id, data.status);
    break;
case "model_info":
    updateModelBadge(data);
    break;
```

### 4.4 Frontend Layout Changes (`fe/pages/chat.py`)

Add model selector component to the chat input area:

```python
# In the chat-input-area div
html.Div([
    dcc.Dropdown(
        id="model-selector",
        className="model-selector",
        clearable=False,
    ),
], className="model-selector-wrapper"),
```

The dropdown is populated client-side via JavaScript (fetch from `/v1/chat/models`), not via Dash callbacks — consistent with how the rest of the chat UI works (pure JS, no Dash callbacks).

### 4.5 System Prompt Changes

Two separate system prompts:

**PLANNING_SYSTEM_INSTRUCTION** — Used only in Phase 1:
```
You are a planning agent for the Perturbation Catalogue AI Explorer.
Given a user question and a list of available tools, create an execution plan.

Output a JSON plan with:
- plan_summary: 1-sentence description of the approach
- steps: ordered list of tool calls with dependencies
- synthesis_instructions: how to combine results into a final answer

Rules:
- Mark steps as independent (depends_on: []) when they don't need prior results
- Use $FROM_STEP_N placeholders for dynamic arguments
- Minimize total steps — don't add unnecessary tool calls
- Always call find_datasets_for_target before dataset-specific tools
- Prefer visualization tools when the user asks to "show" or "compare"
```

**SYSTEM_INSTRUCTION** — Existing prompt, used for execution and single-pass mode (unchanged).

**SYNTHESIS_INSTRUCTION** — Used in Phase 3:
```
You are synthesizing results from multiple tool calls to answer the user's question.
You have the original question, the execution plan, and summaries of all tool results.
Write a clear, comprehensive answer that references the visualizations shown in the dashboard.
Do NOT include raw data — the user can see it in the dashboard visualizations.
Focus on insights, comparisons, and actionable conclusions.
```

---

## 5. Migration Strategy

### Phase 1: Model Selector (1-2 days)
- Add `model` field to `ChatRequest`
- Add `GET /v1/chat/models` endpoint
- Add model dropdown to frontend
- Pass selected model to backend
- **No behavior change** — all models still use single-pass flow

### Phase 2: Planning Infrastructure (2-3 days)
- Implement `_generate_plan()` with planning system prompt
- Implement `_execute_plan()` with dependency resolution and parallel execution
- Implement `_synthesize()` for final answer generation
- Add `plan` and `plan_step` SSE events
- Wire up: if model `supports_planning`, use new flow

### Phase 3: Frontend Plan UI (1-2 days)
- Implement plan card renderer
- Implement real-time step status updates
- Style the plan card
- Handle edge cases (plan errors, step failures)

### Phase 4: Refinement (ongoing)
- Tune planning system prompt based on real usage
- Add plan caching/reuse for similar queries
- Consider letting users edit/approve the plan before execution
- Add cost estimation display for Pro model queries

---

## 6. Feasibility Assessment

### Is this feasible? **Yes, strongly.**

The current architecture already has all the building blocks:

1. **Tool-calling loop already exists** (`ai_chat.py:4212-4529`) — we're refactoring it, not building from scratch
2. **SSE streaming already works** — we just add new event types (`plan`, `plan_step`)
3. **Visualization rendering already works** — plan execution reuses existing tool dispatch and viz emission
4. **Session persistence already works** — plans are just another part of the conversation
5. **`google-genai` SDK supports all models** — same `client.aio.models.generate_content()` call, different model string
6. **Gemini 3.1 Pro is built for this** — 80.6% SWE-Bench, thinking budget support, optimized for multi-step tool use

### Risks and mitigations

| Risk | Mitigation |
|---|---|
| Planning phase adds latency (3.1 Pro is slower than Flash) | User explicitly opts in by selecting Pro model. Show plan card immediately so they see progress. |
| Plan may be wrong / hallucinate tool args | Validate plan structure before execution. Fall back to single-pass if plan is malformed. |
| Cost of 3.1 Pro ($2/$12 per 1M tokens) | Show model selection clearly. Consider usage limits per user. |
| Gemini 3.1 Pro is still in preview | Fall back gracefully. Model list is server-configured, easy to swap. |
| Parallel tool execution may hit rate limits | Use `asyncio.Semaphore` to limit concurrency (e.g., max 3 parallel tool calls). |
| Large plans may exceed context window | Cap plans at 10 steps. Summarize tool results before feeding to synthesis. |

---

## 7. Model Selector: Should We Add It?

**Yes.** Even independent of the agentic planning feature, a model selector is valuable:

- **User empowerment:** Power users can choose Pro for complex analysis, casual users stay on Flash
- **Cost transparency:** Users understand that Pro costs more / is slower
- **Future-proofing:** New models can be added to the dropdown without code changes (just update the `/models` endpoint)
- **A/B comparison:** Users can ask the same question with different models and compare

**Implementation is minimal:** one new endpoint, one dropdown, one field added to `ChatRequest`.

---

## 8. Open Questions

1. **Should users be able to edit/approve the plan before execution?** This would add a confirmation step where the plan card shows "Execute this plan?" with approve/modify buttons. More control, but adds friction.

2. **Should we store plans in the database?** Currently tool calls aren't persisted individually. Storing plans would enable analytics on what types of queries trigger planning, and allow resuming interrupted plans.

3. **Should the synthesis phase also have access to tool calling?** Sometimes the synthesis might need one more data point. Could allow 1-2 additional tool calls during synthesis. But this adds complexity.

4. **Rate limiting per model tier?** Pro queries are expensive. Should we limit Pro usage (e.g., N queries/day per user) or let infrastructure costs determine this?

5. **Fallback strategy for Gemini 3.1 Pro preview instability?** If 3.1 Pro fails, should we retry with 3 Flash, or just error? Recommend: error with clear message, don't silently downgrade.

---

## 9. References

- [Gemini API Models](https://ai.google.dev/gemini-api/docs/models) — full model list and capabilities
- [Gemini 3.1 Pro Preview docs](https://ai.google.dev/gemini-api/docs/models/gemini-3.1-pro-preview) — model card and API details
- [Gemini 3.1 Pro on Vertex AI](https://docs.cloud.google.com/vertex-ai/generative-ai/docs/models/gemini/3-1-pro) — Vertex AI integration
- [Gemini Tools & Agents guide](https://ai.google.dev/gemini-api/docs/tools) — function calling and agentic patterns
- [Gemini 3 Developer Guide](https://ai.google.dev/gemini-api/docs/gemini-3) — Gemini 3 series overview
- [Gemini 3.1 Pro announcement](https://blog.google/innovation-and-ai/models-and-research/gemini-models/gemini-3-1-pro/) — capabilities and benchmarks
- [Gemini API pricing](https://ai.google.dev/gemini-api/docs/pricing) — token pricing for all models
