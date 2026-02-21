# AI Explorer Dashboard Redesign Plan

## Problem Statement

The current AI Explorer has a **vertically stacked layout** — chat panel on top, data portal below — where visualizations simply stack one after another in a single column. As users ask more questions and accumulate charts, the page becomes a long scroll of disconnected visualizations. It doesn't feel like a **data science dashboard**; it feels like a chat log with charts appended.

**Goal:** Transform the Data Portal into a professional, dashboard-style workspace where multiple visualizations are visible simultaneously in a grid layout, with the chat as a persistent sidebar for directing the analysis.

---

## Current State

```
┌────────────────────────────────────────────┐
│  Header                                     │
├────────────────────────────────────────────┤
│                                             │
│  ┌──────────────────────────────────────┐  │
│  │  Chat Panel (full width)             │  │
│  │  ┌────────────────────────────────┐  │  │
│  │  │  Messages (400px max-height)   │  │  │
│  │  │  [user msg]                    │  │  │
│  │  │  [assistant msg]               │  │  │
│  │  └────────────────────────────────┘  │  │
│  │  [Input box] [Send]                  │  │
│  │  Try: [suggestion] [suggestion] ...  │  │
│  └──────────────────────────────────────┘  │
│                                             │
│  ┌──────────────────────────────────────┐  │
│  │  Data Portal (full width)            │  │
│  │                                      │  │
│  │  ┌──────────────────────────────┐   │  │
│  │  │  Viz 1 (full width)          │   │  │
│  │  └──────────────────────────────┘   │  │
│  │  ┌──────────────────────────────┐   │  │
│  │  │  Viz 2 (full width)          │   │  │
│  │  └──────────────────────────────┘   │  │
│  │  ┌──────────────────────────────┐   │  │
│  │  │  Viz 3 (full width)          │   │  │
│  │  └──────────────────────────────┘   │  │
│  │  ... (keeps stacking forever)       │  │
│  └──────────────────────────────────────┘  │
│                                             │
└────────────────────────────────────────────┘
```

**Problems:**
1. Chat occupies premium top-of-page real estate; users scroll past it to see results
2. Visualizations stack vertically — no overview, no ability to compare charts side by side
3. All viz cards are the same full width regardless of content type (a pie chart doesn't need 1200px)
4. No way to expand a chart for detailed inspection
5. No way to download/export a chart
6. The overall aesthetic is functional but not polished — doesn't convey "professional data science tool"
7. Empty Data Portal placeholder is bland

---

## Target State

```
┌────────────────────────────────────────────────────────────┐
│  Header                                                     │
├────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌──────────┐  ┌──────────────────────────────────────┐   │
│  │ Chat      │  │  Dashboard Canvas                    │   │
│  │ Sidebar   │  │                                      │   │
│  │           │  │  ┌────────────┐  ┌────────────────┐  │   │
│  │ [msgs]    │  │  │ Pie chart  │  │ Gene card      │  │   │
│  │           │  │  │            │  │                │  │   │
│  │           │  │  └────────────┘  └────────────────┘  │   │
│  │           │  │  ┌───────────────────────────────┐   │   │
│  │           │  │  │ Volcano plot (full width)      │   │   │
│  │           │  │  │                               │   │   │
│  │           │  │  └───────────────────────────────┘   │   │
│  │           │  │  ┌───────────────────────────────┐   │   │
│  │           │  │  │ STRING network (full width)    │   │   │
│  │ ───────── │  │  │                               │   │   │
│  │ [input]   │  │  └───────────────────────────────┘   │   │
│  │ [send]    │  │  ┌────────────┐  ┌────────────────┐  │   │
│  │           │  │  │ Bar chart  │  │ Protein struct │  │   │
│  │ Try: ...  │  │  └────────────┘  └────────────────┘  │   │
│  └──────────┘  └──────────────────────────────────────┘   │
│                                                             │
└────────────────────────────────────────────────────────────┘
```

**Key changes:**
1. **Side-by-side layout** — Chat becomes a left sidebar; dashboard canvas fills remaining width
2. **CSS Grid dashboard** — 2-column grid with smart sizing per viz type
3. **Redesigned viz cards** — Type icons, expand/download actions, refined styling
4. **Fullscreen mode** — Any viz can expand to a modal overlay for detailed inspection
5. **Better empty state** — Welcoming illustration with context-aware suggestions
6. **Collapsible chat** — Toggle to maximize dashboard when not actively chatting
7. **Professional polish** — Elevated shadows, refined typography, consistent spacing

---

## Implementation Phases

### Phase 1: Side-by-Side Layout

The highest-impact change. Move from stacked to split layout.

**Files to modify:**
- `fe/pages/chat.py` — Restructure Dash layout from stacked to side-by-side
- `fe/assets/custom.css` — New layout CSS (flexbox/grid for the two-panel split)
- `fe/assets/chat.js` — Adjust scroll behavior, portal targeting (minimal changes)

**Layout approach:**

```css
.ai-explorer-layout {
  display: flex;
  height: calc(100vh - <header-height>);
  gap: 0;
}

.chat-sidebar {
  width: 380px;
  min-width: 320px;
  flex-shrink: 0;
  display: flex;
  flex-direction: column;
  border-right: 1px solid #e2e5e9;
  background: #ffffff;
}

.dashboard-canvas {
  flex: 1;
  overflow-y: auto;
  background: #f0f2f5;  /* subtle contrast from sidebar */
  padding: 1.25rem;
}
```

**Chat sidebar structure:**
- Messages area fills available height (flex-grow: 1, overflow-y: auto)
- Input area fixed at bottom of sidebar
- Suggestions below input
- Collapse/expand toggle button in sidebar header

**Responsive behavior:**
- **>=1200px:** Side-by-side, chat 380px
- **>=992px:** Side-by-side, chat 320px
- **<992px:** Stacked layout (current behavior), chat on top with reduced height
- **<768px:** Full-width stacked, chat as expandable drawer

**Collapse toggle:**
- Small button at the top-right of chat sidebar (chevron icon)
- Collapsed state: chat becomes a thin 48px strip with just the expand button + chat icon
- Dashboard canvas expands to fill the reclaimed space
- Input remains accessible via a floating button in collapsed state

---

### Phase 2: CSS Grid Dashboard

Transform the flat list into a responsive grid.

**Files to modify:**
- `fe/assets/custom.css` — Grid styles for dashboard content
- `fe/assets/chat.js` — Add `data-viz-type` attribute to containers for CSS targeting; add size classes

**Grid approach:**

```css
.dashboard-grid {
  display: grid;
  grid-template-columns: repeat(2, 1fr);
  gap: 1rem;
  align-items: start;  /* cards don't stretch to fill row */
}

/* Size classes applied per viz type */
.viz-container[data-size="small"]  { grid-column: span 1; }  /* half width */
.viz-container[data-size="large"]  { grid-column: span 2; }  /* full width */
```

**Viz type → size mapping:**

| Viz Type | Grid Size | Rationale |
|----------|-----------|-----------|
| `pie_chart` | small (1 col) | Compact, 350px height |
| `bar_chart` | small (1 col) | Compact, doesn't need full width |
| `gene_card` | small (1 col) | Text content, compact |
| `table` | large (2 cols) | Tables need horizontal space |
| `volcano_plot` | large (2 cols) | Scatter plots benefit from width |
| `mave_heatmap` | large (2 cols) | Position axis needs width |
| `protein_structure` | small (1 col) | 3D viewer is aspect-ratio agnostic |
| `gene_interaction_network` | large (2 cols) | Networks need space |
| `string_interaction_network` | large (2 cols) | Networks need space |

**Responsive grid:**
- **>=1200px:** 2-column grid
- **<1200px (stacked layout):** 2-column grid (dashboard is full width)
- **<768px:** 1-column grid

**New viz heights (adapt for smaller cards):**
- Small cards: reduce chart height to ~280px (from 350-450px)
- Large cards: keep 400-450px height
- Protein structure: 380px height

---

### Phase 3: Redesigned Viz Cards

Elevate the visual quality of each dashboard tile.

**Files to modify:**
- `fe/assets/custom.css` — New card styles
- `fe/assets/chat.js` — Update `renderVisualization()` to generate new card markup

**New card anatomy:**

```
┌─────────────────────────────────────┐
│ ┌─┐                                │
│ │📊│ Volcano Plot: TP53       ⛶ ⬇ ✕│
│ └─┘                                │
├─────────────────────────────────────┤
│                                     │
│         [Chart content]             │
│                                     │
│                                     │
└─────────────────────────────────────┘
```

- **Type icon:** Small colored icon for viz category (chart, table, network, protein, card)
- **Title:** Clean, medium-weight typography
- **Action buttons (right side of header):**
  - Expand/fullscreen (⛶ `bi-arrows-fullscreen`)
  - Download image (⬇ `bi-download`) — only for Plotly/Cytoscape charts
  - Remove (✕ `bi-x-lg`)
- **Card body:** White background, subtle border, refined shadow

**Card styling:**

```css
.viz-container {
  background: #ffffff;
  border: 1px solid #e2e5e9;
  border-radius: 10px;
  overflow: hidden;
  box-shadow: 0 1px 3px rgba(0,0,0,0.04), 0 1px 2px rgba(0,0,0,0.06);
  transition: box-shadow 0.2s ease;
}

.viz-container:hover {
  box-shadow: 0 4px 12px rgba(0,0,0,0.08), 0 1px 3px rgba(0,0,0,0.06);
}

.viz-header {
  display: flex;
  align-items: center;
  padding: 0.625rem 0.875rem;
  border-bottom: 1px solid #eef0f3;
  gap: 0.5rem;
}

.viz-type-icon {
  width: 28px;
  height: 28px;
  border-radius: 6px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 0.85rem;
  flex-shrink: 0;
}

/* Color-coded by category */
.viz-type-icon--chart   { background: #e8f5e9; color: #007B53; }
.viz-type-icon--table   { background: #e3f2fd; color: #193F90; }
.viz-type-icon--network { background: #f3e5f5; color: #563D82; }
.viz-type-icon--protein { background: #fff3e0; color: #D4A843; }
.viz-type-icon--card    { background: #fce4ec; color: #A6093D; }

.viz-actions {
  display: flex;
  gap: 0.25rem;
  margin-left: auto;
}

.viz-action-btn {
  width: 28px;
  height: 28px;
  border: none;
  border-radius: 6px;
  background: transparent;
  color: #8c939a;
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: pointer;
  transition: all 0.15s;
  font-size: 0.85rem;
}

.viz-action-btn:hover {
  background: #f0f2f5;
  color: #495057;
}

.viz-action-btn--remove:hover {
  background: #fef2f2;
  color: #dc3545;
}
```

---

### Phase 4: Fullscreen Expand Mode

Allow any visualization to be expanded to a large modal overlay for detailed inspection.

**Files to modify:**
- `fe/assets/chat.js` — Fullscreen logic
- `fe/assets/custom.css` — Modal overlay styles

**Approach:**
- Clicking the expand button clones the viz content into a fullscreen overlay
- For Plotly charts: call `Plotly.relayout()` in the cloned container to resize to full dimensions
- For Cytoscape networks: call `cy.resize()` and `cy.fit()` in the cloned container
- For protein structures: simply resize the container (pdbe-molstar handles it)
- Overlay has close button (✕), title bar, and optional download button
- Escape key or click-outside dismisses

**Overlay structure:**

```css
.viz-fullscreen-overlay {
  position: fixed;
  inset: 0;
  background: rgba(0, 0, 0, 0.5);
  z-index: 1050;
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 2rem;
}

.viz-fullscreen-content {
  background: #ffffff;
  border-radius: 12px;
  width: 100%;
  max-width: 1400px;
  max-height: 90vh;
  overflow: auto;
  box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
}
```

**Implementation detail:** Rather than cloning DOM nodes (which breaks Plotly/Cytoscape bindings), the expand approach should:
1. Move the original viz content div into the fullscreen container
2. Trigger resize on Plotly/Cytoscape
3. On close, move the content div back to its original card
4. Trigger resize again to fit back

---

### Phase 5: Chat Sidebar Polish

Refine the chat UI to feel premium in its new sidebar role.

**Files to modify:**
- `fe/pages/chat.py` — Updated Dash layout
- `fe/assets/custom.css` — Chat sidebar styles
- `fe/assets/chat.js` — Collapse toggle, updated suggestion behavior

**Changes:**

1. **Sidebar header:**
   ```
   ┌─────────────────────────┐
   │ 💬 AI Explorer    [◀]   │
   └─────────────────────────┘
   ```
   Clean header with title and collapse toggle.

2. **Messages area:**
   - Remove max-height; let it fill the sidebar via flex-grow
   - Slightly smaller font (0.9rem) to fit the narrower width
   - User messages: keep green bubbles, reduce max-width to 90%
   - Assistant messages: lighter background, slightly rounder
   - Add subtle date/time separators between conversation turns

3. **Input area:**
   - Pin to bottom of sidebar
   - Textarea instead of text input (for multi-line questions)
   - Send button overlaid at bottom-right of textarea (like modern chat apps)
   - Subtle border-top separator

4. **Suggestions:**
   - Reduce to 3-4 most relevant suggestions
   - Show in collapsed pills below the input
   - When chat is empty, show suggestions more prominently in the messages area as "cards" rather than small pills
   - Optionally: update suggestions based on what's already in the dashboard (context-aware)

5. **Collapsed state:**
   - 48px wide strip with chat bubble icon
   - Floating chat input appears as a bottom-right overlay (like Intercom) when collapsed
   - Badge showing unread assistant message count

---

### Phase 6: Empty State & Dashboard Header

**Files to modify:**
- `fe/pages/chat.py` — Empty state markup
- `fe/assets/custom.css` — Empty state and header styles
- `fe/assets/chat.js` — Show/hide logic

**Dashboard header bar (always visible above the grid):**

```
┌──────────────────────────────────────────────────┐
│  Dashboard                   3 panels  │ Clear all│
└──────────────────────────────────────────────────┘
```

- Title: "Dashboard"
- Item count badge: "3 panels"
- Clear all button (outline, requires confirmation on click)
- Possibly: layout toggle (grid/list) as a nice-to-have

**Empty state (when no visualizations):**

```
┌──────────────────────────────────────────────────┐
│                                                   │
│           ┌───────────────────┐                   │
│           │                   │                   │
│           │   [Chart icon     │                   │
│           │    illustration]  │                   │
│           │                   │                   │
│           └───────────────────┘                   │
│                                                   │
│      Your dashboard will build here               │
│                                                   │
│  Ask a question in the chat to start exploring    │
│  perturbation data. Charts, tables, networks,     │
│  and protein structures will appear as tiles.     │
│                                                   │
│     [Try: "Show TP53 data"]  [Try: "BRCA2"]      │
│                                                   │
└──────────────────────────────────────────────────┘
```

- Centered layout with subtle icon/illustration
- Brief instructional text
- Quick-start buttons that populate the chat input

---

### Phase 7: Download/Export (Enhancement)

**Files to modify:**
- `fe/assets/chat.js` — Download handlers per viz type

**Download capabilities per viz type:**
- **Plotly charts** (pie, bar, volcano, heatmap): Use `Plotly.downloadImage(div, {format: 'png', width: 1200, height: 800})` — Plotly has this built in
- **Cytoscape networks**: Use `cy.png({full: true, scale: 2})` to get a PNG data URL, then trigger download
- **Tables**: Convert to CSV and trigger download
- **Gene cards**: No download (text-based, user can copy)
- **Protein structure**: Link to AlphaFold download page (CIF file)

---

## Design System Notes

### Colors (Existing EBI Brand Palette)
- Primary green: `#007B53`
- Deep blue: `#193F90`
- Crimson: `#A6093D`
- Purple: `#563D82`
- Medium blue: `#3B6FB6`
- Slate: `#54585A`
- Dark green: `#0A5032`
- Gold: `#D4A843`

### New Neutral Tones
- Dashboard background: `#f0f2f5` (warmer than current `#f8f9fa`)
- Card background: `#ffffff`
- Card border: `#e2e5e9`
- Subtle hover shadow: `0 4px 12px rgba(0,0,0,0.08)`
- Text primary: `#1a1d21`
- Text secondary: `#6b7280`

### Typography
- Keep existing system font stack: `-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif`
- Dashboard title: 16px, font-weight 600
- Viz card title: 14px, font-weight 500
- Chat messages: 14px, line-height 1.5
- Small labels: 12px, font-weight 500, uppercase tracking

### Spacing Scale
- xs: 4px
- sm: 8px
- md: 12px
- lg: 16px
- xl: 24px

---

## Priority Order

| # | Phase | Impact | Effort | Dependencies |
|---|-------|--------|--------|-------------|
| 1 | Side-by-side layout | **Very High** | Medium | None |
| 2 | CSS Grid dashboard | **Very High** | Low-Medium | Phase 1 |
| 3 | Redesigned viz cards | **High** | Medium | Phase 2 |
| 4 | Fullscreen expand | **High** | Medium | Phase 3 |
| 5 | Chat sidebar polish | **Medium** | Medium | Phase 1 |
| 6 | Empty state & dashboard header | **Medium** | Low | Phase 2 |
| 7 | Download/export | **Medium** | Low | Phase 3 |

**Phases 1+2 together deliver ~70% of the visual transformation.** They should be implemented together as the first PR.

---

## What We're NOT Doing (and Why)

- **Drag-and-drop reordering (gridstack.js):** Adds a ~30KB dependency and significant JS complexity. The CSS Grid approach gives us the dashboard look without the overhead. Can be added later if users request it.
- **Dashboard saving/persistence:** Would require backend session storage changes. Not in scope for the UI redesign.
- **Dashboard tabs/workspaces:** Interesting for power users but adds significant UI complexity. Revisit after initial redesign is validated.
- **Dark mode:** Nice-to-have but orthogonal to the layout/grid changes. Can be layered on later with CSS custom properties.
- **Real-time dashboard refresh:** Charts are already interactive (Plotly hover, Cytoscape pan/zoom). No need for polling or refresh.

---

## Technical Risks & Mitigations

1. **Plotly resize in grid layout:** Plotly charts may not auto-resize when their container changes width (e.g., sidebar collapse). **Mitigation:** Call `Plotly.Plots.resize(div)` on container resize events using a ResizeObserver.

2. **Cytoscape resize in grid:** Similar issue. **Mitigation:** Call `cy.resize(); cy.fit()` on container resize.

3. **PDbe-molstar in smaller containers:** The 3D viewer needs minimum dimensions. **Mitigation:** Set min-height 350px and test in the smaller (1-col) card size.

4. **Dash layout constraints:** Dash expects a specific layout structure. Moving to a side-by-side layout within `dbc.Container` may need a wrapper div that breaks out of the container's max-width. **Mitigation:** Use `fluid=True` on the container for the chat page, or use a custom wrapper outside the container.

5. **Mobile chat drawer:** Complex touch interaction. **Mitigation:** Keep simple stacked layout on mobile (current behavior) rather than building a drawer. The sidebar collapses to stacked below 992px.

---

## Files to Touch (Summary)

| File | Changes |
|------|---------|
| `fe/pages/chat.py` | Full layout restructure (sidebar + canvas) |
| `fe/assets/custom.css` | Major: new layout, grid, card, fullscreen, empty state styles |
| `fe/assets/chat.js` | Medium: card markup generation, fullscreen logic, resize observers, collapse toggle, download handlers |
| `fe/app.py` | Minor: possibly adjust container for chat page (fluid layout) |
