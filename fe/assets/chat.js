(function () {
  "use strict";

  // Project color palette
  var COLORS = ["#007B53", "#193F90", "#A6093D", "#563D82", "#3B6FB6", "#54585A", "#0A5032", "#D4A843"];

  var sessionId = null;
  var currentController = null;

  function init() {
    var input = document.getElementById("chat-input");
    if (!input) return;

    // Prevent double init
    if (input.dataset.chatInit) return;
    input.dataset.chatInit = "1";

    var sendBtn = document.getElementById("chat-send-btn");
    var clearBtn = document.getElementById("chat-clear-portal-btn");

    input.addEventListener("keydown", function (e) {
      if (e.key === "Enter" && !e.shiftKey) {
        e.preventDefault();
        sendCurrentMessage();
      }
    });

    // Auto-grow textarea height
    input.addEventListener("input", function () {
      this.style.height = "auto";
      var newHeight = Math.min(this.scrollHeight, 120);
      this.style.height = newHeight + "px";
      this.style.overflowY = newHeight >= 120 ? "auto" : "hidden";
    });

    sendBtn.addEventListener("click", function () {
      sendCurrentMessage();
    });

    if (clearBtn) {
      clearBtn.addEventListener("click", function () {
        var portal = document.getElementById("chat-data-portal");
        if (portal) portal.innerHTML = "";
        updatePortalVisibility();
      });
    }

    // Suggestion card buttons (welcome section)
    var suggestionCards = document.querySelectorAll(".chat-suggestion-card");
    suggestionCards.forEach(function (btn) {
      btn.addEventListener("click", function () {
        var query = btn.getAttribute("data-query") || btn.textContent;
        input.value = query;
        sendCurrentMessage();
      });
    });

    // Sidebar collapse / expand
    var layout = document.getElementById("ai-explorer-layout");
    var collapseBtn = document.getElementById("chat-collapse-btn");
    var expandBtn = document.getElementById("chat-expand-btn");

    if (collapseBtn) {
      collapseBtn.addEventListener("click", function () {
        if (layout) layout.classList.add("sidebar-collapsed");
        if (expandBtn) expandBtn.style.display = "";
        resizeAllViz();
      });
    }

    if (expandBtn) {
      expandBtn.addEventListener("click", function () {
        if (layout) layout.classList.remove("sidebar-collapsed");
        expandBtn.style.display = "none";
        input.focus();
        resetUnreadBadge();
        resizeAllViz();
      });
    }

    // Collapsed strip expand button
    var stripBtn = document.getElementById("chat-strip-expand-btn");
    if (stripBtn) {
      stripBtn.addEventListener("click", function () {
        if (layout) layout.classList.remove("sidebar-collapsed");
        if (expandBtn) expandBtn.style.display = "none";
        input.focus();
        resetUnreadBadge();
        resizeAllViz();
      });
    }
  }

  function sendCurrentMessage() {
    var input = document.getElementById("chat-input");
    var text = (input.value || "").trim();
    if (!text) return;
    // Use native setter + event to sync React/Dash controlled state.
    var nativeSetter = Object.getOwnPropertyDescriptor(
      window.HTMLTextAreaElement.prototype, "value"
    ).set;
    nativeSetter.call(input, "");
    input.style.height = "auto";
    input.dispatchEvent(new Event("input", { bubbles: true }));
    hideSuggestions();
    sendMessage(text);
  }

  function hideSuggestions() {
    var welcome = document.getElementById("chat-welcome");
    if (welcome) welcome.style.display = "none";
  }

  function sendMessage(text) {
    appendUserMessage(text);
    showSpinner("Understanding your question...");
    setInputEnabled(false);

    // Abort any in-flight request
    if (currentController) {
      currentController.abort();
    }
    currentController = new AbortController();

    var urlEl = document.getElementById("chat-backend-url");
    var baseUrl = urlEl ? urlEl.getAttribute("data-url") : "";
    var url = baseUrl + "/v1/chat/stream";

    var body = JSON.stringify({ message: text, session_id: sessionId });
    var assistantStarted = false;

    fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: body,
      signal: currentController.signal
    })
      .then(function (response) {
        if (!response.ok) {
          throw new Error("HTTP " + response.status);
        }
        var reader = response.body.getReader();
        var decoder = new TextDecoder();
        var buffer = "";

        function read() {
          return reader.read().then(function (result) {
            if (result.done) {
              hideSpinner();
              if (assistantStarted) finalizeAssistantMessage();
              setInputEnabled(true);
              currentController = null;
              return;
            }

            buffer += decoder.decode(result.value, { stream: true });
            var lines = buffer.split("\n");
            buffer = lines.pop(); // keep incomplete line

            var currentEvent = null;
            for (var i = 0; i < lines.length; i++) {
              var line = lines[i];
              if (line.startsWith("event: ")) {
                currentEvent = line.substring(7).trim();
              } else if (line.startsWith("data: ") && currentEvent) {
                try {
                  var data = JSON.parse(line.substring(6));
                  handleSSEEvent(currentEvent, data);
                  if (currentEvent === "text" && !assistantStarted) {
                    assistantStarted = true;
                  }
                } catch (e) {
                  // ignore parse errors
                }
                currentEvent = null;
              }
            }

            return read();
          });
        }

        return read();
      })
      .catch(function (err) {
        if (err.name === "AbortError") return;
        hideSpinner();
        if (assistantStarted) finalizeAssistantMessage();
        appendErrorMessage("Connection error: " + err.message);
        setInputEnabled(true);
        currentController = null;
      });
  }

  function handleSSEEvent(event, data) {
    switch (event) {
      case "thinking":
        showSpinner(data.status || "Thinking...");
        break;
      case "tool_call":
        showSpinner(data.description || "Calling tool...");
        break;
      case "text":
        hideSpinner();
        appendAssistantText(data.content || "");
        break;
      case "visualization":
        renderVisualization(data);
        break;
      case "error":
        hideSpinner();
        appendErrorMessage(data.message || "An error occurred");
        break;
      case "done":
        if (data.session_id) sessionId = data.session_id;
        hideSpinner();
        break;
    }
  }

  // --- Spinner with cycling brightness icons ---

  var spinnerIconInterval = null;
  var SPINNER_ICONS = [
    "bi bi-flower1 spinner-icon",
    "bi bi-flower2 spinner-icon",
    "bi bi-flower3 spinner-icon"
  ];
  var spinnerIconIndex = 0;

  function showSpinner(status) {
    var el = document.getElementById("chat-spinner");
    if (!el) return;
    el.style.display = "flex";

    var statusEl = el.querySelector(".spinner-status");
    if (statusEl && status) statusEl.textContent = status;

    // Start icon cycling
    if (!spinnerIconInterval) {
      spinnerIconIndex = 0;
      var iconEl = el.querySelector(".spinner-icon");
      if (iconEl) {
        spinnerIconInterval = setInterval(function () {
          spinnerIconIndex = (spinnerIconIndex + 1) % SPINNER_ICONS.length;
          iconEl.className = SPINNER_ICONS[spinnerIconIndex];
        }, 400);
      }
    }
  }

  function hideSpinner() {
    var el = document.getElementById("chat-spinner");
    if (el) el.style.display = "none";

    if (spinnerIconInterval) {
      clearInterval(spinnerIconInterval);
      spinnerIconInterval = null;
    }
  }

  // --- Chat messages ---

  var currentAssistantBubble = null;

  function appendUserMessage(text) {
    var container = document.getElementById("chat-messages");
    if (!container) return;
    // Add turn separator if there are existing messages
    if (container.querySelector(".chat-message")) {
      var sep = document.createElement("div");
      sep.className = "chat-turn-separator";
      var now = new Date();
      sep.textContent = now.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
      container.appendChild(sep);
    }
    var msg = document.createElement("div");
    msg.className = "chat-message chat-message-user";
    msg.textContent = text;
    container.appendChild(msg);
    scrollToBottom(container);
  }

  function appendAssistantText(chunk) {
    var container = document.getElementById("chat-messages");
    if (!container) return;

    if (!currentAssistantBubble) {
      currentAssistantBubble = document.createElement("div");
      currentAssistantBubble.className = "chat-message chat-message-assistant";
      container.appendChild(currentAssistantBubble);
      incrementUnread();
    }
    currentAssistantBubble.innerHTML += formatMarkdown(chunk);
    scrollToBottom(container);
  }

  function finalizeAssistantMessage() {
    currentAssistantBubble = null;
  }

  function appendErrorMessage(text) {
    var container = document.getElementById("chat-messages");
    if (!container) return;
    var msg = document.createElement("div");
    msg.className = "chat-message chat-message-error";
    msg.textContent = text;
    container.appendChild(msg);
    scrollToBottom(container);
  }

  function scrollToBottom(el) {
    el.scrollTop = el.scrollHeight;
  }

  // --- Unread badge tracking ---

  var unreadCount = 0;

  function incrementUnread() {
    var layout = document.getElementById("ai-explorer-layout");
    if (!layout || !layout.classList.contains("sidebar-collapsed")) return;
    unreadCount++;
    var badge = document.getElementById("chat-unread-badge");
    if (badge) {
      badge.textContent = unreadCount;
      badge.classList.add("visible");
    }
  }

  function resetUnreadBadge() {
    unreadCount = 0;
    var badge = document.getElementById("chat-unread-badge");
    if (badge) {
      badge.textContent = "";
      badge.classList.remove("visible");
    }
  }

  function setInputEnabled(enabled) {
    var input = document.getElementById("chat-input");
    var btn = document.getElementById("chat-send-btn");
    if (input) input.disabled = !enabled;
    if (btn) btn.disabled = !enabled;
    if (enabled && input) input.focus();
  }

  // --- Simple markdown formatting ---

  function formatMarkdown(text) {
    // Bullet lists (before bold/italic so "* " at line start isn't treated as emphasis)
    text = text.replace(/(^|\n)\* /g, "$1\u2022 ");
    text = text.replace(/(^|\n)- /g, "$1\u2022 ");
    // Bold
    text = text.replace(/\*\*(.*?)\*\*/g, "<strong>$1</strong>");
    // Italic
    text = text.replace(/\*(.*?)\*/g, "<em>$1</em>");
    // Inline code
    text = text.replace(/`(.*?)`/g, "<code>$1</code>");
    // Line breaks
    text = text.replace(/\n/g, "<br>");
    return text;
  }

  // --- Visualizations ---

  // Viz type → grid size mapping (small = 1 col, large = 2 cols)
  var VIZ_SIZE_MAP = {
    pie_chart: "small",
    bar_chart: "small",
    gene_card: "small",
    protein_structure: "small",
    table: "large",
    volcano_plot: "large",
    mave_heatmap: "large",
    gene_interaction_network: "large",
    string_interaction_network: "large"
  };

  // Viz type → category for icon coloring
  var VIZ_CATEGORY_MAP = {
    pie_chart: "chart",
    bar_chart: "chart",
    volcano_plot: "chart",
    mave_heatmap: "chart",
    table: "table",
    gene_interaction_network: "network",
    string_interaction_network: "network",
    protein_structure: "protein",
    gene_card: "card"
  };

  // Category → Bootstrap icon class
  var VIZ_ICON_MAP = {
    chart: "bi bi-bar-chart-line",
    table: "bi bi-table",
    network: "bi bi-diagram-3",
    protein: "bi bi-box",
    card: "bi bi-card-text"
  };

  // Viz types that support image download (and their renderer type)
  var DOWNLOADABLE_TYPES = {
    pie_chart: "plotly",
    bar_chart: "plotly",
    volcano_plot: "plotly",
    mave_heatmap: "plotly",
    gene_interaction_network: "cytoscape",
    string_interaction_network: "cytoscape"
  };

  function updatePortalVisibility() {
    var portal = document.getElementById("chat-data-portal");
    var clearBtn = document.getElementById("chat-clear-portal-btn");
    var ph = document.getElementById("chat-portal-placeholder");
    var countEl = document.getElementById("chat-panel-count");
    var count = portal ? portal.children.length : 0;
    var hasItems = count > 0;
    if (clearBtn) clearBtn.style.display = hasItems ? "" : "none";
    if (ph) ph.style.display = hasItems ? "none" : "";
    if (countEl) countEl.textContent = hasItems ? count + (count === 1 ? " panel" : " panels") : "";
  }

  // Resize Plotly charts and Cytoscape instances after layout changes
  function resizeAllViz() {
    // Small delay to let CSS transition finish
    setTimeout(function () {
      var portal = document.getElementById("chat-data-portal");
      if (!portal) return;

      // Resize all Plotly charts
      if (typeof Plotly !== "undefined") {
        var plotDivs = portal.querySelectorAll(".js-plotly-plot");
        plotDivs.forEach(function (div) {
          Plotly.Plots.resize(div);
        });
      }

      // Resize all Cytoscape instances
      portal.querySelectorAll(".viz-container").forEach(function (container) {
        if (container._cyInstance) {
          container._cyInstance.resize();
          container._cyInstance.fit();
        }
      });
    }, 300);
  }

  function renderVisualization(data) {
    var portal = document.getElementById("chat-data-portal");
    if (!portal) return;

    // Hide placeholder
    var ph = document.getElementById("chat-portal-placeholder");
    if (ph) ph.style.display = "none";

    // Show clear button
    var clearBtn = document.getElementById("chat-clear-portal-btn");
    if (clearBtn) clearBtn.style.display = "";

    var vizType = data.type || "";
    var category = VIZ_CATEGORY_MAP[vizType] || "chart";

    var wrapper = document.createElement("div");
    wrapper.className = "viz-container";
    wrapper.setAttribute("data-viz-type", vizType);
    wrapper.setAttribute("data-size", VIZ_SIZE_MAP[vizType] || "large");

    // ── Card header with type icon, title, and action buttons ──
    var header = document.createElement("div");
    header.className = "viz-header";

    // Type icon (color-coded by category)
    var typeIcon = document.createElement("div");
    typeIcon.className = "viz-type-icon viz-type-icon--" + category;
    var iconEl = document.createElement("i");
    iconEl.className = VIZ_ICON_MAP[category] || "bi bi-bar-chart-line";
    typeIcon.appendChild(iconEl);
    header.appendChild(typeIcon);

    // Title
    var title = document.createElement("span");
    title.className = "viz-title";
    title.textContent = data.title || "Visualization";
    header.appendChild(title);

    // Action buttons group
    var actions = document.createElement("div");
    actions.className = "viz-actions";

    // Expand button
    var fullscreenBtn = document.createElement("button");
    fullscreenBtn.className = "viz-action-btn";
    fullscreenBtn.title = "Expand";
    fullscreenBtn.innerHTML = '<i class="bi bi-arrows-fullscreen"></i>';
    fullscreenBtn.addEventListener("click", function () {
      openFullscreen(wrapper, data.title || "Visualization");
    });
    actions.appendChild(fullscreenBtn);

    // Download button (only for Plotly/Cytoscape viz types)
    if (DOWNLOADABLE_TYPES[vizType]) {
      var dlBtn = document.createElement("button");
      dlBtn.className = "viz-action-btn";
      dlBtn.title = "Download image";
      dlBtn.innerHTML = '<i class="bi bi-download"></i>';
      dlBtn.addEventListener("click", function () {
        downloadViz(wrapper, vizType, data.title || "visualization");
      });
      actions.appendChild(dlBtn);
    }

    // Remove button
    var removeBtn = document.createElement("button");
    removeBtn.className = "viz-action-btn viz-action-btn--remove";
    removeBtn.title = "Remove";
    removeBtn.innerHTML = '<i class="bi bi-x-lg"></i>';
    removeBtn.addEventListener("click", function () {
      wrapper.remove();
      updatePortalVisibility();
    });
    actions.appendChild(removeBtn);

    header.appendChild(actions);
    wrapper.appendChild(header);

    var content = document.createElement("div");
    content.className = "viz-content";
    wrapper.appendChild(content);

    portal.appendChild(wrapper);

    switch (data.type) {
      case "table":
        renderTable(content, data.data);
        break;
      case "pie_chart":
        renderPieChart(content, data.data);
        break;
      case "bar_chart":
        renderBarChart(content, data.data);
        break;
      case "volcano_plot":
        renderVolcanoPlot(content, data.data);
        break;
      case "mave_heatmap":
        renderMaveHeatmap(content, data.data);
        break;
      case "protein_structure":
        renderProteinStructure(content, data.data);
        break;
      case "gene_card":
        renderGeneCard(content, data.data);
        break;
      case "gene_interaction_network":
        renderGeneInteractionNetwork(content, data.data);
        break;
      case "string_interaction_network":
        renderStringInteractionNetwork(content, data.data);
        break;
      default:
        content.textContent = "Unknown visualization type: " + data.type;
    }

    // Scroll portal into view
    wrapper.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  function renderTable(container, data) {
    var headers = data.headers || [];
    var rows = data.rows || [];

    var table = document.createElement("table");
    table.className = "table table-sm table-striped table-hover viz-table";

    // Header
    var thead = document.createElement("thead");
    var headerRow = document.createElement("tr");
    headers.forEach(function (h) {
      var th = document.createElement("th");
      th.textContent = h;
      headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);
    table.appendChild(thead);

    // Body
    var tbody = document.createElement("tbody");
    rows.forEach(function (row) {
      var tr = document.createElement("tr");
      row.forEach(function (cell) {
        var td = document.createElement("td");
        td.textContent = cell != null ? String(cell) : "";
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
    table.appendChild(tbody);

    var wrapper = document.createElement("div");
    wrapper.className = "table-responsive";
    wrapper.appendChild(table);
    container.appendChild(wrapper);
  }

  function renderPieChart(container, data) {
    var labels = data.labels || [];
    var values = data.values || [];

    var chartDiv = document.createElement("div");
    chartDiv.style.width = "100%";
    container.appendChild(chartDiv);

    if (typeof Plotly === "undefined") {
      chartDiv.textContent = "Chart library not loaded";
      return;
    }

    var chartColors = labels.map(function (_, i) {
      return COLORS[i % COLORS.length];
    });

    Plotly.newPlot(
      chartDiv,
      [{
        type: "pie",
        labels: labels,
        values: values,
        marker: { colors: chartColors },
        textinfo: "label+percent",
        hoverinfo: "label+value+percent",
        hole: 0.35
      }],
      {
        margin: { t: 10, b: 10, l: 10, r: 10 },
        showlegend: true,
        legend: { orientation: "h", y: -0.1 },
        font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif" },
        height: 280
      },
      { responsive: true, displayModeBar: false }
    );
  }

  function renderBarChart(container, data) {
    var labels = data.labels || [];
    var values = data.values || [];
    var xlabel = data.xlabel || "";
    var ylabel = data.ylabel || "";

    var chartDiv = document.createElement("div");
    chartDiv.style.width = "100%";
    container.appendChild(chartDiv);

    if (typeof Plotly === "undefined") {
      chartDiv.textContent = "Chart library not loaded";
      return;
    }

    var chartColors = labels.map(function (_, i) {
      return COLORS[i % COLORS.length];
    });

    Plotly.newPlot(
      chartDiv,
      [{
        type: "bar",
        x: labels,
        y: values,
        marker: { color: chartColors }
      }],
      {
        margin: { t: 10, b: 60, l: 60, r: 20 },
        xaxis: { title: xlabel, tickangle: labels.length > 6 ? -45 : 0 },
        yaxis: { title: ylabel },
        font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif" },
        height: 280
      },
      { responsive: true, displayModeBar: false }
    );
  }

  function renderVolcanoPlot(container, data) {
    var rawGenes = data.genes || [];
    var rawLog2fc = data.log2fc || [];
    var rawPadj = data.padj || [];
    var fcThreshold = data.fc_threshold != null ? parseFloat(data.fc_threshold) : 1.0;
    var padjThreshold = data.padj_threshold != null ? parseFloat(data.padj_threshold) : 0.05;
    var perturbedGene = data.perturbed_gene || "";

    var chartDiv = document.createElement("div");
    chartDiv.style.width = "100%";
    container.appendChild(chartDiv);

    if (typeof Plotly === "undefined") {
      chartDiv.textContent = "Chart library not loaded";
      return;
    }

    // Coerce to numbers and filter out invalid rows
    var genes = [], log2fc = [], padj = [];
    var len = Math.min(rawGenes.length, rawLog2fc.length, rawPadj.length);
    for (var i = 0; i < len; i++) {
      var fc = parseFloat(rawLog2fc[i]);
      var pv = parseFloat(rawPadj[i]);
      if (isFinite(fc) && isFinite(pv) && pv >= 0) {
        genes.push(String(rawGenes[i] || ""));
        log2fc.push(fc);
        padj.push(pv);
      }
    }

    if (genes.length === 0) {
      chartDiv.textContent = "No valid data points for volcano plot";
      return;
    }

    // Transform padj to -log10(padj), clamping minimum padj to 1e-300
    var negLog10Padj = padj.map(function (p) {
      if (p <= 0) return 300;
      return -Math.log10(Math.max(p, 1e-300));
    });
    var negLog10Threshold = -Math.log10(padjThreshold);

    // Classify points: up (significant + positive FC), down (significant + negative FC), ns (not significant)
    var upIdx = [], downIdx = [], nsIdx = [];
    for (var j = 0; j < genes.length; j++) {
      var sig = negLog10Padj[j] >= negLog10Threshold;
      if (sig && log2fc[j] > fcThreshold) {
        upIdx.push(j);
      } else if (sig && log2fc[j] < -fcThreshold) {
        downIdx.push(j);
      } else {
        nsIdx.push(j);
      }
    }

    function subset(arr, indices) {
      return indices.map(function (k) { return arr[k]; });
    }

    function makeHoverText(geneArr, fcArr, pArr) {
      return geneArr.map(function (g, k) {
        return g + "<br>log2FC: " + fcArr[k].toFixed(3) + "<br>padj: " + pArr[k].toExponential(2);
      });
    }

    function makeTrace(name, indices, color, opacity) {
      var tGenes = subset(genes, indices);
      var tFc = subset(log2fc, indices);
      var tPadj = subset(padj, indices);
      return {
        type: "scatter",
        mode: "markers",
        name: name + " (" + indices.length + ")",
        x: tFc,
        y: subset(negLog10Padj, indices),
        text: makeHoverText(tGenes, tFc, tPadj),
        hoverinfo: "text",
        marker: { size: 5, color: color, opacity: opacity }
      };
    }

    var traces = [];
    if (nsIdx.length > 0) traces.push(makeTrace("Not significant", nsIdx, "#B0B0B0", 0.5));
    if (upIdx.length > 0) traces.push(makeTrace("Up", upIdx, "#A6093D", 0.7));
    if (downIdx.length > 0) traces.push(makeTrace("Down", downIdx, "#193F90", 0.7));

    // Threshold lines
    var fcAbsMax = 0;
    for (var m = 0; m < log2fc.length; m++) {
      var absVal = Math.abs(log2fc[m]);
      if (absVal > fcAbsMax) fcAbsMax = absVal;
    }
    var xMax = Math.max(fcAbsMax, (fcThreshold || 1) + 0.5);

    var yMax = 0;
    for (var n = 0; n < negLog10Padj.length; n++) {
      if (negLog10Padj[n] > yMax) yMax = negLog10Padj[n];
    }

    // Horizontal padj threshold line (always shown)
    var shapes = [
      { type: "line", x0: -xMax - 0.5, x1: xMax + 0.5, y0: negLog10Threshold, y1: negLog10Threshold, line: { color: "#54585A", width: 1, dash: "dash" } }
    ];
    // Vertical FC threshold lines (only if fc_threshold > 0)
    if (fcThreshold > 0) {
      shapes.push({ type: "line", x0: fcThreshold, x1: fcThreshold, y0: 0, y1: yMax * 1.05, line: { color: "#54585A", width: 1, dash: "dash" } });
      shapes.push({ type: "line", x0: -fcThreshold, x1: -fcThreshold, y0: 0, y1: yMax * 1.05, line: { color: "#54585A", width: 1, dash: "dash" } });
    }

    var subtitle = perturbedGene ? "Perturbation: " + perturbedGene : "";

    Plotly.newPlot(
      chartDiv,
      traces,
      {
        margin: { t: subtitle ? 30 : 10, b: 60, l: 60, r: 20 },
        xaxis: { title: "log<sub>2</sub> Fold Change", zeroline: true, zerolinecolor: "#ddd" },
        yaxis: { title: "-log<sub>10</sub>(p<sub>adj</sub>)", rangemode: "tozero" },
        shapes: shapes,
        font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif" },
        legend: { orientation: "h", y: -0.2, x: 0.5, xanchor: "center" },
        annotations: subtitle ? [{
          text: subtitle,
          xref: "paper", yref: "paper",
          x: 0.5, y: 1.02,
          xanchor: "center", yanchor: "bottom",
          showarrow: false,
          font: { size: 12, color: "#54585A" }
        }] : [],
        height: 450,
        hovermode: "closest"
      },
      { responsive: true, displayModeBar: false }
    );
  }

  function renderMaveHeatmap(container, data) {
    var z = data.z || [];
    var aminoAcids = data.amino_acids || [];
    var positions = data.positions || [];
    var wtAnnotations = data.wt_annotations || [];
    var geneName = data.gene_name || "";
    var posStart = data.position_start;
    var posEnd = data.position_end;

    var chartDiv = document.createElement("div");
    chartDiv.style.width = "100%";
    container.appendChild(chartDiv);

    if (typeof Plotly === "undefined") {
      chartDiv.textContent = "Chart library not loaded";
      return;
    }

    if (z.length === 0 || positions.length === 0 || aminoAcids.length === 0) {
      chartDiv.textContent = "No valid data for MAVE heatmap";
      return;
    }

    // Build annotation list for WT residues
    var annotations = [];
    for (var i = 0; i < wtAnnotations.length; i++) {
      var wt = wtAnnotations[i];
      annotations.push({
        x: positions[wt.col],
        y: aminoAcids[wt.row],
        text: "WT",
        showarrow: false,
        font: { color: "#000", size: 9 }
      });
    }

    var subtitle = "";
    if (posStart != null && posEnd != null) {
      subtitle = "Positions " + posStart + "–" + posEnd;
    }

    var trace = {
      type: "heatmap",
      z: z,
      x: positions,
      y: aminoAcids,
      colorscale: "RdYlGn",
      hoverongaps: false,
      hovertemplate: "Position: %{x}<br>AA: %{y}<br>Score: %{z:.3f}<extra></extra>",
      colorbar: {
        title: { text: "Score", side: "right" },
        thickness: 15,
        len: 0.9
      }
    };

    Plotly.newPlot(
      chartDiv,
      [trace],
      {
        margin: { t: subtitle ? 30 : 10, b: 60, l: 50, r: 80 },
        xaxis: {
          title: "Position",
          tickmode: "linear",
          dtick: positions.length > 40 ? 5 : 1,
          tickangle: positions.length > 20 ? -45 : 0
        },
        yaxis: {
          title: "Amino Acid",
          tickmode: "linear",
          autorange: "reversed"
        },
        font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif" },
        annotations: subtitle
          ? annotations.concat([{
              text: subtitle,
              xref: "paper", yref: "paper",
              x: 0.5, y: 1.02,
              xanchor: "center", yanchor: "bottom",
              showarrow: false,
              font: { size: 12, color: "#54585A" }
            }])
          : annotations,
        height: 450,
        hovermode: "closest"
      },
      { responsive: true, displayModeBar: false }
    );
  }

  function renderProteinStructure(container, data) {
    var viewerDiv = document.createElement("div");
    viewerDiv.style.width = "100%";
    viewerDiv.style.height = "380px";
    viewerDiv.style.position = "relative";
    viewerDiv.style.overflow = "hidden";
    container.appendChild(viewerDiv);

    var entryId = data.entry_id || "";
    var uniprotId = data.uniprot_id || "";

    // Use pdbe-molstar web component
    var viewer = document.createElement("pdbe-molstar");
    viewer.setAttribute("alphafold-view", "true");
    viewer.setAttribute("hide-water", "true");
    viewer.setAttribute("hide-controls", "true");
    viewer.setAttribute("bg-color-r", "255");
    viewer.setAttribute("bg-color-g", "255");
    viewer.setAttribute("bg-color-b", "255");
    viewer.style.width = "100%";
    viewer.style.height = "100%";
    viewer.style.display = "block";

    // Load from AlphaFold using custom-data URL for the CIF file
    var cifUrl = data.cif_url || "";
    if (cifUrl) {
      viewer.setAttribute("custom-data-url", cifUrl);
      viewer.setAttribute("custom-data-format", "cif");
    } else if (entryId) {
      // Fallback: construct the AlphaFold CIF URL from entry ID
      viewer.setAttribute("custom-data-url",
        "https://alphafold.ebi.ac.uk/files/" + entryId + "-model_v4.cif");
      viewer.setAttribute("custom-data-format", "cif");
    }

    viewerDiv.appendChild(viewer);

    // Info row below viewer
    var info = document.createElement("div");
    info.style.cssText = "padding: 8px 4px; font-size: 0.85rem; color: #54585A;";
    var parts = [];
    if (data.gene) parts.push("<strong>" + escapeHtml(data.gene) + "</strong>");
    if (data.organism) parts.push(escapeHtml(data.organism));
    if (uniprotId) {
      parts.push('<a href="https://www.uniprot.org/uniprot/' + encodeURIComponent(uniprotId) +
        '" target="_blank" rel="noopener">UniProt: ' + escapeHtml(uniprotId) + '</a>');
    }
    if (entryId) {
      parts.push('<a href="https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(uniprotId || entryId) +
        '" target="_blank" rel="noopener">AlphaFold</a>');
    }
    info.innerHTML = parts.join(" &middot; ");
    container.appendChild(info);
  }

  function renderGeneCard(container, data) {
    var card = document.createElement("div");
    card.className = "gene-card";
    card.style.cssText = "font-size: 0.9rem;";

    // Header
    var header = document.createElement("div");
    header.style.cssText = "margin-bottom: 10px;";
    header.innerHTML = '<div style="font-size: 1.1rem; font-weight: 600; color: #193F90;">' +
      escapeHtml(data.gene_name || "") + '</div>' +
      '<div style="color: #54585A; font-size: 0.85rem;">' +
      escapeHtml(data.protein_name || "") + '</div>';
    card.appendChild(header);

    // Function
    if (data.function) {
      var funcDiv = document.createElement("div");
      funcDiv.style.cssText = "margin-bottom: 10px;";
      var funcText = data.function;
      var truncated = funcText.length > 200;
      var funcId = "gene-card-func-" + Date.now();
      funcDiv.innerHTML = '<div style="font-weight: 600; margin-bottom: 2px;">Function</div>' +
        '<div id="' + funcId + '">' +
        escapeHtml(truncated ? funcText.substring(0, 200) + "..." : funcText) +
        (truncated ? ' <a href="#" style="color: #007B53;" onclick="' +
          "this.parentElement.textContent='" + escapeAttr(funcText) + "'; return false;" +
          '">show more</a>' : '') +
        '</div>';
      card.appendChild(funcDiv);
    }

    // Key facts grid
    var facts = [];
    if (data.subcellular_location) {
      facts.push({label: "Location", value: data.subcellular_location});
    }
    if (data.domains && data.domains.length) {
      facts.push({label: "Domains", value: data.domains.join(", ")});
    }
    if (data.diseases && data.diseases.length) {
      facts.push({label: "Disease associations", value: data.diseases.join(", ")});
    }
    if (data.go_terms && data.go_terms.length) {
      facts.push({label: "GO terms", value: data.go_terms.slice(0, 8).join(", ")});
    }

    if (facts.length) {
      var factsDiv = document.createElement("div");
      factsDiv.style.cssText = "margin-bottom: 10px;";
      facts.forEach(function (fact) {
        var row = document.createElement("div");
        row.style.cssText = "margin-bottom: 4px;";
        row.innerHTML = '<span style="font-weight: 600;">' + escapeHtml(fact.label) + ':</span> ' +
          escapeHtml(fact.value);
        factsDiv.appendChild(row);
      });
      card.appendChild(factsDiv);
    }

    // External links
    var links = [];
    if (data.uniprot_id) {
      links.push('<a href="https://www.uniprot.org/uniprot/' + encodeURIComponent(data.uniprot_id) +
        '" target="_blank" rel="noopener" style="color: #007B53; text-decoration: none; margin-right: 12px;">' +
        '<i class="bi bi-box-arrow-up-right"></i> UniProt</a>');
    }
    if (data.alphafold_id || data.uniprot_id) {
      var afId = data.uniprot_id || data.alphafold_id;
      links.push('<a href="https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(afId) +
        '" target="_blank" rel="noopener" style="color: #007B53; text-decoration: none; margin-right: 12px;">' +
        '<i class="bi bi-box-arrow-up-right"></i> AlphaFold</a>');
    }
    if (data.gene_name) {
      links.push('<a href="https://platform.opentargets.org/search?q=' + encodeURIComponent(data.gene_name) +
        '&page=1&entities=target" target="_blank" rel="noopener" style="color: #007B53; text-decoration: none;">' +
        '<i class="bi bi-box-arrow-up-right"></i> Open Targets</a>');
    }

    if (links.length) {
      var linksDiv = document.createElement("div");
      linksDiv.style.cssText = "padding-top: 8px; border-top: 1px solid #e9ecef;";
      linksDiv.innerHTML = links.join("");
      card.appendChild(linksDiv);
    }

    container.appendChild(card);
  }

  function renderGeneInteractionNetwork(container, data) {
    var nodes = data.nodes || [];
    var edges = data.edges || [];
    var perturbedGene = data.perturbed_gene || "";

    var cyDiv = document.createElement("div");
    cyDiv.style.width = "100%";
    cyDiv.style.height = "420px";
    cyDiv.style.border = "1px solid #e9ecef";
    cyDiv.style.borderRadius = "6px";
    cyDiv.style.background = "#ffffff";
    container.appendChild(cyDiv);

    if (typeof cytoscape === "undefined") {
      cyDiv.textContent = "Cytoscape.js library not loaded";
      return;
    }

    if (nodes.length === 0) {
      cyDiv.textContent = "No data for gene interaction network";
      return;
    }

    // Build Cytoscape elements
    var cyNodes = nodes.map(function (n) {
      return {
        data: {
          id: n.id,
          label: n.label,
          type: n.type,
          log2fc: n.log2fc,
          padj: n.padj
        }
      };
    });

    var cyEdges = edges.map(function (e, idx) {
      return {
        data: {
          id: "e" + idx,
          source: e.source,
          target: e.target,
          weight: e.weight
        }
      };
    });

    // Compute max weight for edge width scaling
    var maxWeight = 0;
    for (var i = 0; i < edges.length; i++) {
      if (edges[i].weight > maxWeight) maxWeight = edges[i].weight;
    }
    if (maxWeight === 0) maxWeight = 1;

    var cy = cytoscape({
      container: cyDiv,
      elements: cyNodes.concat(cyEdges),
      style: [
        {
          selector: "node[type='perturbed']",

          style: {
            "background-color": "#193F90",
            "label": "data(label)",
            "text-valign": "center",
            "text-halign": "center",
            "color": "#fff",
            "font-size": "12px",
            "font-weight": "bold",
            "width": 50,
            "height": 50,
            "border-width": 3,
            "border-color": "#0d2456",
            "text-outline-width": 0
          }
        },
        {
          selector: "node[type='up']",
          style: {
            "background-color": "#A6093D",
            "label": "data(label)",
            "text-valign": "bottom",
            "text-margin-y": 4,
            "color": "#333",
            "font-size": "10px",
            "width": "mapData(log2fc, 0, " + maxWeight + ", 20, 40)",
            "height": "mapData(log2fc, 0, " + maxWeight + ", 20, 40)"
          }
        },
        {
          selector: "node[type='down']",
          style: {
            "background-color": "#193F90",
            "label": "data(label)",
            "text-valign": "bottom",
            "text-margin-y": 4,
            "color": "#333",
            "font-size": "10px",
            "width": function (ele) {
              return Math.max(20, Math.min(40, 20 + (Math.abs(ele.data("log2fc")) / maxWeight) * 20));
            },
            "height": function (ele) {
              return Math.max(20, Math.min(40, 20 + (Math.abs(ele.data("log2fc")) / maxWeight) * 20));
            }
          }
        },
        {
          selector: "edge",
          style: {
            "width": function (ele) {
              return Math.max(1, (ele.data("weight") / maxWeight) * 6);
            },
            "line-color": "#ccc",
            "target-arrow-color": "#ccc",
            "curve-style": "bezier",
            "opacity": 0.6
          }
        },
        {
          selector: "node:active",
          style: {
            "overlay-opacity": 0
          }
        }
      ],
      layout: {
        name: "concentric",
        concentric: function (node) {
          return node.data("type") === "perturbed" ? 10 : 1;
        },
        levelWidth: function () { return 1; },
        minNodeSpacing: 30,
        padding: 20,
        animate: false
      },
      userZoomingEnabled: true,
      userPanningEnabled: true,
      boxSelectionEnabled: false
    });

    // Store cy instance for resize handling
    container.closest(".viz-container")._cyInstance = cy;

    // Tooltip on tap
    var tooltip = document.createElement("div");
    tooltip.style.cssText =
      "position: absolute; background: rgba(0,0,0,0.85); color: #fff; padding: 8px 12px; " +
      "border-radius: 6px; font-size: 12px; pointer-events: none; display: none; z-index: 10; " +
      "max-width: 220px; line-height: 1.4;";
    cyDiv.style.position = "relative";
    cyDiv.appendChild(tooltip);

    cy.on("mouseover", "node", function (evt) {
      var node = evt.target;
      var d = node.data();
      var lines = ["<strong>" + escapeHtml(d.label) + "</strong>"];
      if (d.type === "perturbed") {
        lines.push("Perturbed gene (center)");
      } else {
        lines.push("log2FC: " + d.log2fc.toFixed(3));
        lines.push("padj: " + (d.padj < 0.001 ? d.padj.toExponential(2) : d.padj.toFixed(4)));
        lines.push(d.type === "up" ? "Upregulated" : "Downregulated");
      }
      tooltip.innerHTML = lines.join("<br>");
      tooltip.style.display = "block";
      var pos = node.renderedPosition();
      tooltip.style.left = (pos.x + 15) + "px";
      tooltip.style.top = (pos.y - 10) + "px";
    });

    cy.on("mouseout", "node", function () {
      tooltip.style.display = "none";
    });

    cy.on("pan zoom", function () {
      tooltip.style.display = "none";
    });

    // Legend
    var legend = document.createElement("div");
    legend.style.cssText = "padding: 8px 4px; font-size: 0.8rem; color: #54585A; display: flex; gap: 16px; align-items: center; flex-wrap: wrap;";
    legend.innerHTML =
      '<span style="display:inline-flex; align-items:center; gap:4px;">' +
        '<span style="display:inline-block; width:12px; height:12px; border-radius:50%; background:#193F90; border: 2px solid #0d2456;"></span> ' +
        escapeHtml(perturbedGene) + ' (perturbed)</span>' +
      '<span style="display:inline-flex; align-items:center; gap:4px;">' +
        '<span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:#A6093D;"></span> Upregulated</span>' +
      '<span style="display:inline-flex; align-items:center; gap:4px;">' +
        '<span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:#193F90;"></span> Downregulated</span>' +
      '<span style="display:inline-flex; align-items:center; gap:4px;">' +
        '<span style="display:inline-block; width:20px; height:3px; background:#999;"></span> Edge width = |log2FC|</span>';
    container.appendChild(legend);
  }

  function renderStringInteractionNetwork(container, data) {
    var nodes = data.nodes || [];
    var edges = data.edges || [];
    var queryGene = data.query_gene || "";
    var requiredScore = data.required_score || 400;

    var cyDiv = document.createElement("div");
    cyDiv.style.width = "100%";
    cyDiv.style.height = "420px";
    cyDiv.style.border = "1px solid #e9ecef";
    cyDiv.style.borderRadius = "6px";
    cyDiv.style.background = "#fafafa";
    container.appendChild(cyDiv);

    if (typeof cytoscape === "undefined") {
      cyDiv.textContent = "Cytoscape.js library not loaded";
      return;
    }

    if (nodes.length === 0) {
      cyDiv.textContent = "No data for STRING interaction network";
      return;
    }

    var COLOR_QUERY = "#007B53";
    var COLOR_PARTNER = "#563D82";
    var COLOR_BORDER = "#004d34";

    var cyNodes = nodes.map(function (n) {
      return {
        data: {
          id: n.id,
          label: n.label,
          type: n.type,
          score: n.score || 0,
          degree: n.degree || 0,
          evidence: n.evidence || {}
        }
      };
    });

    var cyEdges = edges.map(function (e, idx) {
      return {
        data: {
          id: "se" + idx,
          source: e.source,
          target: e.target,
          weight: e.weight || 0
        }
      };
    });

    var cy = cytoscape({
      container: cyDiv,
      elements: cyNodes.concat(cyEdges),
      style: [
        {
          selector: "node[type='query']",
          style: {
            "background-color": COLOR_QUERY,
            "label": "data(label)",
            "text-valign": "center",
            "text-halign": "center",
            "color": "#fff",
            "font-size": "12px",
            "font-weight": "bold",
            "width": 54,
            "height": 54,
            "border-width": 3,
            "border-color": COLOR_BORDER,
            "text-outline-width": 0
          }
        },
        {
          selector: "node[type='partner']",
          style: {
            "background-color": COLOR_PARTNER,
            "label": "data(label)",
            "text-valign": "bottom",
            "text-margin-y": 4,
            "color": "#333",
            "font-size": "10px",
            "width": function (ele) {
              return 22 + ele.data("score") * 18;
            },
            "height": function (ele) {
              return 22 + ele.data("score") * 18;
            },
            "border-width": 1,
            "border-color": "#3d2a60"
          }
        },
        {
          selector: "edge",
          style: {
            "width": function (ele) {
              return Math.max(1, ele.data("weight") * 6);
            },
            "line-color": "#9b89c0",
            "curve-style": "bezier",
            "opacity": 0.65
          }
        },
        {
          selector: "node:active",
          style: { "overlay-opacity": 0 }
        }
      ],
      layout: {
        name: "concentric",
        concentric: function (node) {
          return node.data("type") === "query" ? 10 : 1;
        },
        levelWidth: function () { return 1; },
        minNodeSpacing: 28,
        padding: 20,
        animate: false
      },
      userZoomingEnabled: true,
      userPanningEnabled: true,
      boxSelectionEnabled: false
    });

    // Store cy instance for resize handling
    container.closest(".viz-container")._cyInstance = cy;

    // Tooltip
    var tooltip = document.createElement("div");
    tooltip.style.cssText =
      "position: absolute; background: rgba(0,0,0,0.85); color: #fff; padding: 8px 12px; " +
      "border-radius: 6px; font-size: 12px; pointer-events: none; display: none; z-index: 10; " +
      "max-width: 250px; line-height: 1.5;";
    cyDiv.style.position = "relative";
    cyDiv.appendChild(tooltip);

    cy.on("mouseover", "node", function (evt) {
      var node = evt.target;
      var d = node.data();
      var lines = ["<strong>" + escapeHtml(d.label) + "</strong>"];

      if (d.type === "query") {
        lines.push("Query gene (center)");
        lines.push("Partners shown: " + d.degree);
      } else {
        lines.push("Confidence: " + (d.score * 1000).toFixed(0) + " / 1000");
        var ev = d.evidence || {};
        if (ev.experimental > 0) lines.push("Experimental: " + (ev.experimental * 1000).toFixed(0));
        if (ev.database > 0) lines.push("Database: " + (ev.database * 1000).toFixed(0));
        if (ev.coexpression > 0) lines.push("Coexpression: " + (ev.coexpression * 1000).toFixed(0));
        if (ev.textmining > 0) lines.push("Textmining: " + (ev.textmining * 1000).toFixed(0));
      }

      tooltip.innerHTML = lines.join("<br>");
      tooltip.style.display = "block";
      var pos = node.renderedPosition();
      tooltip.style.left = (pos.x + 15) + "px";
      tooltip.style.top = (pos.y - 10) + "px";
    });

    cy.on("mouseout", "node", function () {
      tooltip.style.display = "none";
    });

    cy.on("pan zoom", function () {
      tooltip.style.display = "none";
    });

    // Legend
    var legend = document.createElement("div");
    legend.style.cssText =
      "padding: 8px 4px; font-size: 0.8rem; color: #54585A; " +
      "display: flex; gap: 14px; align-items: center; flex-wrap: wrap;";
    legend.innerHTML =
      '<span style="display:inline-flex; align-items:center; gap:4px;">' +
        '<span style="display:inline-block; width:12px; height:12px; border-radius:50%; background:' + COLOR_QUERY + '; border:2px solid ' + COLOR_BORDER + ';"></span> ' +
        escapeHtml(queryGene) + " (query)</span>" +
      '<span style="display:inline-flex; align-items:center; gap:4px;">' +
        '<span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:' + COLOR_PARTNER + ';"></span> Interaction partner</span>' +
      '<span style="display:inline-flex; align-items:center; gap:4px;">' +
        '<span style="display:inline-block; width:20px; height:3px; background:#9b89c0;"></span> Edge width = confidence</span>' +
      '<span style="color:#888;">Min score: ' + requiredScore + '/1000</span>' +
      '<span style="margin-left:auto; color:#007B53;">' +
        '<a href="https://string-db.org" target="_blank" rel="noopener" ' +
        'style="color:#007B53; text-decoration:none;">Data: STRING-DB</a></span>';
    container.appendChild(legend);
  }

  // --- Download and Fullscreen helpers ---

  function downloadViz(container, vizType, title) {
    var downloadType = DOWNLOADABLE_TYPES[vizType];
    var filename = title.replace(/[^a-z0-9]/gi, "_").toLowerCase();

    if (downloadType === "plotly" && typeof Plotly !== "undefined") {
      var plotDiv = container.querySelector(".js-plotly-plot");
      if (plotDiv) {
        Plotly.downloadImage(plotDiv, {
          format: "png",
          width: 1200,
          height: 800,
          filename: filename
        });
      }
    } else if (downloadType === "cytoscape" && container._cyInstance) {
      var pngData = container._cyInstance.png({ full: true, scale: 2 });
      var link = document.createElement("a");
      link.href = pngData;
      link.download = filename + ".png";
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }
  }

  function openFullscreen(container, title) {
    // Create overlay backdrop
    var overlay = document.createElement("div");
    overlay.className = "viz-fullscreen-overlay";

    var panel = document.createElement("div");
    panel.className = "viz-fullscreen-content";

    // ── Header with icon, title, and action buttons ──
    var header = document.createElement("div");
    header.className = "viz-fullscreen-header";

    // Clone the type icon from the card header
    var origIcon = container.querySelector(".viz-type-icon");
    if (origIcon) {
      header.appendChild(origIcon.cloneNode(true));
    }

    var titleEl = document.createElement("span");
    titleEl.className = "viz-fullscreen-title";
    titleEl.textContent = title;
    header.appendChild(titleEl);

    var actions = document.createElement("div");
    actions.className = "viz-actions";

    // Download button (if viz type supports it)
    var vizType = container.getAttribute("data-viz-type");
    if (DOWNLOADABLE_TYPES[vizType]) {
      var dlBtn = document.createElement("button");
      dlBtn.className = "viz-action-btn";
      dlBtn.title = "Download image";
      dlBtn.innerHTML = '<i class="bi bi-download"></i>';
      dlBtn.addEventListener("click", function () {
        downloadViz(container, vizType, title);
      });
      actions.appendChild(dlBtn);
    }

    // Close button
    var closeBtn = document.createElement("button");
    closeBtn.className = "viz-action-btn";
    closeBtn.title = "Close";
    closeBtn.innerHTML = '<i class="bi bi-x-lg"></i>';
    closeBtn.addEventListener("click", close);
    actions.appendChild(closeBtn);

    header.appendChild(actions);
    panel.appendChild(header);

    // ── Body — move the original viz-content into fullscreen ──
    var body = document.createElement("div");
    body.className = "viz-fullscreen-body";
    panel.appendChild(body);

    var vizContent = container.querySelector(".viz-content");
    // Insert a placeholder so we know where to put it back
    var placeholder = document.createElement("div");
    placeholder.className = "viz-content-placeholder";
    vizContent.parentNode.insertBefore(placeholder, vizContent);
    body.appendChild(vizContent);

    overlay.appendChild(panel);
    document.body.appendChild(overlay);
    document.body.style.overflow = "hidden";

    // Expand chart/network heights to fill the modal
    requestAnimationFrame(function () {
      var bodyRect = body.getBoundingClientRect();
      var expandHeight = Math.max(bodyRect.height - 32, 400);

      // Store and expand inline-height children (cytoscape divs, protein viewer, etc.)
      var expandedEls = [];
      var children = vizContent.children;
      for (var i = 0; i < children.length; i++) {
        var child = children[i];
        if (child.style.height && child.style.height.indexOf("px") !== -1) {
          expandedEls.push({ el: child, orig: child.style.height });
          child.style.height = expandHeight + "px";
        }
      }
      vizContent._expandedEls = expandedEls;

      // Expand Plotly chart heights
      if (typeof Plotly !== "undefined") {
        var plots = vizContent.querySelectorAll(".js-plotly-plot");
        plots.forEach(function (div) {
          div._origHeight = (div.layout || {}).height;
          Plotly.relayout(div, { height: expandHeight });
        });
      }

      // Resize Cytoscape to fit new container
      if (container._cyInstance) {
        container._cyInstance.resize();
        container._cyInstance.fit();
      }
    });

    // Click on backdrop closes
    overlay.addEventListener("click", function (e) {
      if (e.target === overlay) close();
    });

    // Escape key closes
    function onKeyDown(e) {
      if (e.key === "Escape") close();
    }
    document.addEventListener("keydown", onKeyDown);

    function close() {
      // Restore inline-height children
      var expandedEls = vizContent._expandedEls || [];
      expandedEls.forEach(function (item) {
        item.el.style.height = item.orig;
      });
      delete vizContent._expandedEls;

      // Restore Plotly chart heights
      if (typeof Plotly !== "undefined") {
        var plots = vizContent.querySelectorAll(".js-plotly-plot");
        plots.forEach(function (div) {
          if (div._origHeight != null) {
            Plotly.relayout(div, { height: div._origHeight });
            delete div._origHeight;
          }
        });
      }

      // Move viz-content back to its original card
      placeholder.parentNode.insertBefore(vizContent, placeholder);
      placeholder.remove();

      // Remove overlay
      overlay.remove();
      document.removeEventListener("keydown", onKeyDown);
      document.body.style.overflow = "";

      // Trigger resize back to card size
      setTimeout(function () {
        if (typeof Plotly !== "undefined") {
          var plots = vizContent.querySelectorAll(".js-plotly-plot");
          plots.forEach(function (div) {
            Plotly.Plots.resize(div);
          });
        }
        if (container._cyInstance) {
          container._cyInstance.resize();
          container._cyInstance.fit();
        }
      }, 50);
    }
  }

  function escapeHtml(str) {
    var div = document.createElement("div");
    div.appendChild(document.createTextNode(str));
    return div.innerHTML;
  }

  function escapeAttr(str) {
    return str.replace(/\\/g, "\\\\").replace(/'/g, "\\'").replace(/\n/g, " ");
  }

  // --- ResizeObserver to keep Plotly/Cytoscape in sync with grid ---

  var resizeObserverActive = false;

  function setupResizeObserver() {
    if (resizeObserverActive) return;
    var portal = document.getElementById("chat-data-portal");
    if (!portal || typeof ResizeObserver === "undefined") return;

    var ro = new ResizeObserver(function (entries) {
      for (var i = 0; i < entries.length; i++) {
        var target = entries[i].target;
        // Resize Plotly charts inside the resized container
        if (typeof Plotly !== "undefined") {
          var plots = target.querySelectorAll(".js-plotly-plot");
          plots.forEach(function (div) {
            Plotly.Plots.resize(div);
          });
        }
        // Resize Cytoscape instances
        var vizContainer = target.closest(".viz-container") || target;
        if (vizContainer._cyInstance) {
          vizContainer._cyInstance.resize();
          vizContainer._cyInstance.fit();
        }
      }
    });

    // Observe each viz-container's content area for size changes
    var containers = portal.querySelectorAll(".viz-content");
    containers.forEach(function (el) { ro.observe(el); });

    // Also observe new containers as they're added
    var portalObserver = new MutationObserver(function (mutations) {
      mutations.forEach(function (m) {
        m.addedNodes.forEach(function (node) {
          if (node.nodeType === 1) {
            var content = node.querySelector ? node.querySelector(".viz-content") : null;
            if (content) ro.observe(content);
          }
        });
      });
    });
    portalObserver.observe(portal, { childList: true });
    resizeObserverActive = true;
  }

  // --- Initialize when chat page is rendered ---

  // Use MutationObserver to detect when Dash renders the chat page
  var observer = new MutationObserver(function () {
    init();
    setupResizeObserver();
  });

  // Start observing
  if (document.body) {
    observer.observe(document.body, { childList: true, subtree: true });
  }

  // Also try immediately in case page is already rendered
  if (document.readyState === "complete" || document.readyState === "interactive") {
    setTimeout(function () { init(); setupResizeObserver(); }, 100);
  } else {
    document.addEventListener("DOMContentLoaded", function () {
      setTimeout(function () { init(); setupResizeObserver(); }, 100);
    });
  }
})();
