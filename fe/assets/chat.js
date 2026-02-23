(function () {
  "use strict";

  // Project color palette
  var COLORS = ["#007B53", "#193F90", "#A6093D", "#563D82", "#3B6FB6", "#54585A", "#0A5032", "#D4A843"];

  var sessionId = sessionStorage.getItem("chat_session_id") || null;
  var currentController = null;

  function init() {
    var input = document.getElementById("chat-input");
    if (!input) return;

    // Auth check: redirect to login if no token
    if (!sessionStorage.getItem("auth_token")) {
      window.location.href = "/perturbation-catalogue/login";
      return;
    }

    // Prevent double init
    if (input.dataset.chatInit) return;
    input.dataset.chatInit = "1";

    // Show logout link in header when authenticated
    var logoutLink = document.getElementById("logout-link");
    if (logoutLink) logoutLink.style.display = "";

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
      var clearConfirmTimer = null;
      clearBtn.addEventListener("click", function () {
        if (clearBtn.dataset.confirming === "1") {
          // Second click — actually clear
          clearBtn.dataset.confirming = "";
          clearBtn.innerHTML = '<i class="bi bi-trash3 me-1"></i>Clear all';
          clearBtn.classList.remove("dashboard-clear-btn--confirm");
          if (clearConfirmTimer) { clearTimeout(clearConfirmTimer); clearConfirmTimer = null; }
          var portal = document.getElementById("chat-data-portal");
          if (portal) portal.innerHTML = "";
          updatePortalVisibility();
          // Bulk-delete visualizations from DB
          if (sessionId) {
            var urlEl = document.getElementById("chat-backend-url");
            var baseUrl = urlEl ? urlEl.getAttribute("data-url") : "";
            var authToken = sessionStorage.getItem("auth_token");
            fetch(baseUrl + "/v1/chat/sessions/" + sessionId + "/visualizations", {
              method: "DELETE",
              headers: { "Authorization": "Bearer " + authToken }
            });
          }
        } else {
          // First click — ask for confirmation
          clearBtn.dataset.confirming = "1";
          clearBtn.innerHTML = '<i class="bi bi-exclamation-triangle me-1"></i>Confirm?';
          clearBtn.classList.add("dashboard-clear-btn--confirm");
          clearConfirmTimer = setTimeout(function () {
            clearBtn.dataset.confirming = "";
            clearBtn.innerHTML = '<i class="bi bi-trash3 me-1"></i>Clear all';
            clearBtn.classList.remove("dashboard-clear-btn--confirm");
          }, 3000);
        }
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

    // Help guide button
    var helpBtn = document.getElementById("chat-help-btn");
    if (helpBtn) {
      helpBtn.addEventListener("click", function () {
        openHelpGuide();
      });
    }

    // Empty-state "Try:" buttons (dashboard canvas)
    var emptyTryBtns = document.querySelectorAll(".empty-state-try-btn");
    emptyTryBtns.forEach(function (btn) {
      btn.addEventListener("click", function () {
        var query = btn.getAttribute("data-query") || btn.textContent;
        input.value = query;
        sendCurrentMessage();
      });
    });

    // Sidebar collapse / expand
    var layout = document.getElementById("ai-explorer-layout");
    var collapseBtn = document.getElementById("chat-collapse-btn");

    if (collapseBtn) {
      collapseBtn.addEventListener("click", function () {
        if (layout) layout.classList.add("sidebar-collapsed");
        resizeAllViz();
      });
    }

    // Collapsed strip expand button
    var stripBtn = document.getElementById("chat-strip-expand-btn");
    if (stripBtn) {
      stripBtn.addEventListener("click", function () {
        if (layout) layout.classList.remove("sidebar-collapsed");
        input.focus();
        resetUnreadBadge();
        resizeAllViz();
      });
    }

    // "New chat" button
    var newSessionBtn = document.getElementById("chat-new-session-btn");
    if (newSessionBtn) {
      newSessionBtn.addEventListener("click", function () {
        sessionId = null;
        sessionStorage.removeItem("chat_session_id");
        clearChat();
        clearDashboard();
        showWelcome();
        loadSessionList();
      });
    }

    // Restore session on page load
    if (sessionId) {
      loadSession(sessionId);
    }
    loadSessionList();
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

    var authToken = sessionStorage.getItem("auth_token");
    var headers = { "Content-Type": "application/json" };
    if (authToken) {
      headers["Authorization"] = "Bearer " + authToken;
    }

    fetch(url, {
      method: "POST",
      headers: headers,
      body: body,
      signal: currentController.signal
    })
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
        if (data.session_id) {
          sessionId = data.session_id;
          sessionStorage.setItem("chat_session_id", data.session_id);
          loadSessionList();
        }
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

  // Viz types that support download (and their renderer type)
  var DOWNLOADABLE_TYPES = {
    pie_chart: "plotly",
    bar_chart: "plotly",
    volcano_plot: "plotly",
    mave_heatmap: "plotly",
    gene_interaction_network: "cytoscape",
    string_interaction_network: "cytoscape",
    table: "csv",
    protein_structure: "protein"
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
    if (data.viz_id) wrapper.setAttribute("data-viz-id", data.viz_id);

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

    // Download button
    if (DOWNLOADABLE_TYPES[vizType]) {
      var dlBtn = document.createElement("button");
      dlBtn.className = "viz-action-btn";
      var dlType = DOWNLOADABLE_TYPES[vizType];
      dlBtn.title = dlType === "csv" ? "Download CSV" :
                    dlType === "protein" ? "Open in AlphaFold" : "Download image";
      dlBtn.innerHTML = dlType === "protein"
        ? '<i class="bi bi-box-arrow-up-right"></i>'
        : '<i class="bi bi-download"></i>';
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
      var vizId = wrapper.getAttribute("data-viz-id");
      wrapper.remove();
      updatePortalVisibility();
      if (vizId && sessionId) {
        var urlEl = document.getElementById("chat-backend-url");
        var baseUrl = urlEl ? urlEl.getAttribute("data-url") : "";
        var authToken = sessionStorage.getItem("auth_token");
        fetch(baseUrl + "/v1/chat/sessions/" + sessionId + "/visualizations/" + vizId, {
          method: "DELETE",
          headers: { "Authorization": "Bearer " + authToken }
        });
      }
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
    } else if (downloadType === "csv") {
      var table = container.querySelector("table");
      if (!table) return;
      var csv = [];
      var rows = table.querySelectorAll("tr");
      for (var i = 0; i < rows.length; i++) {
        var cells = rows[i].querySelectorAll("th, td");
        var row = [];
        for (var j = 0; j < cells.length; j++) {
          var text = cells[j].textContent.replace(/"/g, '""');
          row.push('"' + text + '"');
        }
        csv.push(row.join(","));
      }
      var blob = new Blob([csv.join("\n")], { type: "text/csv;charset=utf-8;" });
      var url = URL.createObjectURL(blob);
      var link = document.createElement("a");
      link.href = url;
      link.download = filename + ".csv";
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      URL.revokeObjectURL(url);
    } else if (downloadType === "protein") {
      // Open AlphaFold entry page where the user can download the CIF file
      var afLink = container.querySelector('a[href*="alphafold.ebi.ac.uk/entry"]');
      if (afLink) {
        window.open(afLink.href, "_blank", "noopener");
      }
    }
  }

  // ── Help Guide Modal ──
  function openHelpGuide() {
    // Prevent duplicate overlays
    if (document.querySelector(".help-guide-overlay")) return;

    var overlay = document.createElement("div");
    overlay.className = "help-guide-overlay";

    var panel = document.createElement("div");
    panel.className = "help-guide-panel";

    // Header
    var header = document.createElement("div");
    header.className = "help-guide-header";
    header.innerHTML =
      '<div class="help-guide-header-left">' +
        '<i class="bi bi-book" style="font-size:1.1rem;color:#007B53"></i>' +
        '<span class="help-guide-title">AI Explorer Guide</span>' +
      '</div>' +
      '<button class="viz-action-btn help-guide-close" title="Close"><i class="bi bi-x-lg"></i></button>';

    panel.appendChild(header);

    // Body
    var body = document.createElement("div");
    body.className = "help-guide-body";

    body.innerHTML =
      // ── Intro ──
      '<section class="help-section">' +
        '<p class="help-intro">' +
          'AI Explorer is your conversational gateway to the Perturbation Catalogue and a suite of biomedical databases. ' +
          'Ask questions in plain English — the AI will query the right services and render interactive visualizations on your dashboard.' +
        '</p>' +
      '</section>' +

      // ── Data Sources ──
      '<section class="help-section">' +
        '<h3 class="help-section-title"><i class="bi bi-database"></i> Integrated Data Sources</h3>' +
        '<div class="help-sources-grid">' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#e8f5e9;color:#007B53"><i class="bi bi-clipboard2-data"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>Perturbation Catalogue</strong>' +
              '<span>Internal database of Perturb-seq, CRISPR screen, and MAVE experiments. Search datasets, query raw results, and get catalogue statistics.</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#e3f2fd;color:#193F90"><i class="bi bi-body-text"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>UniProt</strong>' +
              '<span>Protein function, domains, disease associations, subcellular location, GO terms, known variants, and identifier mapping.</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#fff3e0;color:#D4A843"><i class="bi bi-box"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>AlphaFold</strong>' +
              '<span>AI-predicted 3D protein structures from the EBI AlphaFold database, rendered as interactive molecular viewers.</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#fce4ec;color:#A6093D"><i class="bi bi-capsule"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>Pharos</strong>' +
              '<span>NIH druggability classification (Tclin/Tchem/Tbio/Tdark), existing drugs, ligands, and target development level.</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#f3e5f5;color:#563D82"><i class="bi bi-diagram-3"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>Open Targets</strong>' +
              '<span>Disease associations, drug mechanisms, GWAS evidence, and genetic evidence linking targets to diseases.</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#e8f5e9;color:#0A5032"><i class="bi bi-share"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>STRING</strong>' +
              '<span>Protein-protein interaction networks with confidence scores and evidence channels. Functional enrichment analysis (GO, KEGG, Reactome).</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#e3f2fd;color:#3B6FB6"><i class="bi bi-journal-text"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>Europe PMC</strong>' +
              '<span>Biomedical literature search across PubMed and open-access full text with text-mined entity annotations.</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#fff3e0;color:#D4A843"><i class="bi bi-gear"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>ProtVar</strong>' +
              '<span>Variant molecular consequences: protein stability (FoldX), pathogenicity scores (EVE, ESM-1b, AlphaMissense, CADD), structural context.</span>' +
            '</div>' +
          '</div>' +

          '<div class="help-source-card">' +
            '<div class="help-source-icon" style="background:#fce4ec;color:#A6093D"><i class="bi bi-lightning"></i></div>' +
            '<div class="help-source-info">' +
              '<strong>Ensembl VEP</strong>' +
              '<span>Variant consequence prediction with SpliceAI, LOFTEE, regulatory impact, SIFT, PolyPhen, and batch annotation for up to 200 variants.</span>' +
            '</div>' +
          '</div>' +

        '</div>' +
      '</section>' +

      // ── What You Can Ask ──
      '<section class="help-section">' +
        '<h3 class="help-section-title"><i class="bi bi-chat-left-text"></i> What You Can Ask</h3>' +

        '<div class="help-category">' +
          '<h4 class="help-category-title">Search & explore the catalogue</h4>' +
          '<p class="help-category-desc">Find datasets, query experimental results, and get catalogue statistics.</p>' +
          '<div class="help-examples">' +
            '<button class="help-example-btn" data-query="What perturbation data is available for TP53 across all modalities?">What data is available for TP53?</button>' +
            '<button class="help-example-btn" data-query="Find CRISPR screen datasets studying lung cancer">CRISPR screen datasets for lung cancer</button>' +
            '<button class="help-example-btn" data-query="How many datasets are in the Perturbation Catalogue?">How many datasets are in the catalogue?</button>' +
            '<button class="help-example-btn" data-query="Show me the top differentially expressed genes when TP53 is knocked out in Perturb-seq">Top DEGs when TP53 is knocked out</button>' +
          '</div>' +
        '</div>' +

        '<div class="help-category">' +
          '<h4 class="help-category-title">Protein information & structure</h4>' +
          '<p class="help-category-desc">Look up protein function, domains, and view predicted 3D structures.</p>' +
          '<div class="help-examples">' +
            '<button class="help-example-btn" data-query="Look up the BRCA2 protein and show me a gene card with its function, domains, and disease associations">Look up the BRCA2 protein</button>' +
            '<button class="help-example-btn" data-query="Show me the predicted 3D protein structure for KRAS">Show KRAS 3D structure</button>' +
            '<button class="help-example-btn" data-query="What known variants does TP53 have in UniProt?">Known variants for TP53</button>' +
            '<button class="help-example-btn" data-query="Map BRCA1, BRCA2, TP53 to their UniProt accession IDs">Map gene names to UniProt IDs</button>' +
          '</div>' +
        '</div>' +

        '<div class="help-category">' +
          '<h4 class="help-category-title">Variant analysis</h4>' +
          '<p class="help-category-desc">Annotate variants with pathogenicity scores, stability predictions, and consequence types.</p>' +
          '<div class="help-examples">' +
            '<button class="help-example-btn" data-query="Annotate the variant rs28897696 using ProtVar. What is its pathogenicity and protein stability impact?">Annotate variant rs28897696</button>' +
            '<button class="help-example-btn" data-query="What is the structural context of position 248 in TP53 (UniProt P04637)?">Structural context of TP53 position 248</button>' +
            '<button class="help-example-btn" data-query="Use Ensembl VEP to predict the consequence of rs56116432">VEP consequence for rs56116432</button>' +
            '<button class="help-example-btn" data-query="Batch annotate these variants: rs1042779, rs28897696, rs56116432">Batch annotate multiple variants</button>' +
          '</div>' +
        '</div>' +

        '<div class="help-category">' +
          '<h4 class="help-category-title">Druggability & disease links</h4>' +
          '<p class="help-category-desc">Check if a target is druggable and explore disease associations via Pharos and Open Targets.</p>' +
          '<div class="help-examples">' +
            '<button class="help-example-btn" data-query="Use Pharos to check if EGFR is druggable. What is its target development level?">Is EGFR druggable?</button>' +
            '<button class="help-example-btn" data-query="What diseases are associated with BRCA2 according to Open Targets?">Diseases associated with BRCA2</button>' +
            '<button class="help-example-btn" data-query="What drugs target KRAS and what is their mechanism of action?">Drugs targeting KRAS</button>' +
          '</div>' +
        '</div>' +

        '<div class="help-category">' +
          '<h4 class="help-category-title">Protein interactions & pathways</h4>' +
          '<p class="help-category-desc">Explore protein interaction networks and run functional enrichment analysis.</p>' +
          '<div class="help-examples">' +
            '<button class="help-example-btn" data-query="Show me the STRING protein-protein interaction network for TP53">STRING network for TP53</button>' +
            '<button class="help-example-btn" data-query="Run functional enrichment analysis on TP53, MDM2, CDKN2A, RB1, ATM">Enrichment for a gene set</button>' +
            '<button class="help-example-btn" data-query="Show me a gene interaction network for TP53 based on Perturb-seq differential expression">DEG interaction network for TP53</button>' +
          '</div>' +
        '</div>' +

        '<div class="help-category">' +
          '<h4 class="help-category-title">Literature search</h4>' +
          '<p class="help-category-desc">Search biomedical publications across PubMed and Europe PMC.</p>' +
          '<div class="help-examples">' +
            '<button class="help-example-btn" data-query="Find recent papers about CRISPR screens in acute myeloid leukemia">Papers on CRISPR screens in AML</button>' +
            '<button class="help-example-btn" data-query="Search for publications about BRCA2 perturbation experiments">BRCA2 perturbation papers</button>' +
          '</div>' +
        '</div>' +

        '<div class="help-category">' +
          '<h4 class="help-category-title">Visualizations</h4>' +
          '<p class="help-category-desc">Request specific chart types for your data.</p>' +
          '<div class="help-examples">' +
            '<button class="help-example-btn" data-query="Show a volcano plot for TP53 perturbation in Perturb-seq">Volcano plot for TP53</button>' +
            '<button class="help-example-btn" data-query="Show a MAVE heatmap for BRCA1 variant effects">MAVE heatmap for BRCA1</button>' +
          '</div>' +
        '</div>' +

      '</section>' +

      // ── Dashboard Tiles ──
      '<section class="help-section">' +
        '<h3 class="help-section-title"><i class="bi bi-grid-1x2"></i> Dashboard Tile Types</h3>' +
        '<div class="help-viz-list">' +
          '<div class="help-viz-item"><span class="help-viz-badge" style="background:#e8f5e9;color:#007B53"><i class="bi bi-pie-chart"></i></span> <strong>Pie &amp; Bar Charts</strong> — Distribution breakdowns of categorical data</div>' +
          '<div class="help-viz-item"><span class="help-viz-badge" style="background:#e8f5e9;color:#007B53"><i class="bi bi-table"></i></span> <strong>Tables</strong> — Structured data rows with sorting</div>' +
          '<div class="help-viz-item"><span class="help-viz-badge" style="background:#e8f5e9;color:#007B53"><i class="bi bi-activity"></i></span> <strong>Volcano Plots</strong> — Differential expression (log2FC vs. significance)</div>' +
          '<div class="help-viz-item"><span class="help-viz-badge" style="background:#e8f5e9;color:#007B53"><i class="bi bi-grid-3x3"></i></span> <strong>MAVE Heatmaps</strong> — Variant effect maps (position &times; amino acid)</div>' +
          '<div class="help-viz-item"><span class="help-viz-badge" style="background:#f3e5f5;color:#563D82"><i class="bi bi-diagram-3"></i></span> <strong>Interaction Networks</strong> — Gene/protein interaction graphs (Cytoscape.js)</div>' +
          '<div class="help-viz-item"><span class="help-viz-badge" style="background:#fff3e0;color:#D4A843"><i class="bi bi-box"></i></span> <strong>3D Protein Structures</strong> — Interactive AlphaFold viewer (PDBe-Molstar)</div>' +
          '<div class="help-viz-item"><span class="help-viz-badge" style="background:#fce4ec;color:#A6093D"><i class="bi bi-card-text"></i></span> <strong>Gene Cards</strong> — Summary cards with function, domains, disease links</div>' +
        '</div>' +
        '<p class="help-viz-note">Every tile can be expanded to full screen, and most can be downloaded as PNG or CSV using the buttons in the tile header.</p>' +
      '</section>' +

      // ── Tips ──
      '<section class="help-section">' +
        '<h3 class="help-section-title"><i class="bi bi-lightbulb"></i> Tips</h3>' +
        '<ul class="help-tips">' +
          '<li>Be specific with gene names — use official HGNC symbols (e.g. <em>TP53</em>, <em>BRCA2</em>, <em>KRAS</em>).</li>' +
          '<li>For variants, you can use rsIDs (<em>rs28897696</em>), HGVS notation, or genomic coordinates.</li>' +
          '<li>Ask follow-up questions — the AI remembers your conversation context within a session.</li>' +
          '<li>Request specific visualizations by name: <em>"Show a volcano plot"</em>, <em>"Create a heatmap"</em>.</li>' +
          '<li>Combine queries: <em>"Look up TP53, show its structure, and check if it\'s druggable."</em></li>' +
          '<li>For batch operations, list multiple items: <em>"Annotate variants rs123, rs456, rs789."</em></li>' +
        '</ul>' +
      '</section>';

    panel.appendChild(body);
    overlay.appendChild(panel);
    document.body.appendChild(overlay);

    // Wire up example buttons to populate chat input
    var exampleBtns = overlay.querySelectorAll(".help-example-btn");
    exampleBtns.forEach(function (btn) {
      btn.addEventListener("click", function () {
        var input = document.getElementById("chat-input");
        if (input) {
          input.value = btn.getAttribute("data-query") || btn.textContent;
          input.focus();
        }
        closeGuide();
      });
    });

    // Close logic
    function closeGuide() {
      overlay.remove();
    }

    var closeBtn = overlay.querySelector(".help-guide-close");
    if (closeBtn) {
      closeBtn.addEventListener("click", closeGuide);
    }

    overlay.addEventListener("click", function (e) {
      if (e.target === overlay) closeGuide();
    });

    document.addEventListener("keydown", function onEsc(e) {
      if (e.key === "Escape") {
        closeGuide();
        document.removeEventListener("keydown", onEsc);
      }
    });
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
      var dlType = DOWNLOADABLE_TYPES[vizType];
      var dlBtn = document.createElement("button");
      dlBtn.className = "viz-action-btn";
      dlBtn.title = dlType === "csv" ? "Download CSV" :
                    dlType === "protein" ? "Open in AlphaFold" : "Download image";
      dlBtn.innerHTML = dlType === "protein"
        ? '<i class="bi bi-box-arrow-up-right"></i>'
        : '<i class="bi bi-download"></i>';
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

  // --- Session management ---

  function getBaseUrl() {
    var urlEl = document.getElementById("chat-backend-url");
    return urlEl ? urlEl.getAttribute("data-url") : "";
  }

  function getAuthHeaders() {
    var authToken = sessionStorage.getItem("auth_token");
    return authToken ? { "Authorization": "Bearer " + authToken } : {};
  }

  function clearChat() {
    var container = document.getElementById("chat-messages");
    if (!container) return;
    // Remove all messages and separators, keep welcome section
    var children = container.querySelectorAll(".chat-message, .chat-turn-separator");
    children.forEach(function (el) { el.remove(); });
    currentAssistantBubble = null;
  }

  function clearDashboard() {
    var portal = document.getElementById("chat-data-portal");
    if (portal) portal.innerHTML = "";
    updatePortalVisibility();
  }

  function showWelcome() {
    var welcome = document.getElementById("chat-welcome");
    if (welcome) welcome.style.display = "";
  }

  function loadSession(id) {
    var baseUrl = getBaseUrl();
    fetch(baseUrl + "/v1/chat/sessions/" + id, {
      headers: getAuthHeaders()
    })
    .then(function (r) {
      if (!r.ok) {
        sessionId = null;
        sessionStorage.removeItem("chat_session_id");
        return null;
      }
      return r.json();
    })
    .then(function (session) {
      if (!session) return;

      hideSuggestions();
      session.messages.forEach(function (msg) {
        if (msg.role === "user") {
          appendUserMessage(msg.content);
        } else {
          appendAssistantText(msg.content);
          finalizeAssistantMessage();
        }
      });

      session.visualizations.forEach(function (viz) {
        renderVisualization({
          type: viz.viz_type,
          title: viz.title,
          data: viz.data,
          viz_id: viz.id
        });
      });

      var msgs = document.getElementById("chat-messages");
      if (msgs) msgs.scrollTop = msgs.scrollHeight;
    });
  }

  function loadSessionList() {
    var baseUrl = getBaseUrl();
    var authToken = sessionStorage.getItem("auth_token");
    if (!authToken) return;

    fetch(baseUrl + "/v1/chat/sessions", {
      headers: getAuthHeaders()
    })
    .then(function (r) { return r.ok ? r.json() : { sessions: [] }; })
    .then(function (data) {
      renderSessionList(data.sessions || []);
    })
    .catch(function () {});
  }

  function renderSessionList(sessions) {
    var container = document.getElementById("chat-session-list");
    if (!container) return;
    container.innerHTML = "";

    if (sessions.length === 0) {
      var empty = document.createElement("div");
      empty.className = "session-empty-state";
      empty.textContent = "No previous conversations";
      container.appendChild(empty);
      return;
    }

    var now = new Date();
    var todayStr = now.toDateString();
    var yesterday = new Date(now);
    yesterday.setDate(yesterday.getDate() - 1);
    var yesterdayStr = yesterday.toDateString();

    var lastGroup = "";
    sessions.forEach(function (s) {
      var d = new Date(s.updated_at);
      var group;
      if (d.toDateString() === todayStr) group = "Today";
      else if (d.toDateString() === yesterdayStr) group = "Yesterday";
      else group = "Older";

      if (group !== lastGroup) {
        var label = document.createElement("div");
        label.className = "session-date-group";
        label.textContent = group;
        container.appendChild(label);
        lastGroup = group;
      }

      var item = document.createElement("div");
      item.className = "chat-session-item";
      if (s.id === sessionId) item.classList.add("active");
      item.setAttribute("data-session-id", s.id);

      var icon = document.createElement("i");
      icon.className = "bi bi-chat-left-text";
      icon.style.fontSize = "0.75rem";
      icon.style.flexShrink = "0";
      item.appendChild(icon);

      var titleSpan = document.createElement("span");
      titleSpan.className = "session-title";
      titleSpan.textContent = s.title || "New conversation";
      titleSpan.title = s.title || "New conversation";
      item.appendChild(titleSpan);

      // Delete button
      var delBtn = document.createElement("button");
      delBtn.className = "session-delete-btn";
      delBtn.title = "Delete";
      delBtn.innerHTML = '<i class="bi bi-x"></i>';
      delBtn.addEventListener("click", function (e) {
        e.stopPropagation();
        deleteSession(s.id);
      });
      item.appendChild(delBtn);

      // Click to switch session
      item.addEventListener("click", function () {
        switchToSession(s.id);
      });

      // Double-click to rename
      titleSpan.addEventListener("dblclick", function (e) {
        e.stopPropagation();
        startRenameSession(s.id, titleSpan);
      });

      container.appendChild(item);
    });
  }

  function switchToSession(id) {
    if (id === sessionId) return;
    sessionId = id;
    sessionStorage.setItem("chat_session_id", id);
    clearChat();
    clearDashboard();
    loadSession(id);
    loadSessionList();
  }

  function deleteSession(id) {
    if (!confirm("Delete this conversation?")) return;
    var baseUrl = getBaseUrl();
    fetch(baseUrl + "/v1/chat/sessions/" + id, {
      method: "DELETE",
      headers: getAuthHeaders()
    }).then(function () {
      loadSessionList();
      if (id === sessionId) {
        sessionId = null;
        sessionStorage.removeItem("chat_session_id");
        clearChat();
        clearDashboard();
        showWelcome();
      }
    });
  }

  function startRenameSession(id, titleSpan) {
    var current = titleSpan.textContent;
    var input = document.createElement("input");
    input.type = "text";
    input.value = current;
    input.className = "session-rename-input";
    titleSpan.textContent = "";
    titleSpan.appendChild(input);
    input.focus();
    input.select();

    function finish() {
      var newTitle = (input.value || "").trim();
      if (!newTitle || newTitle === current) {
        titleSpan.textContent = current;
        return;
      }
      titleSpan.textContent = newTitle;
      var baseUrl = getBaseUrl();
      fetch(baseUrl + "/v1/chat/sessions/" + id, {
        method: "PATCH",
        headers: Object.assign({ "Content-Type": "application/json" }, getAuthHeaders()),
        body: JSON.stringify({ title: newTitle })
      });
    }

    input.addEventListener("blur", finish);
    input.addEventListener("keydown", function (e) {
      if (e.key === "Enter") { e.preventDefault(); input.blur(); }
      if (e.key === "Escape") { titleSpan.textContent = current; }
    });
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
