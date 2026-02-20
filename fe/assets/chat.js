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

    // Suggestion buttons
    var suggestions = document.querySelectorAll(".chat-suggestion-btn");
    suggestions.forEach(function (btn) {
      btn.addEventListener("click", function () {
        var query = btn.getAttribute("data-query") || btn.textContent;
        input.value = query;
        sendCurrentMessage();
      });
    });
  }

  function sendCurrentMessage() {
    var input = document.getElementById("chat-input");
    var text = (input.value || "").trim();
    if (!text) return;
    // Use native setter + event to sync React/Dash controlled state.
    var nativeSetter = Object.getOwnPropertyDescriptor(
      window.HTMLInputElement.prototype, "value"
    ).set;
    nativeSetter.call(input, "");
    input.dispatchEvent(new Event("input", { bubbles: true }));
    hideSuggestions();
    sendMessage(text);
  }

  function hideSuggestions() {
    var el = document.getElementById("chat-suggestions");
    if (el) el.style.display = "none";
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

  function setInputEnabled(enabled) {
    var input = document.getElementById("chat-input");
    var btn = document.getElementById("chat-send-btn");
    if (input) input.disabled = !enabled;
    if (btn) btn.disabled = !enabled;
    if (enabled && input) input.focus();
  }

  // --- Simple markdown formatting ---

  function formatMarkdown(text) {
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

  function updatePortalVisibility() {
    var portal = document.getElementById("chat-data-portal");
    var clearBtn = document.getElementById("chat-clear-portal-btn");
    var ph = document.getElementById("chat-portal-placeholder");
    var hasItems = portal && portal.children.length > 0;
    if (clearBtn) clearBtn.style.display = hasItems ? "" : "none";
    if (ph) ph.style.display = hasItems ? "none" : "";
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

    var wrapper = document.createElement("div");
    wrapper.className = "viz-container";

    var header = document.createElement("div");
    header.className = "viz-header";

    var title = document.createElement("h6");
    title.className = "viz-title";
    title.textContent = data.title || "Visualization";
    header.appendChild(title);

    var removeBtn = document.createElement("button");
    removeBtn.className = "viz-remove-btn";
    var removeIcon = document.createElement("i");
    removeIcon.className = "bi bi-x-circle-fill";
    removeBtn.appendChild(removeIcon);
    removeBtn.addEventListener("click", function () {
      wrapper.remove();
      updatePortalVisibility();
    });
    header.appendChild(removeBtn);

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
    chartDiv.style.maxWidth = "500px";
    chartDiv.style.margin = "0 auto";
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
        height: 350
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
        height: 350
      },
      { responsive: true, displayModeBar: false }
    );
  }

  // --- Initialize when chat page is rendered ---

  // Use MutationObserver to detect when Dash renders the chat page
  var observer = new MutationObserver(function () {
    init();
  });

  // Start observing
  if (document.body) {
    observer.observe(document.body, { childList: true, subtree: true });
  }

  // Also try immediately in case page is already rendered
  if (document.readyState === "complete" || document.readyState === "interactive") {
    setTimeout(init, 100);
  } else {
    document.addEventListener("DOMContentLoaded", function () {
      setTimeout(init, 100);
    });
  }
})();
