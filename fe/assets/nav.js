(function () {
  "use strict";

  function normalizePath(p) {
    return p.length > 1 && p.endsWith("/") ? p.slice(0, -1) : p;
  }

  function updateActiveLink() {
    var path = normalizePath(window.location.pathname);
    var links = document.querySelectorAll(".header-link");
    if (!links.length) return;
    links.forEach(function (link) {
      var href = link.getAttribute("href");
      if (!href || link.target === "_blank" || href.startsWith("http")) {
        link.classList.remove("active");
        return;
      }
      var normalizedHref = normalizePath(href);
      var isActive = path === normalizedHref || path.startsWith(normalizedHref + "/");
      link.classList.toggle("active", isActive);
    });
  }

  // Observe the entire document until header links appear, then keep updating on navigation
  var bodyObserver = new MutationObserver(function () {
    if (document.querySelectorAll(".header-link").length) {
      updateActiveLink();
    }
  });

  bodyObserver.observe(document.documentElement, { childList: true, subtree: true });

  window.addEventListener("popstate", updateActiveLink);
})();
