document.addEventListener("DOMContentLoaded", function () {
  const controls = document.getElementById("qsm-controls");
  if (!controls) return;

  const phone = window.matchMedia("(max-width: 767px)");
  function setControlsLayout() { controls.open = !phone.matches; }
  setControlsLayout();
  phone.addEventListener("change", setControlsLayout);

  // Leaflet must remeasure its canvas after a tab or iframe changes size.
  function resizeMaps() {
    if (!window.HTMLWidgets) return;
    ["baseline_map", "results_map"].forEach(function (id) {
      const node = document.getElementById(id);
      if (!node || !node.clientWidth || !node.clientHeight) return;
      const widget = window.HTMLWidgets.find("#" + id);
      const map = widget && widget.getMap ? widget.getMap() : null;
      if (map) map.invalidateSize({ pan: false });
    });
  }
  const observer = new ResizeObserver(function () { window.requestAnimationFrame(resizeMaps); });
  ["baseline_map", "results_map"].forEach(function (id) {
    const node = document.getElementById(id);
    if (node) observer.observe(node);
  });
  if (window.jQuery) {
    window.jQuery(document).on("shown.bs.tab shiny:connected", function () {
      window.requestAnimationFrame(resizeMaps);
    });
  }
});
