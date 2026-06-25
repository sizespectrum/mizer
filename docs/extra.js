// Make a click on the "Articles" navbar entry go to the articles index page
// instead of merely toggling its drop-down menu (the drop-down still opens on
// hover, see extra.css). The toggle is rendered by pkgdown as a <button>, so we
// reuse the relative href of the "More articles..." link inside the menu, which
// already points at articles/index.html with the correct depth for this page.
(function() {
  function setup() {
    var toggle = document.getElementById("dropdown-articles");
    if (!toggle) return;

    var menu = document.querySelector('[aria-labelledby="dropdown-articles"]');
    if (!menu) return;

    var indexHref = null;
    var links = menu.querySelectorAll("a.dropdown-item");
    for (var i = 0; i < links.length; i++) {
      var href = links[i].getAttribute("href");
      if (href && /articles\/index\.html$/.test(href)) {
        indexHref = href;
      }
    }
    if (!indexHref) return;

    toggle.addEventListener("click", function() {
      window.location.href = indexHref;
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", setup);
  } else {
    setup();
  }
})();
