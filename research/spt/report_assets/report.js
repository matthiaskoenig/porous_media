"use strict";
(() => {
  const cards = [...document.querySelectorAll(".figure-card")];
  const groups = [...document.querySelectorAll(".figure-section")];
  const controls = ["section", "family", "geometry", "scale", "pattern", "dose"];
  const search = document.getElementById("search");
  const section = document.getElementById("section");
  const viewer = document.getElementById("viewer");
  const image = document.getElementById("viewer-image");
  const zoom = document.getElementById("zoom");
  let current = null;
  let opener = null;
  const visible = () => cards.filter(card => !card.hidden && card.querySelector("[data-inspect]"));
  function filter() {
    const words = search.value.toLowerCase().trim().split(/\s+/).filter(Boolean);
    cards.forEach(card => {
      card.hidden = !words.every(word => card.dataset.search.includes(word)) || controls.some(key => {
        const value = document.getElementById(key).value;
        if (!value) return false;
        return ["pattern", "dose"].includes(key)
          ? !JSON.parse(card.dataset[key]).includes(value)
          : card.dataset[key] !== value;
      });
    });
    groups.forEach(group => {
      const count = group.querySelectorAll(".figure-card:not([hidden])").length;
      group.hidden = count === 0;
      group.querySelector(".section-count").textContent = `${count} figure${count === 1 ? "" : "s"}`;
    });
    const count = cards.filter(card => !card.hidden).length;
    document.getElementById("result-count").textContent = `${count} of ${cards.length} figures`;
    document.getElementById("empty").hidden = count !== 0;
    document.querySelectorAll("[data-section][href]").forEach(link => {
      if (link.dataset.section === section.value) link.setAttribute("aria-current", "page");
      else link.removeAttribute("aria-current");
    });
  }
  search.addEventListener("input", filter);
  controls.forEach(key => document.getElementById(key).addEventListener("change", filter));
  document.getElementById("reset").addEventListener("click", () => {
    search.value = "";
    controls.forEach(key => { document.getElementById(key).value = ""; });
    filter();
  });
  document.querySelectorAll("[data-section][href]").forEach(link => link.addEventListener("click", () => {
    section.value = link.dataset.section;
    filter();
  }));
  function fit() {
    image.classList.add("fit");
    image.style.width = "";
    zoom.value = "100";
    document.getElementById("canvas").scrollTo(0, 0);
  }
  function magnify() {
    image.classList.remove("fit");
    image.style.width = `${image.naturalWidth * Number(zoom.value) / 100}px`;
  }
  function inspect(card) {
    current = card;
    const all = visible();
    document.getElementById("viewer-position").textContent = `Figure ${all.indexOf(card) + 1} of ${all.length} in selection`;
    document.getElementById("viewer-title").textContent = card.dataset.title;
    document.getElementById("viewer-caption").textContent = card.querySelector(".card-body p").textContent;
    image.alt = card.querySelector(".preview").getAttribute("aria-label");
    image.src = card.querySelector("[data-inspect]").dataset.inspect;
    const exports = document.getElementById("viewer-exports");
    exports.replaceChildren(...[...card.querySelectorAll(".export-links a:not(.permalink)")].map(link => link.cloneNode(true)));
    fit();
    if (!viewer.open) viewer.showModal();
  }
  document.querySelectorAll("[data-inspect]").forEach(link => link.addEventListener("click", event => {
    if (event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    opener = link;
    inspect(link.closest(".figure-card"));
  }));
  function step(direction) {
    const all = visible();
    if (all.length) inspect(all[(all.indexOf(current) + direction + all.length) % all.length]);
  }
  document.getElementById("close").addEventListener("click", () => viewer.close());
  viewer.addEventListener("close", () => {
    image.removeAttribute("src");
    if (opener) opener.focus({preventScroll: true});
  });
  document.getElementById("previous").addEventListener("click", () => step(-1));
  document.getElementById("next").addEventListener("click", () => step(1));
  document.getElementById("fit").addEventListener("click", fit);
  document.getElementById("actual").addEventListener("click", () => { zoom.value = "100"; magnify(); });
  zoom.addEventListener("input", magnify);
  image.addEventListener("load", () => { if (!image.classList.contains("fit")) magnify(); });
  document.addEventListener("keydown", event => {
    const editing = ["INPUT", "SELECT", "TEXTAREA"].includes(event.target.tagName);
    if (viewer.open && !editing && ["ArrowLeft", "ArrowRight"].includes(event.key)) {
      event.preventDefault();
      step(event.key === "ArrowRight" ? 1 : -1);
    } else if (!viewer.open && !editing && event.key === "/") {
      event.preventDefault();
      search.focus();
    }
  });
  // A copied figure permalink must reveal its card even after filtering.
  function revealHash() {
    const target = document.getElementById(location.hash.slice(1));
    if (target?.classList.contains("figure-card") && target.hidden) {
      document.getElementById("reset").click();
      target.scrollIntoView();
    }
  }
  window.addEventListener("hashchange", revealHash);
  filter();
  revealHash();
})();
