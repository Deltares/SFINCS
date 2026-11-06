// Interactive map of countries where SFINCS has been applied.
// Data: applications.json (edit this) + countries.geojson (Natural Earth outlines, static).
(function () {
  var el = document.getElementById("sfincs-map");
  if (!el || typeof L === "undefined") return;
  var base = el.getAttribute("data-base") || "_static/sfincs_map/";

  var color = "#cc3333";
  var map = L.map(el, { worldCopyJump: true, minZoom: 1 }).setView([20, 10], 2);

  L.tileLayer("https://tile.openstreetmap.org/{z}/{x}/{y}.png", {
    maxZoom: 10,
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors | ' +
      'Borders: <a href="https://www.naturalearthdata.com">Natural Earth</a> | ' +
      '<a href="https://www.deltares.nl">Deltares</a>'
  }).addTo(map);

  function getJson(name) {
    return fetch(base + name).then(function (r) { return r.json(); });
  }

  function popupHtml(name, links) {
    var html = "<b>" + name + "</b>";
    if (links && links.length) {
      html += "<ul style='margin:4px 0 0 0;padding-left:16px'>";
      links.forEach(function (l) {
        html += "<li><a href='" + l.url + "' target='_blank' rel='noopener'>" + l.title + "</a></li>";
      });
      html += "</ul>";
    }
    return html;
  }

  Promise.all([getJson("applications.json"), getJson("countries.geojson")]).then(function (res) {
    var apps = {};
    res[0].countries.forEach(function (c) { apps[c.iso] = c; });

    var found = {};
    var style = { color: color, weight: 1, fillColor: color, fillOpacity: 0.5 };
    var layer = L.geoJSON(res[1], {
      filter: function (f) { return f.properties.iso in apps; },
      style: style,
      onEachFeature: function (f, lyr) {
        var c = apps[f.properties.iso];
        var name = c.name || f.properties.name;
        found[f.properties.iso] = true;
        lyr.bindPopup(popupHtml(name, c.links));
        lyr.bindTooltip(name, { sticky: true });
        lyr.on("mouseover", function () { lyr.setStyle({ fillOpacity: 0.8 }); });
        lyr.on("mouseout", function () { lyr.setStyle({ fillOpacity: 0.5 }); });

        // Small island states: add a marker so they are visible and clickable at world zoom
        var b = lyr.getBounds();
        if (b.getNorth() - b.getSouth() < 1.5 && b.getEast() - b.getWest() < 1.5) {
          L.circleMarker(b.getCenter(), { radius: 4, color: color, weight: 1, fillOpacity: 0.8 })
            .bindPopup(popupHtml(name, c.links))
            .bindTooltip(name)
            .addTo(map);
        }
      }
    }).addTo(map);
    map.fitBounds(layer.getBounds(), { padding: [10, 10] });

    Object.keys(apps).forEach(function (iso) {
      if (!found[iso]) console.warn("sfincs_map: ISO code not found in countries.geojson: " + iso);
    });

    var legend = L.control({ position: "bottomleft" });
    legend.onAdd = function () {
      var d = L.DomUtil.create("div", "sfincs-map-legend");
      d.innerHTML = "<span style='background:" + color + "'></span>SFINCS applied (" +
        Object.keys(found).length + " countries)";
      return d;
    };
    legend.addTo(map);
  }).catch(function (e) {
    el.innerHTML = "<p>Map could not be loaded.</p>";
    console.error("sfincs_map:", e);
  });
})();
