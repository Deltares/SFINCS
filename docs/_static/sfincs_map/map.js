// Interactive map of countries where SFINCS has been applied.
// Data: applications.json (edit this) + countries.geojson (Natural Earth outlines, static).
(function () {
  var el = document.getElementById("sfincs-map");
  if (!el || typeof L === "undefined") return;
  var base = el.getAttribute("data-base") || "_static/sfincs_map/";

  var color = "#0a28a3";  // Deltares blue
  var map = L.map(el, {
    worldCopyJump: true,
    minZoom: 1,
    zoomControl: false,
    maxBounds: [[-85, -540], [85, 540]],  // no grey band beyond the poles
    maxBoundsViscosity: 1.0
  }).setView([20, 10], 2);
  L.control.zoom({ position: "topright" }).addTo(map);

  var logo = L.control({ position: "topleft" });
  logo.onAdd = function () {
    var d = L.DomUtil.create("div", "sfincs-map-logo");
    d.innerHTML = "<img src='" + base + "sfincs_logo.png' alt='SFINCS'>";
    return d;
  };
  logo.addTo(map);

  // Added after the zoom buttons; map.css puts it to their left
  var deltares = L.control({ position: "topright" });
  deltares.onAdd = function () {
    var d = L.DomUtil.create("div", "sfincs-map-logo sfincs-map-logo-deltares");
    d.innerHTML = "<a href='https://www.deltares.nl' target='_blank' rel='noopener'>" +
      "<img src='" + base + "deltares_logo.svg' alt='Deltares'></a>";
    return d;
  };
  deltares.addTo(map);

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

  // Land area in square degrees of a polygon (outer ring only)
  function ringArea(ring) {
    var a = 0;
    for (var i = 0, j = ring.length - 1; i < ring.length; j = i++) {
      a += (ring[j][0] + ring[i][0]) * (ring[j][1] - ring[i][1]);
    }
    return Math.abs(a / 2);
  }

  // For countries with little land (< 2 square degrees), return the centre of the
  // largest island as [lat, lon]; otherwise null. Works for island groups spread
  // over a large area or across the date line (Kiribati, French Polynesia).
  function smallIslandPoint(geom) {
    var polys = geom.type === "Polygon" ? [geom.coordinates] : geom.coordinates;
    var total = 0, best = null, bestArea = -1;
    polys.forEach(function (p) {
      var a = ringArea(p[0]);
      total += a;
      if (a > bestArea) { bestArea = a; best = p[0]; }
    });
    if (total >= 2) return null;
    var b = L.latLngBounds(best.map(function (c) { return [c[1], c[0]]; }));
    return b.getCenter();
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
        var small = smallIslandPoint(f.geometry);
        if (small) {
          L.circleMarker(small, { radius: 4, color: color, weight: 1, fillOpacity: 0.8 })
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
        Object.keys(found).length + " countries and territories)";
      return d;
    };
    legend.addTo(map);
  }).catch(function (e) {
    el.innerHTML = "<p>Map could not be loaded.</p>";
    console.error("sfincs_map:", e);
  });
})();
