// Interactive map of countries where SFINCS has been applied, in the Equal Earth projection.
// Data: applications.json (edit this) + countries.geojson (Natural Earth outlines, static).
// Drawn with D3 (d3-geo, d3-zoom); no basemap tiles, as web tiles only exist in Web Mercator.
(function () {
  var el = document.getElementById("sfincs-map");
  if (!el || typeof d3 === "undefined") return;
  var base = el.getAttribute("data-base") || "_static/sfincs_map/";

  var color = "#0a28a3";  // Deltares blue
  var islandRadius = 3.5;

  // Overlays: logos, zoom buttons, legend, tooltip, popup, attribution
  el.innerHTML =
    "<svg role='img' aria-label='World map of countries where SFINCS has been applied'></svg>" +
    "<div class='sfincs-map-logo sfincs-map-logo-sfincs'><img src='" + base + "sfincs_logo.png' alt='SFINCS'></div>" +
    "<div class='sfincs-map-topright'>" +
      "<div class='sfincs-map-logo sfincs-map-logo-deltares'><a href='https://www.deltares.nl' target='_blank' rel='noopener'>" +
        "<img src='" + base + "deltares_logo.svg' alt='Deltares'></a></div>" +
      "<div class='sfincs-map-zoom'><button type='button' aria-label='Zoom in'>+</button>" +
        "<button type='button' aria-label='Zoom out'>&minus;</button></div>" +
    "</div>" +
    "<div class='sfincs-map-legend' hidden></div>" +
    "<div class='sfincs-map-attribution'>Projection: <a href='https://equal-earth.com' target='_blank' rel='noopener'>Equal Earth</a> | " +
      "Borders: <a href='https://www.naturalearthdata.com' target='_blank' rel='noopener'>Natural Earth</a> | " +
      "<a href='https://www.deltares.nl' target='_blank' rel='noopener'>Deltares</a></div>" +
    "<div class='sfincs-map-tooltip' hidden></div>" +
    "<div class='sfincs-map-popup' hidden></div>";

  var svg = d3.select(el).select("svg");
  var g = svg.append("g");
  var tooltip = el.querySelector(".sfincs-map-tooltip");
  var popup = el.querySelector(".sfincs-map-popup");
  var legend = el.querySelector(".sfincs-map-legend");
  var zoomButtons = el.querySelectorAll(".sfincs-map-zoom button");

  function getJson(name) {
    return fetch(base + name).then(function (r) { return r.json(); });
  }

  function linkList(items) {
    var html = "<ul>";
    items.forEach(function (it) {
      html += "<li><a href='" + it.url + "' target='_blank' rel='noopener'>" + it.text + "</a></li>";
    });
    return html + "</ul>";
  }

  // Popup: country name, news/website links first, then papers (newest first)
  function popupHtml(name, c) {
    var html = "<button type='button' class='sfincs-map-popup-close' aria-label='Close'>&times;</button>" +
      "<b>" + name + "</b>";
    if (c.links && c.links.length) {
      html += linkList(c.links.map(function (l) { return { url: l.url, text: l.title }; }));
    }
    if (c.papers && c.papers.length) {
      html += "<div class='sfincs-map-papers'>Papers and abstracts:</div>";
      html += linkList(c.papers.map(function (p) {
        return { url: p.url || "https://doi.org/" + p.doi, text: p.title + " (" + p.year + ")" };
      }));
    }
    return html;
  }

  // Signed planar area of a ring in square degrees (positive = counter-clockwise)
  function signedArea(ring) {
    var a = 0;
    for (var i = 0, j = ring.length - 1; i < ring.length; j = i++) {
      a += (ring[j][0] + ring[i][0]) * (ring[j][1] - ring[i][1]);
    }
    return -a / 2;
  }

  // d3-geo is spherical: an exterior ring must be clockwise, a hole counter-clockwise.
  // countries.geojson follows the GeoJSON convention (the opposite), so reverse where needed;
  // otherwise d3 fills the whole globe except the country.
  function rewind(geom) {
    var polys = geom.type === "Polygon" ? [geom.coordinates] : geom.coordinates;
    polys.forEach(function (p) {
      p.forEach(function (ring, k) {
        var clockwise = signedArea(ring) < 0;
        if ((k === 0) !== clockwise) ring.reverse();
      });
    });
  }

  // For countries with little land (< 2 square degrees), return the centre of the
  // largest island as [lon, lat]; otherwise null. Works for island groups spread
  // over a large area or across the date line (Kiribati, French Polynesia).
  function smallIslandPoint(geom) {
    var polys = geom.type === "Polygon" ? [geom.coordinates] : geom.coordinates;
    var total = 0, best = null, bestArea = -1;
    polys.forEach(function (p) {
      var a = Math.abs(signedArea(p[0]));
      total += a;
      if (a > bestArea) { bestArea = a; best = p[0]; }
    });
    if (total >= 2) return null;
    var lon = d3.extent(best, function (c) { return c[0]; });
    var lat = d3.extent(best, function (c) { return c[1]; });
    return [(lon[0] + lon[1]) / 2, (lat[0] + lat[1]) / 2];
  }

  function showTooltip(ev, name) {
    var b = el.getBoundingClientRect();
    tooltip.textContent = name;
    tooltip.style.left = (ev.clientX - b.left + 12) + "px";
    tooltip.style.top = (ev.clientY - b.top + 12) + "px";
    tooltip.hidden = false;
  }

  function showPopup(ev, html) {
    var b = el.getBoundingClientRect();
    popup.innerHTML = html;
    popup.hidden = false;
    var x = Math.min(ev.clientX - b.left + 10, b.width - popup.offsetWidth - 10);
    var y = Math.min(ev.clientY - b.top + 10, b.height - popup.offsetHeight - 10);
    popup.style.left = Math.max(10, x) + "px";
    popup.style.top = Math.max(10, y) + "px";
    popup.querySelector(".sfincs-map-popup-close").onclick = function () { popup.hidden = true; };
  }

  var zoom = d3.zoom().scaleExtent([1, 12]).on("zoom", function (ev) {
    g.attr("transform", ev.transform);
    g.selectAll(".sfincs-map-island").attr("r", islandRadius / Math.sqrt(ev.transform.k));
  });
  svg.call(zoom).on("click", function () { popup.hidden = true; });
  zoomButtons[0].onclick = function () { svg.transition().duration(250).call(zoom.scaleBy, 1.6); };
  zoomButtons[1].onclick = function () { svg.transition().duration(250).call(zoom.scaleBy, 1 / 1.6); };

  Promise.all([getJson("applications.json"), getJson("countries.geojson")]).then(function (res) {
    var apps = {};
    res[0].countries.forEach(function (c) { apps[c.iso] = c; });

    var world = res[1];
    world.features.forEach(function (f) { rewind(f.geometry); });

    var found = {};
    var applied = world.features.filter(function (f) {
      if (!(f.properties.iso in apps)) return false;
      found[f.properties.iso] = true;
      return true;
    });

    var islands = [];
    applied.forEach(function (f) {
      var p = smallIslandPoint(f.geometry);
      if (p) islands.push({ feature: f, point: p });
    });

    function countryName(f) { return apps[f.properties.iso].name || f.properties.name; }

    function interactive(sel, featureOf) {
      sel.on("mousemove", function (ev, d) { showTooltip(ev, countryName(featureOf(d))); })
        .on("mouseleave", function () { tooltip.hidden = true; })
        .on("click", function (ev, d) {
          var f = featureOf(d);
          showPopup(ev, popupHtml(countryName(f), apps[f.properties.iso]));
          ev.stopPropagation();
        });
    }

    // (Re)draw at the current size of the map element
    function draw() {
      var w = el.clientWidth, h = el.clientHeight;
      svg.attr("viewBox", "0 0 " + w + " " + h);
      var projection = d3.geoEqualEarth().fitExtent([[8, 8], [w - 8, h - 8]], { type: "Sphere" });
      var path = d3.geoPath(projection);

      g.selectAll("*").remove();
      g.append("path").attr("class", "sfincs-map-sphere").attr("d", path({ type: "Sphere" }));
      g.append("path").attr("class", "sfincs-map-graticule").attr("d", path(d3.geoGraticule10()));
      g.append("g").selectAll("path").data(world.features).join("path")
        .attr("class", "sfincs-map-land").attr("d", path);
      g.append("g").selectAll("path").data(applied).join("path")
        .attr("class", "sfincs-map-applied").attr("d", path)
        .call(interactive, function (f) { return f; });

      // Small island states: add a marker so they are visible and clickable at world zoom
      g.append("g").selectAll("circle").data(islands).join("circle")
        .attr("class", "sfincs-map-island")
        .attr("cx", function (d) { return projection(d.point)[0]; })
        .attr("cy", function (d) { return projection(d.point)[1]; })
        .attr("r", islandRadius / Math.sqrt(d3.zoomTransform(svg.node()).k))
        .call(interactive, function (d) { return d.feature; });

      zoom.extent([[0, 0], [w, h]]).translateExtent([[0, 0], [w, h]]);
    }

    draw();
    var resizeTimer;
    window.addEventListener("resize", function () {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () {
        popup.hidden = true;
        svg.call(zoom.transform, d3.zoomIdentity);
        draw();
      }, 150);
    });

    Object.keys(apps).forEach(function (iso) {
      if (!found[iso]) console.warn("sfincs_map: ISO code not found in countries.geojson: " + iso);
    });

    legend.innerHTML = "<span style='background:" + color + "'></span>SFINCS applied (" +
      Object.keys(found).length + " countries and territories)";
    legend.hidden = false;
  }).catch(function (e) {
    el.innerHTML = "<p>Map could not be loaded.</p>";
    console.error("sfincs_map:", e);
  });
})();
