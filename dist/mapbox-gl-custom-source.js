(function(f){if(typeof exports==="object"&&typeof module!=="undefined"){module.exports=f()}else if(typeof define==="function"&&define.amd){define([],f)}else{var g;if(typeof window!=="undefined"){g=window}else if(typeof global!=="undefined"){g=global}else if(typeof self!=="undefined"){g=self}else{g=this}g.CustomGeoJSONSource = f()}})(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
var GeoJSONSource = require('mapbox-gl/src/source/geojson_source');
var webworkify = require('webworkify');

class CustomGeoJSONSource extends GeoJSONSource {
}

CustomGeoJSONSource.workerSourceURL = URL.createObjectURL(webworkify(require('./custom_geojson_worker.js'), {bare: true}));

module.exports = CustomGeoJSONSource;

},{"./custom_geojson_worker.js":2,"mapbox-gl/src/source/geojson_source":24,"webworkify":51}],2:[function(require,module,exports){
var GeoJSONWorkerSource = require('mapbox-gl/src/source/geojson_worker_source')

function customLoadGeoJSON (options, callback) {
  GeoJSONWorkerSource.prototype.loadGeoJSON.call(this, options, function (err, data) {
    // Read the data normally, then call this custom callback, which shows options and data in console
    if (err) { return callback(err) }
    console.log(options, data)
    callback(null, data)
  })
}

class CustomGeoJSONWorkerSource extends GeoJSONWorkerSource {
  constructor(actor, style) {
      super(actor, style, customLoadGeoJSON)
  }
}

module.exports = function(self) {
  self.registerWorkerSource('custom_geojson', CustomGeoJSONWorkerSource)
}

},{"mapbox-gl/src/source/geojson_worker_source":25}],3:[function(require,module,exports){
/*
 * Copyright (C) 2008 Apple Inc. All Rights Reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY APPLE INC. ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL APPLE INC. OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Ported from Webkit
 * http://svn.webkit.org/repository/webkit/trunk/Source/WebCore/platform/graphics/UnitBezier.h
 */

module.exports = UnitBezier;

function UnitBezier(p1x, p1y, p2x, p2y) {
    // Calculate the polynomial coefficients, implicit first and last control points are (0,0) and (1,1).
    this.cx = 3.0 * p1x;
    this.bx = 3.0 * (p2x - p1x) - this.cx;
    this.ax = 1.0 - this.cx - this.bx;

    this.cy = 3.0 * p1y;
    this.by = 3.0 * (p2y - p1y) - this.cy;
    this.ay = 1.0 - this.cy - this.by;

    this.p1x = p1x;
    this.p1y = p2y;
    this.p2x = p2x;
    this.p2y = p2y;
}

UnitBezier.prototype.sampleCurveX = function(t) {
    // `ax t^3 + bx t^2 + cx t' expanded using Horner's rule.
    return ((this.ax * t + this.bx) * t + this.cx) * t;
};

UnitBezier.prototype.sampleCurveY = function(t) {
    return ((this.ay * t + this.by) * t + this.cy) * t;
};

UnitBezier.prototype.sampleCurveDerivativeX = function(t) {
    return (3.0 * this.ax * t + 2.0 * this.bx) * t + this.cx;
};

UnitBezier.prototype.solveCurveX = function(x, epsilon) {
    if (typeof epsilon === 'undefined') epsilon = 1e-6;

    var t0, t1, t2, x2, i;

    // First try a few iterations of Newton's method -- normally very fast.
    for (t2 = x, i = 0; i < 8; i++) {

        x2 = this.sampleCurveX(t2) - x;
        if (Math.abs(x2) < epsilon) return t2;

        var d2 = this.sampleCurveDerivativeX(t2);
        if (Math.abs(d2) < 1e-6) break;

        t2 = t2 - x2 / d2;
    }

    // Fall back to the bisection method for reliability.
    t0 = 0.0;
    t1 = 1.0;
    t2 = x;

    if (t2 < t0) return t0;
    if (t2 > t1) return t1;

    while (t0 < t1) {

        x2 = this.sampleCurveX(t2);
        if (Math.abs(x2 - x) < epsilon) return t2;

        if (x > x2) {
            t0 = t2;
        } else {
            t1 = t2;
        }

        t2 = (t1 - t0) * 0.5 + t0;
    }

    // Failure.
    return t2;
};

UnitBezier.prototype.solve = function(x, epsilon) {
    return this.sampleCurveY(this.solveCurveX(x, epsilon));
};

},{}],4:[function(require,module,exports){
var wgs84 = require('wgs84');

module.exports.geometry = geometry;
module.exports.ring = ringArea;

function geometry(_) {
    if (_.type === 'Polygon') return polygonArea(_.coordinates);
    else if (_.type === 'MultiPolygon') {
        var area = 0;
        for (var i = 0; i < _.coordinates.length; i++) {
            area += polygonArea(_.coordinates[i]);
        }
        return area;
    } else {
        return null;
    }
}

function polygonArea(coords) {
    var area = 0;
    if (coords && coords.length > 0) {
        area += Math.abs(ringArea(coords[0]));
        for (var i = 1; i < coords.length; i++) {
            area -= Math.abs(ringArea(coords[i]));
        }
    }
    return area;
}

/**
 * Calculate the approximate area of the polygon were it projected onto
 *     the earth.  Note that this area will be positive if ring is oriented
 *     clockwise, otherwise it will be negative.
 *
 * Reference:
 * Robert. G. Chamberlain and William H. Duquette, "Some Algorithms for
 *     Polygons on a Sphere", JPL Publication 07-03, Jet Propulsion
 *     Laboratory, Pasadena, CA, June 2007 http://trs-new.jpl.nasa.gov/dspace/handle/2014/40409
 *
 * Returns:
 * {float} The approximate signed geodesic area of the polygon in square
 *     meters.
 */

function ringArea(coords) {
    var area = 0;

    if (coords.length > 2) {
        var p1, p2;
        for (var i = 0; i < coords.length - 1; i++) {
            p1 = coords[i];
            p2 = coords[i + 1];
            area += rad(p2[0] - p1[0]) * (2 + Math.sin(rad(p1[1])) + Math.sin(rad(p2[1])));
        }

        area = area * wgs84.RADIUS * wgs84.RADIUS / 2;
    }

    return area;
}

function rad(_) {
    return _ * Math.PI / 180;
}

},{"wgs84":52}],5:[function(require,module,exports){
var geojsonArea = require('geojson-area');

module.exports = rewind;

function rewind(gj, outer) {
    switch ((gj && gj.type) || null) {
        case 'FeatureCollection':
            gj.features = gj.features.map(curryOuter(rewind, outer));
            return gj;
        case 'Feature':
            gj.geometry = rewind(gj.geometry, outer);
            return gj;
        case 'Polygon':
        case 'MultiPolygon':
            return correct(gj, outer);
        default:
            return gj;
    }
}

function curryOuter(a, b) {
    return function(_) { return a(_, b); };
}

function correct(_, outer) {
    if (_.type === 'Polygon') {
        _.coordinates = correctRings(_.coordinates, outer);
    } else if (_.type === 'MultiPolygon') {
        _.coordinates = _.coordinates.map(curryOuter(correctRings, outer));
    }
    return _;
}

function correctRings(_, outer) {
    outer = !!outer;
    _[0] = wind(_[0], !outer);
    for (var i = 1; i < _.length; i++) {
        _[i] = wind(_[i], outer);
    }
    return _;
}

function wind(_, dir) {
    return cw(_) === dir ? _ : _.reverse();
}

function cw(_) {
    return geojsonArea.ring(_) >= 0;
}

},{"geojson-area":4}],6:[function(require,module,exports){
'use strict';

module.exports = clip;

var createFeature = require('./feature');

/* clip features between two axis-parallel lines:
 *     |        |
 *  ___|___     |     /
 * /   |   \____|____/
 *     |        |
 */

function clip(features, scale, k1, k2, axis, intersect, minAll, maxAll) {

    k1 /= scale;
    k2 /= scale;

    if (minAll >= k1 && maxAll <= k2) return features; // trivial accept
    else if (minAll > k2 || maxAll < k1) return null; // trivial reject

    var clipped = [];

    for (var i = 0; i < features.length; i++) {

        var feature = features[i],
            geometry = feature.geometry,
            type = feature.type,
            min, max;

        min = feature.min[axis];
        max = feature.max[axis];

        if (min >= k1 && max <= k2) { // trivial accept
            clipped.push(feature);
            continue;
        } else if (min > k2 || max < k1) continue; // trivial reject

        var slices = type === 1 ?
                clipPoints(geometry, k1, k2, axis) :
                clipGeometry(geometry, k1, k2, axis, intersect, type === 3);

        if (slices.length) {
            // if a feature got clipped, it will likely get clipped on the next zoom level as well,
            // so there's no need to recalculate bboxes
            clipped.push(createFeature(feature.tags, type, slices, feature.id));
        }
    }

    return clipped.length ? clipped : null;
}

function clipPoints(geometry, k1, k2, axis) {
    var slice = [];

    for (var i = 0; i < geometry.length; i++) {
        var a = geometry[i],
            ak = a[axis];

        if (ak >= k1 && ak <= k2) slice.push(a);
    }
    return slice;
}

function clipGeometry(geometry, k1, k2, axis, intersect, closed) {

    var slices = [];

    for (var i = 0; i < geometry.length; i++) {

        var ak = 0,
            bk = 0,
            b = null,
            points = geometry[i],
            area = points.area,
            dist = points.dist,
            outer = points.outer,
            len = points.length,
            a, j, last;

        var slice = [];

        for (j = 0; j < len - 1; j++) {
            a = b || points[j];
            b = points[j + 1];
            ak = bk || a[axis];
            bk = b[axis];

            if (ak < k1) {

                if ((bk > k2)) { // ---|-----|-->
                    slice.push(intersect(a, b, k1), intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist, outer);

                } else if (bk >= k1) slice.push(intersect(a, b, k1)); // ---|-->  |

            } else if (ak > k2) {

                if ((bk < k1)) { // <--|-----|---
                    slice.push(intersect(a, b, k2), intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist, outer);

                } else if (bk <= k2) slice.push(intersect(a, b, k2)); // |  <--|---

            } else {

                slice.push(a);

                if (bk < k1) { // <--|---  |
                    slice.push(intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist, outer);

                } else if (bk > k2) { // |  ---|-->
                    slice.push(intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist, outer);
                }
                // | --> |
            }
        }

        // add the last point
        a = points[len - 1];
        ak = a[axis];
        if (ak >= k1 && ak <= k2) slice.push(a);

        // close the polygon if its endpoints are not the same after clipping

        last = slice[slice.length - 1];
        if (closed && last && (slice[0][0] !== last[0] || slice[0][1] !== last[1])) slice.push(slice[0]);

        // add the final slice
        newSlice(slices, slice, area, dist, outer);
    }

    return slices;
}

function newSlice(slices, slice, area, dist, outer) {
    if (slice.length) {
        // we don't recalculate the area/length of the unclipped geometry because the case where it goes
        // below the visibility threshold as a result of clipping is rare, so we avoid doing unnecessary work
        slice.area = area;
        slice.dist = dist;
        if (outer !== undefined) slice.outer = outer;

        slices.push(slice);
    }
    return [];
}

},{"./feature":8}],7:[function(require,module,exports){
'use strict';

module.exports = convert;

var simplify = require('./simplify');
var createFeature = require('./feature');

// converts GeoJSON feature into an intermediate projected JSON vector format with simplification data

function convert(data, tolerance) {
    var features = [];

    if (data.type === 'FeatureCollection') {
        for (var i = 0; i < data.features.length; i++) {
            convertFeature(features, data.features[i], tolerance);
        }
    } else if (data.type === 'Feature') {
        convertFeature(features, data, tolerance);

    } else {
        // single geometry or a geometry collection
        convertFeature(features, {geometry: data}, tolerance);
    }
    return features;
}

function convertFeature(features, feature, tolerance) {
    if (feature.geometry === null) {
        // ignore features with null geometry
        return;
    }

    var geom = feature.geometry,
        type = geom.type,
        coords = geom.coordinates,
        tags = feature.properties,
        id = feature.id,
        i, j, rings, projectedRing;

    if (type === 'Point') {
        features.push(createFeature(tags, 1, [projectPoint(coords)], id));

    } else if (type === 'MultiPoint') {
        features.push(createFeature(tags, 1, project(coords), id));

    } else if (type === 'LineString') {
        features.push(createFeature(tags, 2, [project(coords, tolerance)], id));

    } else if (type === 'MultiLineString' || type === 'Polygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            projectedRing = project(coords[i], tolerance);
            if (type === 'Polygon') projectedRing.outer = (i === 0);
            rings.push(projectedRing);
        }
        features.push(createFeature(tags, type === 'Polygon' ? 3 : 2, rings, id));

    } else if (type === 'MultiPolygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            for (j = 0; j < coords[i].length; j++) {
                projectedRing = project(coords[i][j], tolerance);
                projectedRing.outer = (j === 0);
                rings.push(projectedRing);
            }
        }
        features.push(createFeature(tags, 3, rings, id));

    } else if (type === 'GeometryCollection') {
        for (i = 0; i < geom.geometries.length; i++) {
            convertFeature(features, {
                geometry: geom.geometries[i],
                properties: tags
            }, tolerance);
        }

    } else {
        throw new Error('Input data is not a valid GeoJSON object.');
    }
}

function project(lonlats, tolerance) {
    var projected = [];
    for (var i = 0; i < lonlats.length; i++) {
        projected.push(projectPoint(lonlats[i]));
    }
    if (tolerance) {
        simplify(projected, tolerance);
        calcSize(projected);
    }
    return projected;
}

function projectPoint(p) {
    var sin = Math.sin(p[1] * Math.PI / 180),
        x = (p[0] / 360 + 0.5),
        y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);

    y = y < 0 ? 0 :
        y > 1 ? 1 : y;

    return [x, y, 0];
}

// calculate area and length of the poly
function calcSize(points) {
    var area = 0,
        dist = 0;

    for (var i = 0, a, b; i < points.length - 1; i++) {
        a = b || points[i];
        b = points[i + 1];

        area += a[0] * b[1] - b[0] * a[1];

        // use Manhattan distance instead of Euclidian one to avoid expensive square root computation
        dist += Math.abs(b[0] - a[0]) + Math.abs(b[1] - a[1]);
    }
    points.area = Math.abs(area / 2);
    points.dist = dist;
}

},{"./feature":8,"./simplify":10}],8:[function(require,module,exports){
'use strict';

module.exports = createFeature;

function createFeature(tags, type, geom, id) {
    var feature = {
        id: id || null,
        type: type,
        geometry: geom,
        tags: tags || null,
        min: [Infinity, Infinity], // initial bbox values
        max: [-Infinity, -Infinity]
    };
    calcBBox(feature);
    return feature;
}

// calculate the feature bounding box for faster clipping later
function calcBBox(feature) {
    var geometry = feature.geometry,
        min = feature.min,
        max = feature.max;

    if (feature.type === 1) {
        calcRingBBox(min, max, geometry);
    } else {
        for (var i = 0; i < geometry.length; i++) {
            calcRingBBox(min, max, geometry[i]);
        }
    }

    return feature;
}

function calcRingBBox(min, max, points) {
    for (var i = 0, p; i < points.length; i++) {
        p = points[i];
        min[0] = Math.min(p[0], min[0]);
        max[0] = Math.max(p[0], max[0]);
        min[1] = Math.min(p[1], min[1]);
        max[1] = Math.max(p[1], max[1]);
    }
}

},{}],9:[function(require,module,exports){
'use strict';

module.exports = geojsonvt;

var convert = require('./convert'),     // GeoJSON conversion and preprocessing
    transform = require('./transform'), // coordinate transformation
    clip = require('./clip'),           // stripe clipping algorithm
    wrap = require('./wrap'),           // date line processing
    createTile = require('./tile');     // final simplified tile generation


function geojsonvt(data, options) {
    return new GeoJSONVT(data, options);
}

function GeoJSONVT(data, options) {
    options = this.options = extend(Object.create(this.options), options);

    var debug = options.debug;

    if (debug) console.time('preprocess data');

    var z2 = 1 << options.maxZoom, // 2^z
        features = convert(data, options.tolerance / (z2 * options.extent));

    this.tiles = {};
    this.tileCoords = [];

    if (debug) {
        console.timeEnd('preprocess data');
        console.log('index: maxZoom: %d, maxPoints: %d', options.indexMaxZoom, options.indexMaxPoints);
        console.time('generate tiles');
        this.stats = {};
        this.total = 0;
    }

    features = wrap(features, options.buffer / options.extent, intersectX);

    // start slicing from the top tile down
    if (features.length) this.splitTile(features, 0, 0, 0);

    if (debug) {
        if (features.length) console.log('features: %d, points: %d', this.tiles[0].numFeatures, this.tiles[0].numPoints);
        console.timeEnd('generate tiles');
        console.log('tiles generated:', this.total, JSON.stringify(this.stats));
    }
}

GeoJSONVT.prototype.options = {
    maxZoom: 14,            // max zoom to preserve detail on
    indexMaxZoom: 5,        // max zoom in the tile index
    indexMaxPoints: 100000, // max number of points per tile in the tile index
    solidChildren: false,   // whether to tile solid square tiles further
    tolerance: 3,           // simplification tolerance (higher means simpler)
    extent: 4096,           // tile extent
    buffer: 64,             // tile buffer on each side
    debug: 0                // logging level (0, 1 or 2)
};

GeoJSONVT.prototype.splitTile = function (features, z, x, y, cz, cx, cy) {

    var stack = [features, z, x, y],
        options = this.options,
        debug = options.debug,
        solid = null;

    // avoid recursion by using a processing queue
    while (stack.length) {
        y = stack.pop();
        x = stack.pop();
        z = stack.pop();
        features = stack.pop();

        var z2 = 1 << z,
            id = toID(z, x, y),
            tile = this.tiles[id],
            tileTolerance = z === options.maxZoom ? 0 : options.tolerance / (z2 * options.extent);

        if (!tile) {
            if (debug > 1) console.time('creation');

            tile = this.tiles[id] = createTile(features, z2, x, y, tileTolerance, z === options.maxZoom);
            this.tileCoords.push({z: z, x: x, y: y});

            if (debug) {
                if (debug > 1) {
                    console.log('tile z%d-%d-%d (features: %d, points: %d, simplified: %d)',
                        z, x, y, tile.numFeatures, tile.numPoints, tile.numSimplified);
                    console.timeEnd('creation');
                }
                var key = 'z' + z;
                this.stats[key] = (this.stats[key] || 0) + 1;
                this.total++;
            }
        }

        // save reference to original geometry in tile so that we can drill down later if we stop now
        tile.source = features;

        // if it's the first-pass tiling
        if (!cz) {
            // stop tiling if we reached max zoom, or if the tile is too simple
            if (z === options.indexMaxZoom || tile.numPoints <= options.indexMaxPoints) continue;

        // if a drilldown to a specific tile
        } else {
            // stop tiling if we reached base zoom or our target tile zoom
            if (z === options.maxZoom || z === cz) continue;

            // stop tiling if it's not an ancestor of the target tile
            var m = 1 << (cz - z);
            if (x !== Math.floor(cx / m) || y !== Math.floor(cy / m)) continue;
        }

        // stop tiling if the tile is solid clipped square
        if (!options.solidChildren && isClippedSquare(tile, options.extent, options.buffer)) {
            if (cz) solid = z; // and remember the zoom if we're drilling down
            continue;
        }

        // if we slice further down, no need to keep source geometry
        tile.source = null;

        if (debug > 1) console.time('clipping');

        // values we'll use for clipping
        var k1 = 0.5 * options.buffer / options.extent,
            k2 = 0.5 - k1,
            k3 = 0.5 + k1,
            k4 = 1 + k1,
            tl, bl, tr, br, left, right;

        tl = bl = tr = br = null;

        left  = clip(features, z2, x - k1, x + k3, 0, intersectX, tile.min[0], tile.max[0]);
        right = clip(features, z2, x + k2, x + k4, 0, intersectX, tile.min[0], tile.max[0]);

        if (left) {
            tl = clip(left, z2, y - k1, y + k3, 1, intersectY, tile.min[1], tile.max[1]);
            bl = clip(left, z2, y + k2, y + k4, 1, intersectY, tile.min[1], tile.max[1]);
        }

        if (right) {
            tr = clip(right, z2, y - k1, y + k3, 1, intersectY, tile.min[1], tile.max[1]);
            br = clip(right, z2, y + k2, y + k4, 1, intersectY, tile.min[1], tile.max[1]);
        }

        if (debug > 1) console.timeEnd('clipping');

        if (features.length) {
            stack.push(tl || [], z + 1, x * 2,     y * 2);
            stack.push(bl || [], z + 1, x * 2,     y * 2 + 1);
            stack.push(tr || [], z + 1, x * 2 + 1, y * 2);
            stack.push(br || [], z + 1, x * 2 + 1, y * 2 + 1);
        }
    }

    return solid;
};

GeoJSONVT.prototype.getTile = function (z, x, y) {
    var options = this.options,
        extent = options.extent,
        debug = options.debug;

    var z2 = 1 << z;
    x = ((x % z2) + z2) % z2; // wrap tile x coordinate

    var id = toID(z, x, y);
    if (this.tiles[id]) return transform.tile(this.tiles[id], extent);

    if (debug > 1) console.log('drilling down to z%d-%d-%d', z, x, y);

    var z0 = z,
        x0 = x,
        y0 = y,
        parent;

    while (!parent && z0 > 0) {
        z0--;
        x0 = Math.floor(x0 / 2);
        y0 = Math.floor(y0 / 2);
        parent = this.tiles[toID(z0, x0, y0)];
    }

    if (!parent || !parent.source) return null;

    // if we found a parent tile containing the original geometry, we can drill down from it
    if (debug > 1) console.log('found parent tile z%d-%d-%d', z0, x0, y0);

    // it parent tile is a solid clipped square, return it instead since it's identical
    if (isClippedSquare(parent, extent, options.buffer)) return transform.tile(parent, extent);

    if (debug > 1) console.time('drilling down');
    var solid = this.splitTile(parent.source, z0, x0, y0, z, x, y);
    if (debug > 1) console.timeEnd('drilling down');

    // one of the parent tiles was a solid clipped square
    if (solid !== null) {
        var m = 1 << (z - solid);
        id = toID(solid, Math.floor(x / m), Math.floor(y / m));
    }

    return this.tiles[id] ? transform.tile(this.tiles[id], extent) : null;
};

function toID(z, x, y) {
    return (((1 << z) * y + x) * 32) + z;
}

function intersectX(a, b, x) {
    return [x, (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0]) + a[1], 1];
}
function intersectY(a, b, y) {
    return [(y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) + a[0], y, 1];
}

function extend(dest, src) {
    for (var i in src) dest[i] = src[i];
    return dest;
}

// checks whether a tile is a whole-area fill after clipping; if it is, there's no sense slicing it further
function isClippedSquare(tile, extent, buffer) {

    var features = tile.source;
    if (features.length !== 1) return false;

    var feature = features[0];
    if (feature.type !== 3 || feature.geometry.length > 1) return false;

    var len = feature.geometry[0].length;
    if (len !== 5) return false;

    for (var i = 0; i < len; i++) {
        var p = transform.point(feature.geometry[0][i], extent, tile.z2, tile.x, tile.y);
        if ((p[0] !== -buffer && p[0] !== extent + buffer) ||
            (p[1] !== -buffer && p[1] !== extent + buffer)) return false;
    }

    return true;
}

},{"./clip":6,"./convert":7,"./tile":11,"./transform":12,"./wrap":13}],10:[function(require,module,exports){
'use strict';

module.exports = simplify;

// calculate simplification data using optimized Douglas-Peucker algorithm

function simplify(points, tolerance) {

    var sqTolerance = tolerance * tolerance,
        len = points.length,
        first = 0,
        last = len - 1,
        stack = [],
        i, maxSqDist, sqDist, index;

    // always retain the endpoints (1 is the max value)
    points[first][2] = 1;
    points[last][2] = 1;

    // avoid recursion by using a stack
    while (last) {

        maxSqDist = 0;

        for (i = first + 1; i < last; i++) {
            sqDist = getSqSegDist(points[i], points[first], points[last]);

            if (sqDist > maxSqDist) {
                index = i;
                maxSqDist = sqDist;
            }
        }

        if (maxSqDist > sqTolerance) {
            points[index][2] = maxSqDist; // save the point importance in squared pixels as a z coordinate
            stack.push(first);
            stack.push(index);
            first = index;

        } else {
            last = stack.pop();
            first = stack.pop();
        }
    }
}

// square distance from a point to a segment
function getSqSegDist(p, a, b) {

    var x = a[0], y = a[1],
        bx = b[0], by = b[1],
        px = p[0], py = p[1],
        dx = bx - x,
        dy = by - y;

    if (dx !== 0 || dy !== 0) {

        var t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy);

        if (t > 1) {
            x = bx;
            y = by;

        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = px - x;
    dy = py - y;

    return dx * dx + dy * dy;
}

},{}],11:[function(require,module,exports){
'use strict';

module.exports = createTile;

function createTile(features, z2, tx, ty, tolerance, noSimplify) {
    var tile = {
        features: [],
        numPoints: 0,
        numSimplified: 0,
        numFeatures: 0,
        source: null,
        x: tx,
        y: ty,
        z2: z2,
        transformed: false,
        min: [2, 1],
        max: [-1, 0]
    };
    for (var i = 0; i < features.length; i++) {
        tile.numFeatures++;
        addFeature(tile, features[i], tolerance, noSimplify);

        var min = features[i].min,
            max = features[i].max;

        if (min[0] < tile.min[0]) tile.min[0] = min[0];
        if (min[1] < tile.min[1]) tile.min[1] = min[1];
        if (max[0] > tile.max[0]) tile.max[0] = max[0];
        if (max[1] > tile.max[1]) tile.max[1] = max[1];
    }
    return tile;
}

function addFeature(tile, feature, tolerance, noSimplify) {

    var geom = feature.geometry,
        type = feature.type,
        simplified = [],
        sqTolerance = tolerance * tolerance,
        i, j, ring, p;

    if (type === 1) {
        for (i = 0; i < geom.length; i++) {
            simplified.push(geom[i]);
            tile.numPoints++;
            tile.numSimplified++;
        }

    } else {

        // simplify and transform projected coordinates for tile geometry
        for (i = 0; i < geom.length; i++) {
            ring = geom[i];

            // filter out tiny polylines & polygons
            if (!noSimplify && ((type === 2 && ring.dist < tolerance) ||
                                (type === 3 && ring.area < sqTolerance))) {
                tile.numPoints += ring.length;
                continue;
            }

            var simplifiedRing = [];

            for (j = 0; j < ring.length; j++) {
                p = ring[j];
                // keep points with importance > tolerance
                if (noSimplify || p[2] > sqTolerance) {
                    simplifiedRing.push(p);
                    tile.numSimplified++;
                }
                tile.numPoints++;
            }

            if (type === 3) rewind(simplifiedRing, ring.outer);

            simplified.push(simplifiedRing);
        }
    }

    if (simplified.length) {
        var tileFeature = {
            geometry: simplified,
            type: type,
            tags: feature.tags || null
        };
        if (feature.id !== null) {
            tileFeature.id = feature.id;
        }
        tile.features.push(tileFeature);
    }
}

function rewind(ring, clockwise) {
    var area = signedArea(ring);
    if (area < 0 === clockwise) ring.reverse();
}

function signedArea(ring) {
    var sum = 0;
    for (var i = 0, len = ring.length, j = len - 1, p1, p2; i < len; j = i++) {
        p1 = ring[i];
        p2 = ring[j];
        sum += (p2[0] - p1[0]) * (p1[1] + p2[1]);
    }
    return sum;
}

},{}],12:[function(require,module,exports){
'use strict';

exports.tile = transformTile;
exports.point = transformPoint;

// Transforms the coordinates of each feature in the given tile from
// mercator-projected space into (extent x extent) tile space.
function transformTile(tile, extent) {
    if (tile.transformed) return tile;

    var z2 = tile.z2,
        tx = tile.x,
        ty = tile.y,
        i, j, k;

    for (i = 0; i < tile.features.length; i++) {
        var feature = tile.features[i],
            geom = feature.geometry,
            type = feature.type;

        if (type === 1) {
            for (j = 0; j < geom.length; j++) geom[j] = transformPoint(geom[j], extent, z2, tx, ty);

        } else {
            for (j = 0; j < geom.length; j++) {
                var ring = geom[j];
                for (k = 0; k < ring.length; k++) ring[k] = transformPoint(ring[k], extent, z2, tx, ty);
            }
        }
    }

    tile.transformed = true;

    return tile;
}

function transformPoint(p, extent, z2, tx, ty) {
    var x = Math.round(extent * (p[0] * z2 - tx)),
        y = Math.round(extent * (p[1] * z2 - ty));
    return [x, y];
}

},{}],13:[function(require,module,exports){
'use strict';

var clip = require('./clip');
var createFeature = require('./feature');

module.exports = wrap;

function wrap(features, buffer, intersectX) {
    var merged = features,
        left  = clip(features, 1, -1 - buffer, buffer,     0, intersectX, -1, 2), // left world copy
        right = clip(features, 1,  1 - buffer, 2 + buffer, 0, intersectX, -1, 2); // right world copy

    if (left || right) {
        merged = clip(features, 1, -buffer, 1 + buffer, 0, intersectX, -1, 2) || []; // center world copy

        if (left) merged = shiftFeatureCoords(left, 1).concat(merged); // merge left into center
        if (right) merged = merged.concat(shiftFeatureCoords(right, -1)); // merge right into center
    }

    return merged;
}

function shiftFeatureCoords(features, offset) {
    var newFeatures = [];

    for (var i = 0; i < features.length; i++) {
        var feature = features[i],
            type = feature.type;

        var newGeometry;

        if (type === 1) {
            newGeometry = shiftCoords(feature.geometry, offset);
        } else {
            newGeometry = [];
            for (var j = 0; j < feature.geometry.length; j++) {
                newGeometry.push(shiftCoords(feature.geometry[j], offset));
            }
        }

        newFeatures.push(createFeature(feature.tags, type, newGeometry, feature.id));
    }

    return newFeatures;
}

function shiftCoords(points, offset) {
    var newPoints = [];
    newPoints.area = points.area;
    newPoints.dist = points.dist;

    for (var i = 0; i < points.length; i++) {
        newPoints.push([points[i][0] + offset, points[i][1], points[i][2]]);
    }
    return newPoints;
}

},{"./clip":6,"./feature":8}],14:[function(require,module,exports){
'use strict';

module.exports = GridIndex;

var NUM_PARAMS = 3;

function GridIndex(extent, n, padding) {
    var cells = this.cells = [];

    if (extent instanceof ArrayBuffer) {
        this.arrayBuffer = extent;
        var array = new Int32Array(this.arrayBuffer);
        extent = array[0];
        n = array[1];
        padding = array[2];

        this.d = n + 2 * padding;
        for (var k = 0; k < this.d * this.d; k++) {
            var start = array[NUM_PARAMS + k];
            var end = array[NUM_PARAMS + k + 1];
            cells.push(start === end ?
                    null :
                    array.subarray(start, end));
        }
        var keysOffset = array[NUM_PARAMS + cells.length];
        var bboxesOffset = array[NUM_PARAMS + cells.length + 1];
        this.keys = array.subarray(keysOffset, bboxesOffset);
        this.bboxes = array.subarray(bboxesOffset);

        this.insert = this._insertReadonly;

    } else {
        this.d = n + 2 * padding;
        for (var i = 0; i < this.d * this.d; i++) {
            cells.push([]);
        }
        this.keys = [];
        this.bboxes = [];
    }

    this.n = n;
    this.extent = extent;
    this.padding = padding;
    this.scale = n / extent;
    this.uid = 0;

    var p = (padding / n) * extent;
    this.min = -p;
    this.max = extent + p;
}


GridIndex.prototype.insert = function(key, x1, y1, x2, y2) {
    this._forEachCell(x1, y1, x2, y2, this._insertCell, this.uid++);
    this.keys.push(key);
    this.bboxes.push(x1);
    this.bboxes.push(y1);
    this.bboxes.push(x2);
    this.bboxes.push(y2);
};

GridIndex.prototype._insertReadonly = function() {
    throw 'Cannot insert into a GridIndex created from an ArrayBuffer.';
};

GridIndex.prototype._insertCell = function(x1, y1, x2, y2, cellIndex, uid) {
    this.cells[cellIndex].push(uid);
};

GridIndex.prototype.query = function(x1, y1, x2, y2) {
    var min = this.min;
    var max = this.max;
    if (x1 <= min && y1 <= min && max <= x2 && max <= y2) {
        // We use `Array#slice` because `this.keys` may be a `Int32Array` and
        // some browsers (Safari and IE) do not support `TypedArray#slice`
        // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/TypedArray/slice#Browser_compatibility
        return Array.prototype.slice.call(this.keys);

    } else {
        var result = [];
        var seenUids = {};
        this._forEachCell(x1, y1, x2, y2, this._queryCell, result, seenUids);
        return result;
    }
};

GridIndex.prototype._queryCell = function(x1, y1, x2, y2, cellIndex, result, seenUids) {
    var cell = this.cells[cellIndex];
    if (cell !== null) {
        var keys = this.keys;
        var bboxes = this.bboxes;
        for (var u = 0; u < cell.length; u++) {
            var uid = cell[u];
            if (seenUids[uid] === undefined) {
                var offset = uid * 4;
                if ((x1 <= bboxes[offset + 2]) &&
                    (y1 <= bboxes[offset + 3]) &&
                    (x2 >= bboxes[offset + 0]) &&
                    (y2 >= bboxes[offset + 1])) {
                    seenUids[uid] = true;
                    result.push(keys[uid]);
                } else {
                    seenUids[uid] = false;
                }
            }
        }
    }
};

GridIndex.prototype._forEachCell = function(x1, y1, x2, y2, fn, arg1, arg2) {
    var cx1 = this._convertToCellCoord(x1);
    var cy1 = this._convertToCellCoord(y1);
    var cx2 = this._convertToCellCoord(x2);
    var cy2 = this._convertToCellCoord(y2);
    for (var x = cx1; x <= cx2; x++) {
        for (var y = cy1; y <= cy2; y++) {
            var cellIndex = this.d * y + x;
            if (fn.call(this, x1, y1, x2, y2, cellIndex, arg1, arg2)) return;
        }
    }
};

GridIndex.prototype._convertToCellCoord = function(x) {
    return Math.max(0, Math.min(this.d - 1, Math.floor(x * this.scale) + this.padding));
};

GridIndex.prototype.toArrayBuffer = function() {
    if (this.arrayBuffer) return this.arrayBuffer;

    var cells = this.cells;

    var metadataLength = NUM_PARAMS + this.cells.length + 1 + 1;
    var totalCellLength = 0;
    for (var i = 0; i < this.cells.length; i++) {
        totalCellLength += this.cells[i].length;
    }

    var array = new Int32Array(metadataLength + totalCellLength + this.keys.length + this.bboxes.length);
    array[0] = this.extent;
    array[1] = this.n;
    array[2] = this.padding;

    var offset = metadataLength;
    for (var k = 0; k < cells.length; k++) {
        var cell = cells[k];
        array[NUM_PARAMS + k] = offset;
        array.set(cell, offset);
        offset += cell.length;
    }

    array[NUM_PARAMS + cells.length] = offset;
    array.set(this.keys, offset);
    offset += this.keys.length;

    array[NUM_PARAMS + cells.length + 1] = offset;
    array.set(this.bboxes, offset);
    offset += this.bboxes.length;

    return array.buffer;
};

},{}],15:[function(require,module,exports){
exports.read = function (buffer, offset, isLE, mLen, nBytes) {
  var e, m
  var eLen = nBytes * 8 - mLen - 1
  var eMax = (1 << eLen) - 1
  var eBias = eMax >> 1
  var nBits = -7
  var i = isLE ? (nBytes - 1) : 0
  var d = isLE ? -1 : 1
  var s = buffer[offset + i]

  i += d

  e = s & ((1 << (-nBits)) - 1)
  s >>= (-nBits)
  nBits += eLen
  for (; nBits > 0; e = e * 256 + buffer[offset + i], i += d, nBits -= 8) {}

  m = e & ((1 << (-nBits)) - 1)
  e >>= (-nBits)
  nBits += mLen
  for (; nBits > 0; m = m * 256 + buffer[offset + i], i += d, nBits -= 8) {}

  if (e === 0) {
    e = 1 - eBias
  } else if (e === eMax) {
    return m ? NaN : ((s ? -1 : 1) * Infinity)
  } else {
    m = m + Math.pow(2, mLen)
    e = e - eBias
  }
  return (s ? -1 : 1) * m * Math.pow(2, e - mLen)
}

exports.write = function (buffer, value, offset, isLE, mLen, nBytes) {
  var e, m, c
  var eLen = nBytes * 8 - mLen - 1
  var eMax = (1 << eLen) - 1
  var eBias = eMax >> 1
  var rt = (mLen === 23 ? Math.pow(2, -24) - Math.pow(2, -77) : 0)
  var i = isLE ? 0 : (nBytes - 1)
  var d = isLE ? 1 : -1
  var s = value < 0 || (value === 0 && 1 / value < 0) ? 1 : 0

  value = Math.abs(value)

  if (isNaN(value) || value === Infinity) {
    m = isNaN(value) ? 1 : 0
    e = eMax
  } else {
    e = Math.floor(Math.log(value) / Math.LN2)
    if (value * (c = Math.pow(2, -e)) < 1) {
      e--
      c *= 2
    }
    if (e + eBias >= 1) {
      value += rt / c
    } else {
      value += rt * Math.pow(2, 1 - eBias)
    }
    if (value * c >= 2) {
      e++
      c /= 2
    }

    if (e + eBias >= eMax) {
      m = 0
      e = eMax
    } else if (e + eBias >= 1) {
      m = (value * c - 1) * Math.pow(2, mLen)
      e = e + eBias
    } else {
      m = value * Math.pow(2, eBias - 1) * Math.pow(2, mLen)
      e = 0
    }
  }

  for (; mLen >= 8; buffer[offset + i] = m & 0xff, i += d, m /= 256, mLen -= 8) {}

  e = (e << mLen) | m
  eLen += mLen
  for (; eLen > 0; buffer[offset + i] = e & 0xff, i += d, e /= 256, eLen -= 8) {}

  buffer[offset + i - d] |= s * 128
}

},{}],16:[function(require,module,exports){
'use strict';

var sort = require('./sort');
var range = require('./range');
var within = require('./within');

module.exports = kdbush;

function kdbush(points, getX, getY, nodeSize, ArrayType) {
    return new KDBush(points, getX, getY, nodeSize, ArrayType);
}

function KDBush(points, getX, getY, nodeSize, ArrayType) {
    getX = getX || defaultGetX;
    getY = getY || defaultGetY;
    ArrayType = ArrayType || Array;

    this.nodeSize = nodeSize || 64;
    this.points = points;

    this.ids = new ArrayType(points.length);
    this.coords = new ArrayType(points.length * 2);

    for (var i = 0; i < points.length; i++) {
        this.ids[i] = i;
        this.coords[2 * i] = getX(points[i]);
        this.coords[2 * i + 1] = getY(points[i]);
    }

    sort(this.ids, this.coords, this.nodeSize, 0, this.ids.length - 1, 0);
}

KDBush.prototype = {
    range: function (minX, minY, maxX, maxY) {
        return range(this.ids, this.coords, minX, minY, maxX, maxY, this.nodeSize);
    },

    within: function (x, y, r) {
        return within(this.ids, this.coords, x, y, r, this.nodeSize);
    }
};

function defaultGetX(p) { return p[0]; }
function defaultGetY(p) { return p[1]; }

},{"./range":17,"./sort":18,"./within":19}],17:[function(require,module,exports){
'use strict';

module.exports = range;

function range(ids, coords, minX, minY, maxX, maxY, nodeSize) {
    var stack = [0, ids.length - 1, 0];
    var result = [];
    var x, y;

    while (stack.length) {
        var axis = stack.pop();
        var right = stack.pop();
        var left = stack.pop();

        if (right - left <= nodeSize) {
            for (var i = left; i <= right; i++) {
                x = coords[2 * i];
                y = coords[2 * i + 1];
                if (x >= minX && x <= maxX && y >= minY && y <= maxY) result.push(ids[i]);
            }
            continue;
        }

        var m = Math.floor((left + right) / 2);

        x = coords[2 * m];
        y = coords[2 * m + 1];

        if (x >= minX && x <= maxX && y >= minY && y <= maxY) result.push(ids[m]);

        var nextAxis = (axis + 1) % 2;

        if (axis === 0 ? minX <= x : minY <= y) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(nextAxis);
        }
        if (axis === 0 ? maxX >= x : maxY >= y) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(nextAxis);
        }
    }

    return result;
}

},{}],18:[function(require,module,exports){
'use strict';

module.exports = sortKD;

function sortKD(ids, coords, nodeSize, left, right, depth) {
    if (right - left <= nodeSize) return;

    var m = Math.floor((left + right) / 2);

    select(ids, coords, m, left, right, depth % 2);

    sortKD(ids, coords, nodeSize, left, m - 1, depth + 1);
    sortKD(ids, coords, nodeSize, m + 1, right, depth + 1);
}

function select(ids, coords, k, left, right, inc) {

    while (right > left) {
        if (right - left > 600) {
            var n = right - left + 1;
            var m = k - left + 1;
            var z = Math.log(n);
            var s = 0.5 * Math.exp(2 * z / 3);
            var sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);
            var newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
            var newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));
            select(ids, coords, k, newLeft, newRight, inc);
        }

        var t = coords[2 * k + inc];
        var i = left;
        var j = right;

        swapItem(ids, coords, left, k);
        if (coords[2 * right + inc] > t) swapItem(ids, coords, left, right);

        while (i < j) {
            swapItem(ids, coords, i, j);
            i++;
            j--;
            while (coords[2 * i + inc] < t) i++;
            while (coords[2 * j + inc] > t) j--;
        }

        if (coords[2 * left + inc] === t) swapItem(ids, coords, left, j);
        else {
            j++;
            swapItem(ids, coords, j, right);
        }

        if (j <= k) left = j + 1;
        if (k <= j) right = j - 1;
    }
}

function swapItem(ids, coords, i, j) {
    swap(ids, i, j);
    swap(coords, 2 * i, 2 * j);
    swap(coords, 2 * i + 1, 2 * j + 1);
}

function swap(arr, i, j) {
    var tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

},{}],19:[function(require,module,exports){
'use strict';

module.exports = within;

function within(ids, coords, qx, qy, r, nodeSize) {
    var stack = [0, ids.length - 1, 0];
    var result = [];
    var r2 = r * r;

    while (stack.length) {
        var axis = stack.pop();
        var right = stack.pop();
        var left = stack.pop();

        if (right - left <= nodeSize) {
            for (var i = left; i <= right; i++) {
                if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2) result.push(ids[i]);
            }
            continue;
        }

        var m = Math.floor((left + right) / 2);

        var x = coords[2 * m];
        var y = coords[2 * m + 1];

        if (sqDist(x, y, qx, qy) <= r2) result.push(ids[m]);

        var nextAxis = (axis + 1) % 2;

        if (axis === 0 ? qx - r <= x : qy - r <= y) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(nextAxis);
        }
        if (axis === 0 ? qx + r >= x : qy + r >= y) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(nextAxis);
        }
    }

    return result;
}

function sqDist(ax, ay, bx, by) {
    var dx = ax - bx;
    var dy = ay - by;
    return dx * dx + dy * dy;
}

},{}],20:[function(require,module,exports){
'use strict';
module.exports = 8192;
},{}],21:[function(require,module,exports){
'use strict';
var Point = require('point-geometry');
var loadGeometry = require('./load_geometry');
var EXTENT = require('./extent');
var featureFilter = require('../style-spec/feature_filter');
var createStructArrayType = require('../util/struct_array');
var Grid = require('grid-index');
var DictionaryCoder = require('../util/dictionary_coder');
var vt = require('vector-tile');
var Protobuf = require('pbf');
var GeoJSONFeature = require('../util/vectortile_to_geojson');
var arraysIntersect = require('../util/util').arraysIntersect;
var intersection = require('../util/intersection_tests');
var multiPolygonIntersectsBufferedMultiPoint = intersection.multiPolygonIntersectsBufferedMultiPoint;
var multiPolygonIntersectsMultiPolygon = intersection.multiPolygonIntersectsMultiPolygon;
var multiPolygonIntersectsBufferedMultiLine = intersection.multiPolygonIntersectsBufferedMultiLine;
var FeatureIndexArray = createStructArrayType({
    members: [
        {
            type: 'Uint32',
            name: 'featureIndex'
        },
        {
            type: 'Uint16',
            name: 'sourceLayerIndex'
        },
        {
            type: 'Uint16',
            name: 'bucketIndex'
        }
    ]
});
var FeatureIndex = function FeatureIndex(coord, overscaling, collisionTile) {
    if (coord.grid) {
        var serialized = coord;
        var rawTileData = overscaling;
        coord = serialized.coord;
        overscaling = serialized.overscaling;
        this.grid = new Grid(serialized.grid);
        this.featureIndexArray = new FeatureIndexArray(serialized.featureIndexArray);
        this.rawTileData = rawTileData;
        this.bucketLayerIDs = serialized.bucketLayerIDs;
        this.paintPropertyStatistics = serialized.paintPropertyStatistics;
    } else {
        this.grid = new Grid(EXTENT, 16, 0);
        this.featureIndexArray = new FeatureIndexArray();
    }
    this.coord = coord;
    this.overscaling = overscaling;
    this.x = coord.x;
    this.y = coord.y;
    this.z = coord.z - Math.log(overscaling) / Math.LN2;
    this.setCollisionTile(collisionTile);
};
FeatureIndex.prototype.insert = function insert(feature, bucketIndex) {
    var this$1 = this;
    var key = this.featureIndexArray.length;
    this.featureIndexArray.emplaceBack(feature.index, feature.sourceLayerIndex, bucketIndex);
    var geometry = loadGeometry(feature);
    for (var r = 0; r < geometry.length; r++) {
        var ring = geometry[r];
        var bbox = [
            Infinity,
            Infinity,
            -Infinity,
            -Infinity
        ];
        for (var i = 0; i < ring.length; i++) {
            var p = ring[i];
            bbox[0] = Math.min(bbox[0], p.x);
            bbox[1] = Math.min(bbox[1], p.y);
            bbox[2] = Math.max(bbox[2], p.x);
            bbox[3] = Math.max(bbox[3], p.y);
        }
        this$1.grid.insert(key, bbox[0], bbox[1], bbox[2], bbox[3]);
    }
};
FeatureIndex.prototype.setCollisionTile = function setCollisionTile(collisionTile) {
    this.collisionTile = collisionTile;
};
FeatureIndex.prototype.serialize = function serialize(transferables) {
    var grid = this.grid.toArrayBuffer();
    if (transferables) {
        transferables.push(grid);
    }
    return {
        coord: this.coord,
        overscaling: this.overscaling,
        grid: grid,
        featureIndexArray: this.featureIndexArray.serialize(transferables),
        bucketLayerIDs: this.bucketLayerIDs,
        paintPropertyStatistics: this.paintPropertyStatistics
    };
};
FeatureIndex.prototype.query = function query(args, styleLayers) {
    var this$1 = this;
    if (!this.vtLayers) {
        this.vtLayers = new vt.VectorTile(new Protobuf(this.rawTileData)).layers;
        this.sourceLayerCoder = new DictionaryCoder(this.vtLayers ? Object.keys(this.vtLayers).sort() : ['_geojsonTileLayer']);
    }
    var result = {};
    var params = args.params || {}, pixelsToTileUnits = EXTENT / args.tileSize / args.scale, filter = featureFilter(params.filter);
    var additionalRadius = 0;
    for (var id in styleLayers) {
        if (!this$1.hasLayer(id)) {
            continue;
        }
        var styleLayer = styleLayers[id];
        var styleLayerDistance = 0;
        if (styleLayer.type === 'line') {
            var width = getLineWidth(this$1.getPaintValue('line-width', styleLayer), this$1.getPaintValue('line-gap-width', styleLayer));
            var offset = this$1.getPaintValue('line-offset', styleLayer);
            var translate = this$1.getPaintValue('line-translate', styleLayer);
            styleLayerDistance = width / 2 + Math.abs(offset) + translateDistance(translate);
        } else if (styleLayer.type === 'fill') {
            styleLayerDistance = translateDistance(this$1.getPaintValue('fill-translate', styleLayer));
        } else if (styleLayer.type === 'fill-extrusion') {
            styleLayerDistance = translateDistance(this$1.getPaintValue('fill-extrusion-translate', styleLayer));
        } else if (styleLayer.type === 'circle') {
            styleLayerDistance = this$1.getPaintValue('circle-radius', styleLayer) + translateDistance(this$1.getPaintValue('circle-translate', styleLayer));
        }
        additionalRadius = Math.max(additionalRadius, styleLayerDistance * pixelsToTileUnits);
    }
    var queryGeometry = args.queryGeometry.map(function (q) {
        return q.map(function (p) {
            return new Point(p.x, p.y);
        });
    });
    var minX = Infinity;
    var minY = Infinity;
    var maxX = -Infinity;
    var maxY = -Infinity;
    for (var i = 0; i < queryGeometry.length; i++) {
        var ring = queryGeometry[i];
        for (var k = 0; k < ring.length; k++) {
            var p = ring[k];
            minX = Math.min(minX, p.x);
            minY = Math.min(minY, p.y);
            maxX = Math.max(maxX, p.x);
            maxY = Math.max(maxY, p.y);
        }
    }
    var matching = this.grid.query(minX - additionalRadius, minY - additionalRadius, maxX + additionalRadius, maxY + additionalRadius);
    matching.sort(topDownFeatureComparator);
    this.filterMatching(result, matching, this.featureIndexArray, queryGeometry, filter, params.layers, styleLayers, args.bearing, pixelsToTileUnits);
    var matchingSymbols = this.collisionTile.queryRenderedSymbols(queryGeometry, args.scale);
    matchingSymbols.sort();
    this.filterMatching(result, matchingSymbols, this.collisionTile.collisionBoxArray, queryGeometry, filter, params.layers, styleLayers, args.bearing, pixelsToTileUnits);
    return result;
};
FeatureIndex.prototype.filterMatching = function filterMatching(result, matching, array, queryGeometry, filter, filterLayerIDs, styleLayers, bearing, pixelsToTileUnits) {
    var this$1 = this;
    var previousIndex;
    for (var k = 0; k < matching.length; k++) {
        var index = matching[k];
        if (index === previousIndex) {
            continue;
        }
        previousIndex = index;
        var match = array.get(index);
        var layerIDs = this$1.bucketLayerIDs[match.bucketIndex];
        if (filterLayerIDs && !arraysIntersect(filterLayerIDs, layerIDs)) {
            continue;
        }
        var sourceLayerName = this$1.sourceLayerCoder.decode(match.sourceLayerIndex);
        var sourceLayer = this$1.vtLayers[sourceLayerName];
        var feature = sourceLayer.feature(match.featureIndex);
        if (!filter(feature)) {
            continue;
        }
        var geometry = null;
        for (var l = 0; l < layerIDs.length; l++) {
            var layerID = layerIDs[l];
            if (filterLayerIDs && filterLayerIDs.indexOf(layerID) < 0) {
                continue;
            }
            var styleLayer = styleLayers[layerID];
            if (!styleLayer) {
                continue;
            }
            var translatedPolygon = void 0;
            if (styleLayer.type !== 'symbol') {
                if (!geometry) {
                    geometry = loadGeometry(feature);
                }
                if (styleLayer.type === 'line') {
                    translatedPolygon = translate(queryGeometry, this$1.getPaintValue('line-translate', styleLayer, feature), this$1.getPaintValue('line-translate-anchor', styleLayer, feature), bearing, pixelsToTileUnits);
                    var halfWidth = pixelsToTileUnits / 2 * getLineWidth(this$1.getPaintValue('line-width', styleLayer, feature), this$1.getPaintValue('line-gap-width', styleLayer, feature));
                    var lineOffset = this$1.getPaintValue('line-offset', styleLayer, feature);
                    if (lineOffset) {
                        geometry = offsetLine(geometry, lineOffset * pixelsToTileUnits);
                    }
                    if (!multiPolygonIntersectsBufferedMultiLine(translatedPolygon, geometry, halfWidth)) {
                        continue;
                    }
                } else if (styleLayer.type === 'fill' || styleLayer.type === 'fill-extrusion') {
                    var typePrefix = styleLayer.type;
                    translatedPolygon = translate(queryGeometry, this$1.getPaintValue(typePrefix + '-translate', styleLayer, feature), this$1.getPaintValue(typePrefix + '-translate-anchor', styleLayer, feature), bearing, pixelsToTileUnits);
                    if (!multiPolygonIntersectsMultiPolygon(translatedPolygon, geometry)) {
                        continue;
                    }
                } else if (styleLayer.type === 'circle') {
                    translatedPolygon = translate(queryGeometry, this$1.getPaintValue('circle-translate', styleLayer, feature), this$1.getPaintValue('circle-translate-anchor', styleLayer, feature), bearing, pixelsToTileUnits);
                    var circleRadius = this$1.getPaintValue('circle-radius', styleLayer, feature) * pixelsToTileUnits;
                    if (!multiPolygonIntersectsBufferedMultiPoint(translatedPolygon, geometry, circleRadius)) {
                        continue;
                    }
                }
            }
            var geojsonFeature = new GeoJSONFeature(feature, this$1.z, this$1.x, this$1.y);
            geojsonFeature.layer = styleLayer.serialize();
            var layerResult = result[layerID];
            if (layerResult === undefined) {
                layerResult = result[layerID] = [];
            }
            layerResult.push(geojsonFeature);
        }
    }
};
FeatureIndex.prototype.hasLayer = function hasLayer(id) {
    var this$1 = this;
    for (var index in this$1.bucketLayerIDs) {
        for (var i = 0, list = this$1.bucketLayerIDs[index]; i < list.length; i += 1) {
            var layerID = list[i];
            if (id === layerID) {
                return true;
            }
        }
    }
    return false;
};
FeatureIndex.prototype.getPaintValue = function getPaintValue(property, layer, feature) {
    var featureConstant = layer.isPaintValueFeatureConstant(property);
    if (featureConstant || feature) {
        var featureProperties = feature ? feature.properties : {};
        return layer.getPaintValue(property, { zoom: this.z }, featureProperties);
    }
    return this.paintPropertyStatistics[layer.id][property].max;
};
module.exports = FeatureIndex;
function translateDistance(translate) {
    return Math.sqrt(translate[0] * translate[0] + translate[1] * translate[1]);
}
function topDownFeatureComparator(a, b) {
    return b - a;
}
function getLineWidth(lineWidth, lineGapWidth) {
    if (lineGapWidth > 0) {
        return lineGapWidth + 2 * lineWidth;
    } else {
        return lineWidth;
    }
}
function translate(queryGeometry, translate, translateAnchor, bearing, pixelsToTileUnits) {
    if (!translate[0] && !translate[1]) {
        return queryGeometry;
    }
    translate = Point.convert(translate);
    if (translateAnchor === 'viewport') {
        translate._rotate(-bearing);
    }
    var translated = [];
    for (var i = 0; i < queryGeometry.length; i++) {
        var ring = queryGeometry[i];
        var translatedRing = [];
        for (var k = 0; k < ring.length; k++) {
            translatedRing.push(ring[k].sub(translate._mult(pixelsToTileUnits)));
        }
        translated.push(translatedRing);
    }
    return translated;
}
function offsetLine(rings, offset) {
    var newRings = [];
    var zero = new Point(0, 0);
    for (var k = 0; k < rings.length; k++) {
        var ring = rings[k];
        var newRing = [];
        for (var i = 0; i < ring.length; i++) {
            var a = ring[i - 1];
            var b = ring[i];
            var c = ring[i + 1];
            var aToB = i === 0 ? zero : b.sub(a)._unit()._perp();
            var bToC = i === ring.length - 1 ? zero : c.sub(b)._unit()._perp();
            var extrude = aToB._add(bToC)._unit();
            var cosHalfAngle = extrude.x * bToC.x + extrude.y * bToC.y;
            extrude._mult(1 / cosHalfAngle);
            newRing.push(extrude._mult(offset)._add(b));
        }
        newRings.push(newRing);
    }
    return newRings;
}
},{"../style-spec/feature_filter":29,"../util/dictionary_coder":34,"../util/intersection_tests":36,"../util/struct_array":37,"../util/util":38,"../util/vectortile_to_geojson":39,"./extent":20,"./load_geometry":22,"grid-index":14,"pbf":41,"point-geometry":42,"vector-tile":44}],22:[function(require,module,exports){
'use strict';
var util = require('../util/util');
var EXTENT = require('./extent');
function createBounds(bits) {
    return {
        min: -1 * Math.pow(2, bits - 1),
        max: Math.pow(2, bits - 1) - 1
    };
}
var boundsLookup = {
    15: createBounds(15),
    16: createBounds(16)
};
module.exports = function loadGeometry(feature, bits) {
    var bounds = boundsLookup[bits || 16];
    var scale = EXTENT / feature.extent;
    var geometry = feature.loadGeometry();
    for (var r = 0; r < geometry.length; r++) {
        var ring = geometry[r];
        for (var p = 0; p < ring.length; p++) {
            var point = ring[p];
            point.x = Math.round(point.x * scale);
            point.y = Math.round(point.y * scale);
            if (point.x < bounds.min || point.x > bounds.max || point.y < bounds.min || point.y > bounds.max) {
                util.warnOnce('Geometry exceeds allowed extent, reduce your vector tile buffer size');
            }
        }
    }
    return geometry;
};
},{"../util/util":38,"./extent":20}],23:[function(require,module,exports){
'use strict';
var Coordinate = function Coordinate(column, row, zoom) {
    this.column = column;
    this.row = row;
    this.zoom = zoom;
};
Coordinate.prototype.clone = function clone() {
    return new Coordinate(this.column, this.row, this.zoom);
};
Coordinate.prototype.zoomTo = function zoomTo(zoom) {
    return this.clone()._zoomTo(zoom);
};
Coordinate.prototype.sub = function sub(c) {
    return this.clone()._sub(c);
};
Coordinate.prototype._zoomTo = function _zoomTo(zoom) {
    var scale = Math.pow(2, zoom - this.zoom);
    this.column *= scale;
    this.row *= scale;
    this.zoom = zoom;
    return this;
};
Coordinate.prototype._sub = function _sub(c) {
    c = c.zoomTo(this.zoom);
    this.column -= c.column;
    this.row -= c.row;
    return this;
};
module.exports = Coordinate;
},{}],24:[function(require,module,exports){
'use strict';
var Evented = require('../util/evented');
var util = require('../util/util');
var window = require('../util/window');
var EXTENT = require('../data/extent');
var GeoJSONSource = function (Evented) {
    function GeoJSONSource(id, options, dispatcher, eventedParent) {
        Evented.call(this);
        options = options || {};
        this.id = id;
        this.type = 'geojson';
        this.minzoom = 0;
        this.maxzoom = 18;
        this.tileSize = 512;
        this.isTileClipped = true;
        this.reparseOverscaled = true;
        this.dispatcher = dispatcher;
        this.setEventedParent(eventedParent);
        this._data = options.data;
        if (options.maxzoom !== undefined) {
            this.maxzoom = options.maxzoom;
        }
        if (options.type) {
            this.type = options.type;
        }
        var scale = EXTENT / this.tileSize;
        this.workerOptions = util.extend({
            source: this.id,
            cluster: options.cluster || false,
            geojsonVtOptions: {
                buffer: (options.buffer !== undefined ? options.buffer : 128) * scale,
                tolerance: (options.tolerance !== undefined ? options.tolerance : 0.375) * scale,
                extent: EXTENT,
                maxZoom: this.maxzoom
            },
            superclusterOptions: {
                maxZoom: Math.min(options.clusterMaxZoom, this.maxzoom - 1) || this.maxzoom - 1,
                extent: EXTENT,
                radius: (options.clusterRadius || 50) * scale,
                log: false
            }
        }, options.workerOptions);
    }
    if (Evented)
        GeoJSONSource.__proto__ = Evented;
    GeoJSONSource.prototype = Object.create(Evented && Evented.prototype);
    GeoJSONSource.prototype.constructor = GeoJSONSource;
    GeoJSONSource.prototype.load = function load() {
        var this$1 = this;
        this.fire('dataloading', { dataType: 'source' });
        this._updateWorkerData(function (err) {
            if (err) {
                this$1.fire('error', { error: err });
                return;
            }
            this$1.fire('data', {
                dataType: 'source',
                sourceDataType: 'metadata'
            });
        });
    };
    GeoJSONSource.prototype.onAdd = function onAdd(map) {
        this.load();
        this.map = map;
    };
    GeoJSONSource.prototype.setData = function setData(data) {
        var this$1 = this;
        this._data = data;
        this.fire('dataloading', { dataType: 'source' });
        this._updateWorkerData(function (err) {
            if (err) {
                return this$1.fire('error', { error: err });
            }
            this$1.fire('data', {
                dataType: 'source',
                sourceDataType: 'content'
            });
        });
        return this;
    };
    GeoJSONSource.prototype._updateWorkerData = function _updateWorkerData(callback) {
        var this$1 = this;
        var options = util.extend({}, this.workerOptions);
        var data = this._data;
        if (typeof data === 'string') {
            options.url = resolveURL(data);
        } else {
            options.data = JSON.stringify(data);
        }
        this.workerID = this.dispatcher.send(this.type + '.loadData', options, function (err) {
            this$1._loaded = true;
            callback(err);
        });
    };
    GeoJSONSource.prototype.loadTile = function loadTile(tile, callback) {
        var this$1 = this;
        var overscaling = tile.coord.z > this.maxzoom ? Math.pow(2, tile.coord.z - this.maxzoom) : 1;
        var params = {
            type: this.type,
            uid: tile.uid,
            coord: tile.coord,
            zoom: tile.coord.z,
            maxZoom: this.maxzoom,
            tileSize: this.tileSize,
            source: this.id,
            overscaling: overscaling,
            angle: this.map.transform.angle,
            pitch: this.map.transform.pitch,
            showCollisionBoxes: this.map.showCollisionBoxes
        };
        tile.workerID = this.dispatcher.send('loadTile', params, function (err, data) {
            tile.unloadVectorData();
            if (tile.aborted) {
                return;
            }
            if (err) {
                return callback(err);
            }
            tile.loadVectorData(data, this$1.map.painter);
            if (tile.redoWhenDone) {
                tile.redoWhenDone = false;
                tile.redoPlacement(this$1);
            }
            return callback(null);
        }, this.workerID);
    };
    GeoJSONSource.prototype.abortTile = function abortTile(tile) {
        tile.aborted = true;
    };
    GeoJSONSource.prototype.unloadTile = function unloadTile(tile) {
        tile.unloadVectorData();
        this.dispatcher.send('removeTile', {
            uid: tile.uid,
            type: this.type,
            source: this.id
        }, function () {
        }, tile.workerID);
    };
    GeoJSONSource.prototype.onRemove = function onRemove() {
        this.dispatcher.broadcast('removeSource', {
            type: this.type,
            source: this.id
        }, function () {
        });
    };
    GeoJSONSource.prototype.serialize = function serialize() {
        return {
            type: this.type,
            data: this._data
        };
    };
    return GeoJSONSource;
}(Evented);
function resolveURL(url) {
    var a = window.document.createElement('a');
    a.href = url;
    return a.href;
}
module.exports = GeoJSONSource;
},{"../data/extent":20,"../util/evented":35,"../util/util":38,"../util/window":33}],25:[function(require,module,exports){
'use strict';
var ajax = require('../util/ajax');
var rewind = require('geojson-rewind');
var GeoJSONWrapper = require('./geojson_wrapper');
var vtpbf = require('vt-pbf');
var supercluster = require('supercluster');
var geojsonvt = require('geojson-vt');
var VectorTileWorkerSource = require('./vector_tile_worker_source');
var GeoJSONWorkerSource = function (VectorTileWorkerSource) {
    function GeoJSONWorkerSource(actor, layerIndex, loadGeoJSON) {
        VectorTileWorkerSource.call(this, actor, layerIndex);
        if (loadGeoJSON) {
            this.loadGeoJSON = loadGeoJSON;
        }
        this._geoJSONIndexes = {};
    }
    if (VectorTileWorkerSource)
        GeoJSONWorkerSource.__proto__ = VectorTileWorkerSource;
    GeoJSONWorkerSource.prototype = Object.create(VectorTileWorkerSource && VectorTileWorkerSource.prototype);
    GeoJSONWorkerSource.prototype.constructor = GeoJSONWorkerSource;
    GeoJSONWorkerSource.prototype.loadVectorData = function loadVectorData(params, callback) {
        var source = params.source, coord = params.coord;
        if (!this._geoJSONIndexes[source]) {
            return callback(null, null);
        }
        var geoJSONTile = this._geoJSONIndexes[source].getTile(Math.min(coord.z, params.maxZoom), coord.x, coord.y);
        if (!geoJSONTile) {
            return callback(null, null);
        }
        var geojsonWrapper = new GeoJSONWrapper(geoJSONTile.features);
        geojsonWrapper.name = '_geojsonTileLayer';
        var pbf = vtpbf({ layers: { '_geojsonTileLayer': geojsonWrapper } });
        if (pbf.byteOffset !== 0 || pbf.byteLength !== pbf.buffer.byteLength) {
            pbf = new Uint8Array(pbf);
        }
        geojsonWrapper.rawData = pbf.buffer;
        callback(null, geojsonWrapper);
    };
    GeoJSONWorkerSource.prototype.loadData = function loadData(params, callback) {
        var handleData = function (err, data) {
            var this$1 = this;
            if (err) {
                return callback(err);
            }
            if (typeof data != 'object') {
                return callback(new Error('Input data is not a valid GeoJSON object.'));
            }
            rewind(data, true);
            this._indexData(data, params, function (err, indexed) {
                if (err) {
                    return callback(err);
                }
                this$1._geoJSONIndexes[params.source] = indexed;
                callback(null);
            });
        }.bind(this);
        this.loadGeoJSON(params, handleData);
    };
    GeoJSONWorkerSource.prototype.loadGeoJSON = function loadGeoJSON(params, callback) {
        if (params.url) {
            ajax.getJSON(params.url, callback);
        } else if (typeof params.data === 'string') {
            try {
                return callback(null, JSON.parse(params.data));
            } catch (e) {
                return callback(new Error('Input data is not a valid GeoJSON object.'));
            }
        } else {
            return callback(new Error('Input data is not a valid GeoJSON object.'));
        }
    };
    GeoJSONWorkerSource.prototype.removeSource = function removeSource(params) {
        if (this._geoJSONIndexes[params.source]) {
            delete this._geoJSONIndexes[params.source];
        }
    };
    GeoJSONWorkerSource.prototype._indexData = function _indexData(data, params, callback) {
        try {
            if (params.cluster) {
                callback(null, supercluster(params.superclusterOptions).load(data.features));
            } else {
                callback(null, geojsonvt(data, params.geojsonVtOptions));
            }
        } catch (err) {
            return callback(err);
        }
    };
    return GeoJSONWorkerSource;
}(VectorTileWorkerSource);
module.exports = GeoJSONWorkerSource;
},{"../util/ajax":32,"./geojson_wrapper":26,"./vector_tile_worker_source":27,"geojson-rewind":5,"geojson-vt":9,"supercluster":43,"vt-pbf":48}],26:[function(require,module,exports){
'use strict';
var Point = require('point-geometry');
var VectorTileFeature = require('vector-tile').VectorTileFeature;
var EXTENT = require('../data/extent');
var FeatureWrapper = function FeatureWrapper(feature) {
    var this$1 = this;
    this.type = feature.type;
    if (feature.type === 1) {
        this.rawGeometry = [];
        for (var i = 0; i < feature.geometry.length; i++) {
            this$1.rawGeometry.push([feature.geometry[i]]);
        }
    } else {
        this.rawGeometry = feature.geometry;
    }
    this.properties = feature.tags;
    if ('id' in feature && !isNaN(feature.id)) {
        this.id = parseInt(feature.id, 10);
    }
    this.extent = EXTENT;
};
FeatureWrapper.prototype.loadGeometry = function loadGeometry() {
    var this$1 = this;
    var rings = this.rawGeometry;
    this.geometry = [];
    for (var i = 0; i < rings.length; i++) {
        var ring = rings[i], newRing = [];
        for (var j = 0; j < ring.length; j++) {
            newRing.push(new Point(ring[j][0], ring[j][1]));
        }
        this$1.geometry.push(newRing);
    }
    return this.geometry;
};
FeatureWrapper.prototype.bbox = function bbox() {
    if (!this.geometry) {
        this.loadGeometry();
    }
    var rings = this.geometry;
    var x1 = Infinity, x2 = -Infinity, y1 = Infinity, y2 = -Infinity;
    for (var i = 0; i < rings.length; i++) {
        var ring = rings[i];
        for (var j = 0; j < ring.length; j++) {
            var coord = ring[j];
            x1 = Math.min(x1, coord.x);
            x2 = Math.max(x2, coord.x);
            y1 = Math.min(y1, coord.y);
            y2 = Math.max(y2, coord.y);
        }
    }
    return [
        x1,
        y1,
        x2,
        y2
    ];
};
FeatureWrapper.prototype.toGeoJSON = function toGeoJSON() {
    VectorTileFeature.prototype.toGeoJSON.call(this);
};
var GeoJSONWrapper = function GeoJSONWrapper(features) {
    this.features = features;
    this.length = features.length;
    this.extent = EXTENT;
};
GeoJSONWrapper.prototype.feature = function feature(i) {
    return new FeatureWrapper(this.features[i]);
};
module.exports = GeoJSONWrapper;
},{"../data/extent":20,"point-geometry":42,"vector-tile":44}],27:[function(require,module,exports){
'use strict';
var ajax = require('../util/ajax');
var vt = require('vector-tile');
var Protobuf = require('pbf');
var WorkerTile = require('./worker_tile');
var util = require('../util/util');
var VectorTileWorkerSource = function VectorTileWorkerSource(actor, layerIndex, loadVectorData) {
    this.actor = actor;
    this.layerIndex = layerIndex;
    if (loadVectorData) {
        this.loadVectorData = loadVectorData;
    }
    this.loading = {};
    this.loaded = {};
};
VectorTileWorkerSource.prototype.loadTile = function loadTile(params, callback) {
    var source = params.source, uid = params.uid;
    if (!this.loading[source]) {
        this.loading[source] = {};
    }
    var workerTile = this.loading[source][uid] = new WorkerTile(params);
    workerTile.abort = this.loadVectorData(params, done.bind(this));
    function done(err, vectorTile) {
        delete this.loading[source][uid];
        if (err) {
            return callback(err);
        }
        if (!vectorTile) {
            return callback(null, null);
        }
        workerTile.vectorTile = vectorTile;
        workerTile.parse(vectorTile, this.layerIndex, this.actor, function (err, result, transferrables) {
            if (err) {
                return callback(err);
            }
            var cacheControl = {};
            if (vectorTile.expires) {
                cacheControl.expires = vectorTile.expires;
            }
            if (vectorTile.cacheControl) {
                cacheControl.cacheControl = vectorTile.cacheControl;
            }
            callback(null, util.extend({ rawTileData: vectorTile.rawData }, result, cacheControl), transferrables);
        });
        this.loaded[source] = this.loaded[source] || {};
        this.loaded[source][uid] = workerTile;
    }
};
VectorTileWorkerSource.prototype.reloadTile = function reloadTile(params, callback) {
    var loaded = this.loaded[params.source], uid = params.uid, vtSource = this;
    if (loaded && loaded[uid]) {
        var workerTile = loaded[uid];
        if (workerTile.status === 'parsing') {
            workerTile.reloadCallback = callback;
        } else if (workerTile.status === 'done') {
            workerTile.parse(workerTile.vectorTile, this.layerIndex, this.actor, done.bind(workerTile));
        }
    }
    function done(err, data) {
        if (this.reloadCallback) {
            var reloadCallback = this.reloadCallback;
            delete this.reloadCallback;
            this.parse(this.vectorTile, vtSource.layerIndex, vtSource.actor, reloadCallback);
        }
        callback(err, data);
    }
};
VectorTileWorkerSource.prototype.abortTile = function abortTile(params) {
    var loading = this.loading[params.source], uid = params.uid;
    if (loading && loading[uid] && loading[uid].abort) {
        loading[uid].abort();
        delete loading[uid];
    }
};
VectorTileWorkerSource.prototype.removeTile = function removeTile(params) {
    var loaded = this.loaded[params.source], uid = params.uid;
    if (loaded && loaded[uid]) {
        delete loaded[uid];
    }
};
VectorTileWorkerSource.prototype.loadVectorData = function loadVectorData(params, callback) {
    var xhr = ajax.getArrayBuffer(params.url, done.bind(this));
    return function abort() {
        xhr.abort();
    };
    function done(err, response) {
        if (err) {
            return callback(err);
        }
        var vectorTile = new vt.VectorTile(new Protobuf(response.data));
        vectorTile.rawData = response.data;
        vectorTile.cacheControl = response.cacheControl;
        vectorTile.expires = response.expires;
        callback(err, vectorTile);
    }
};
VectorTileWorkerSource.prototype.redoPlacement = function redoPlacement(params, callback) {
    var loaded = this.loaded[params.source], loading = this.loading[params.source], uid = params.uid;
    if (loaded && loaded[uid]) {
        var workerTile = loaded[uid];
        var result = workerTile.redoPlacement(params.angle, params.pitch, params.showCollisionBoxes);
        if (result.result) {
            callback(null, result.result, result.transferables);
        }
    } else if (loading && loading[uid]) {
        loading[uid].angle = params.angle;
    }
};
module.exports = VectorTileWorkerSource;
},{"../util/ajax":32,"../util/util":38,"./worker_tile":28,"pbf":41,"vector-tile":44}],28:[function(require,module,exports){
'use strict';
var FeatureIndex = require('../data/feature_index');
var CollisionTile = require('../symbol/collision_tile');
var CollisionBoxArray = require('../symbol/collision_box');
var DictionaryCoder = require('../util/dictionary_coder');
var util = require('../util/util');
var WorkerTile = function WorkerTile(params) {
    this.coord = params.coord;
    this.uid = params.uid;
    this.zoom = params.zoom;
    this.tileSize = params.tileSize;
    this.source = params.source;
    this.overscaling = params.overscaling;
    this.angle = params.angle;
    this.pitch = params.pitch;
    this.showCollisionBoxes = params.showCollisionBoxes;
};
WorkerTile.prototype.parse = function parse(data, layerIndex, actor, callback) {
    var this$1 = this;
    if (!data.layers) {
        data = { layers: { '_geojsonTileLayer': data } };
    }
    this.status = 'parsing';
    this.data = data;
    this.collisionBoxArray = new CollisionBoxArray();
    var sourceLayerCoder = new DictionaryCoder(Object.keys(data.layers).sort());
    var featureIndex = new FeatureIndex(this.coord, this.overscaling);
    featureIndex.bucketLayerIDs = {};
    var buckets = {};
    var bucketIndex = 0;
    var options = {
        featureIndex: featureIndex,
        iconDependencies: {},
        glyphDependencies: {}
    };
    var layerFamilies = layerIndex.familiesBySource[this.source];
    for (var sourceLayerId in layerFamilies) {
        var sourceLayer = data.layers[sourceLayerId];
        if (!sourceLayer) {
            continue;
        }
        if (sourceLayer.version === 1) {
            util.warnOnce('Vector tile source "' + this$1.source + '" layer "' + sourceLayerId + '" does not use vector tile spec v2 ' + 'and therefore may have some rendering errors.');
        }
        var sourceLayerIndex = sourceLayerCoder.encode(sourceLayerId);
        var features = [];
        for (var i = 0; i < sourceLayer.length; i++) {
            var feature = sourceLayer.feature(i);
            feature.index = i;
            feature.sourceLayerIndex = sourceLayerIndex;
            features.push(feature);
        }
        for (var i$2 = 0, list = layerFamilies[sourceLayerId]; i$2 < list.length; i$2 += 1) {
            var family = list[i$2];
            var layer = family[0];
            if (layer.minzoom && this$1.zoom < layer.minzoom) {
                continue;
            }
            if (layer.maxzoom && this$1.zoom >= layer.maxzoom) {
                continue;
            }
            if (layer.layout && layer.layout.visibility === 'none') {
                continue;
            }
            for (var i$3 = 0, list$1 = family; i$3 < list$1.length; i$3 += 1) {
                var layer$1 = list$1[i$3];
                layer$1.recalculate(this$1.zoom);
            }
            var bucket = buckets[layer.id] = layer.createBucket({
                index: bucketIndex,
                layers: family,
                zoom: this$1.zoom,
                overscaling: this$1.overscaling,
                collisionBoxArray: this$1.collisionBoxArray
            });
            bucket.populate(features, options);
            featureIndex.bucketLayerIDs[bucketIndex] = family.map(function (l) {
                return l.id;
            });
            bucketIndex++;
        }
    }
    var done = function (collisionTile) {
        this$1.status = 'done';
        featureIndex.paintPropertyStatistics = {};
        for (var id in buckets) {
            util.extend(featureIndex.paintPropertyStatistics, buckets[id].getPaintPropertyStatistics());
        }
        var transferables = [];
        callback(null, {
            buckets: serializeBuckets(util.values(buckets), transferables),
            featureIndex: featureIndex.serialize(transferables),
            collisionTile: collisionTile.serialize(transferables),
            collisionBoxArray: this$1.collisionBoxArray.serialize()
        }, transferables);
    };
    this.symbolBuckets = [];
    for (var i$1 = layerIndex.symbolOrder.length - 1; i$1 >= 0; i$1--) {
        var bucket$1 = buckets[layerIndex.symbolOrder[i$1]];
        if (bucket$1) {
            this$1.symbolBuckets.push(bucket$1);
        }
    }
    if (this.symbolBuckets.length === 0) {
        return done(new CollisionTile(this.angle, this.pitch, this.collisionBoxArray));
    }
    var deps = 0;
    var icons = Object.keys(options.iconDependencies);
    var stacks = util.mapObject(options.glyphDependencies, function (glyphs) {
        return Object.keys(glyphs).map(Number);
    });
    var gotDependency = function (err) {
        if (err) {
            return callback(err);
        }
        deps++;
        if (deps === 2) {
            var collisionTile = new CollisionTile(this$1.angle, this$1.pitch, this$1.collisionBoxArray);
            for (var i = 0, list = this$1.symbolBuckets; i < list.length; i += 1) {
                var bucket = list[i];
                recalculateLayers(bucket, this$1.zoom);
                bucket.prepare(stacks, icons);
                bucket.place(collisionTile, this$1.showCollisionBoxes);
            }
            done(collisionTile);
        }
    };
    if (Object.keys(stacks).length) {
        actor.send('getGlyphs', {
            uid: this.uid,
            stacks: stacks
        }, function (err, newStacks) {
            stacks = newStacks;
            gotDependency(err);
        });
    } else {
        gotDependency();
    }
    if (icons.length) {
        actor.send('getIcons', { icons: icons }, function (err, newIcons) {
            icons = newIcons;
            gotDependency(err);
        });
    } else {
        gotDependency();
    }
};
WorkerTile.prototype.redoPlacement = function redoPlacement(angle, pitch, showCollisionBoxes) {
    var this$1 = this;
    this.angle = angle;
    this.pitch = pitch;
    if (this.status !== 'done') {
        return {};
    }
    var collisionTile = new CollisionTile(this.angle, this.pitch, this.collisionBoxArray);
    for (var i = 0, list = this$1.symbolBuckets; i < list.length; i += 1) {
        var bucket = list[i];
        recalculateLayers(bucket, this$1.zoom);
        bucket.place(collisionTile, showCollisionBoxes);
    }
    var transferables = [];
    return {
        result: {
            buckets: serializeBuckets(this.symbolBuckets, transferables),
            collisionTile: collisionTile.serialize(transferables)
        },
        transferables: transferables
    };
};
function recalculateLayers(bucket, zoom) {
    for (var i = 0, list = bucket.layers; i < list.length; i += 1) {
        var layer = list[i];
        layer.recalculate(zoom);
    }
}
function serializeBuckets(buckets, transferables) {
    return buckets.filter(function (b) {
        return !b.isEmpty();
    }).map(function (b) {
        return b.serialize(transferables);
    });
}
module.exports = WorkerTile;
},{"../data/feature_index":21,"../symbol/collision_box":30,"../symbol/collision_tile":31,"../util/dictionary_coder":34,"../util/util":38}],29:[function(require,module,exports){
'use strict';
module.exports = createFilter;
var types = [
    'Unknown',
    'Point',
    'LineString',
    'Polygon'
];
function createFilter(filter) {
    return new Function('f', 'var p = (f && f.properties || {}); return ' + compile(filter));
}
function compile(filter) {
    if (!filter) {
        return 'true';
    }
    var op = filter[0];
    if (filter.length <= 1) {
        return op === 'any' ? 'false' : 'true';
    }
    var str = op === '==' ? compileComparisonOp(filter[1], filter[2], '===', false) : op === '!=' ? compileComparisonOp(filter[1], filter[2], '!==', false) : op === '<' || op === '>' || op === '<=' || op === '>=' ? compileComparisonOp(filter[1], filter[2], op, true) : op === 'any' ? compileLogicalOp(filter.slice(1), '||') : op === 'all' ? compileLogicalOp(filter.slice(1), '&&') : op === 'none' ? compileNegation(compileLogicalOp(filter.slice(1), '||')) : op === 'in' ? compileInOp(filter[1], filter.slice(2)) : op === '!in' ? compileNegation(compileInOp(filter[1], filter.slice(2))) : op === 'has' ? compileHasOp(filter[1]) : op === '!has' ? compileNegation(compileHasOp(filter[1])) : 'true';
    return '(' + str + ')';
}
function compilePropertyReference(property) {
    return property === '$type' ? 'f.type' : property === '$id' ? 'f.id' : 'p[' + JSON.stringify(property) + ']';
}
function compileComparisonOp(property, value, op, checkType) {
    var left = compilePropertyReference(property);
    var right = property === '$type' ? types.indexOf(value) : JSON.stringify(value);
    return (checkType ? 'typeof ' + left + '=== typeof ' + right + '&&' : '') + left + op + right;
}
function compileLogicalOp(expressions, op) {
    return expressions.map(compile).join(op);
}
function compileInOp(property, values) {
    if (property === '$type') {
        values = values.map(function (value) {
            return types.indexOf(value);
        });
    }
    var left = JSON.stringify(values.sort(compare));
    var right = compilePropertyReference(property);
    if (values.length <= 200) {
        return left + '.indexOf(' + right + ') !== -1';
    }
    return '' + ('function(v, a, i, j) {' + 'while (i <= j) { var m = (i + j) >> 1;' + '    if (a[m] === v) return true; if (a[m] > v) j = m - 1; else i = m + 1;' + '}' + 'return false; }(') + right + ', ' + left + ',0,' + (values.length - 1) + ')';
}
function compileHasOp(property) {
    return property === '$id' ? '"id" in f' : JSON.stringify(property) + ' in p';
}
function compileNegation(expression) {
    return '!(' + expression + ')';
}
function compare(a, b) {
    return a < b ? -1 : a > b ? 1 : 0;
}
},{}],30:[function(require,module,exports){
'use strict';
var createStructArrayType = require('../util/struct_array');
var Point = require('point-geometry');
var CollisionBoxArray = createStructArrayType({
    members: [
        {
            type: 'Int16',
            name: 'anchorPointX'
        },
        {
            type: 'Int16',
            name: 'anchorPointY'
        },
        {
            type: 'Int16',
            name: 'x1'
        },
        {
            type: 'Int16',
            name: 'y1'
        },
        {
            type: 'Int16',
            name: 'x2'
        },
        {
            type: 'Int16',
            name: 'y2'
        },
        {
            type: 'Float32',
            name: 'maxScale'
        },
        {
            type: 'Uint32',
            name: 'featureIndex'
        },
        {
            type: 'Uint16',
            name: 'sourceLayerIndex'
        },
        {
            type: 'Uint16',
            name: 'bucketIndex'
        },
        {
            type: 'Int16',
            name: 'bbox0'
        },
        {
            type: 'Int16',
            name: 'bbox1'
        },
        {
            type: 'Int16',
            name: 'bbox2'
        },
        {
            type: 'Int16',
            name: 'bbox3'
        },
        {
            type: 'Float32',
            name: 'placementScale'
        }
    ]
});
Object.defineProperty(CollisionBoxArray.prototype.StructType.prototype, 'anchorPoint', {
    get: function get() {
        return new Point(this.anchorPointX, this.anchorPointY);
    }
});
module.exports = CollisionBoxArray;
},{"../util/struct_array":37,"point-geometry":42}],31:[function(require,module,exports){
'use strict';
var Point = require('point-geometry');
var EXTENT = require('../data/extent');
var Grid = require('grid-index');
var intersectionTests = require('../util/intersection_tests');
var CollisionTile = function CollisionTile(angle, pitch, collisionBoxArray) {
    if (typeof angle === 'object') {
        var serialized = angle;
        collisionBoxArray = pitch;
        angle = serialized.angle;
        pitch = serialized.pitch;
        this.grid = new Grid(serialized.grid);
        this.ignoredGrid = new Grid(serialized.ignoredGrid);
    } else {
        this.grid = new Grid(EXTENT, 12, 6);
        this.ignoredGrid = new Grid(EXTENT, 12, 0);
    }
    this.minScale = 0.5;
    this.maxScale = 2;
    this.angle = angle;
    this.pitch = pitch;
    var sin = Math.sin(angle), cos = Math.cos(angle);
    this.rotationMatrix = [
        cos,
        -sin,
        sin,
        cos
    ];
    this.reverseRotationMatrix = [
        cos,
        sin,
        -sin,
        cos
    ];
    this.yStretch = 1 / Math.cos(pitch / 180 * Math.PI);
    this.yStretch = Math.pow(this.yStretch, 1.3);
    this.collisionBoxArray = collisionBoxArray;
    if (collisionBoxArray.length === 0) {
        collisionBoxArray.emplaceBack();
        var maxInt16 = 32767;
        collisionBoxArray.emplaceBack(0, 0, 0, -maxInt16, 0, maxInt16, maxInt16, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        collisionBoxArray.emplaceBack(EXTENT, 0, 0, -maxInt16, 0, maxInt16, maxInt16, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        collisionBoxArray.emplaceBack(0, 0, -maxInt16, 0, maxInt16, 0, maxInt16, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        collisionBoxArray.emplaceBack(0, EXTENT, -maxInt16, 0, maxInt16, 0, maxInt16, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }
    this.tempCollisionBox = collisionBoxArray.get(0);
    this.edges = [
        collisionBoxArray.get(1),
        collisionBoxArray.get(2),
        collisionBoxArray.get(3),
        collisionBoxArray.get(4)
    ];
};
CollisionTile.prototype.serialize = function serialize(transferables) {
    var grid = this.grid.toArrayBuffer();
    var ignoredGrid = this.ignoredGrid.toArrayBuffer();
    if (transferables) {
        transferables.push(grid);
        transferables.push(ignoredGrid);
    }
    return {
        angle: this.angle,
        pitch: this.pitch,
        grid: grid,
        ignoredGrid: ignoredGrid
    };
};
CollisionTile.prototype.placeCollisionFeature = function placeCollisionFeature(collisionFeature, allowOverlap, avoidEdges) {
    var this$1 = this;
    var collisionBoxArray = this.collisionBoxArray;
    var minPlacementScale = this.minScale;
    var rotationMatrix = this.rotationMatrix;
    var yStretch = this.yStretch;
    for (var b = collisionFeature.boxStartIndex; b < collisionFeature.boxEndIndex; b++) {
        var box = collisionBoxArray.get(b);
        var anchorPoint = box.anchorPoint._matMult(rotationMatrix);
        var x = anchorPoint.x;
        var y = anchorPoint.y;
        var x1 = x + box.x1;
        var y1 = y + box.y1 * yStretch;
        var x2 = x + box.x2;
        var y2 = y + box.y2 * yStretch;
        box.bbox0 = x1;
        box.bbox1 = y1;
        box.bbox2 = x2;
        box.bbox3 = y2;
        if (!allowOverlap) {
            var blockingBoxes = this$1.grid.query(x1, y1, x2, y2);
            for (var i = 0; i < blockingBoxes.length; i++) {
                var blocking = collisionBoxArray.get(blockingBoxes[i]);
                var blockingAnchorPoint = blocking.anchorPoint._matMult(rotationMatrix);
                minPlacementScale = this$1.getPlacementScale(minPlacementScale, anchorPoint, box, blockingAnchorPoint, blocking);
                if (minPlacementScale >= this$1.maxScale) {
                    return minPlacementScale;
                }
            }
        }
        if (avoidEdges) {
            var rotatedCollisionBox = void 0;
            if (this$1.angle) {
                var reverseRotationMatrix = this$1.reverseRotationMatrix;
                var tl = new Point(box.x1, box.y1).matMult(reverseRotationMatrix);
                var tr = new Point(box.x2, box.y1).matMult(reverseRotationMatrix);
                var bl = new Point(box.x1, box.y2).matMult(reverseRotationMatrix);
                var br = new Point(box.x2, box.y2).matMult(reverseRotationMatrix);
                rotatedCollisionBox = this$1.tempCollisionBox;
                rotatedCollisionBox.anchorPointX = box.anchorPoint.x;
                rotatedCollisionBox.anchorPointY = box.anchorPoint.y;
                rotatedCollisionBox.x1 = Math.min(tl.x, tr.x, bl.x, br.x);
                rotatedCollisionBox.y1 = Math.min(tl.y, tr.x, bl.x, br.x);
                rotatedCollisionBox.x2 = Math.max(tl.x, tr.x, bl.x, br.x);
                rotatedCollisionBox.y2 = Math.max(tl.y, tr.x, bl.x, br.x);
                rotatedCollisionBox.maxScale = box.maxScale;
            } else {
                rotatedCollisionBox = box;
            }
            for (var k = 0; k < this.edges.length; k++) {
                var edgeBox = this$1.edges[k];
                minPlacementScale = this$1.getPlacementScale(minPlacementScale, box.anchorPoint, rotatedCollisionBox, edgeBox.anchorPoint, edgeBox);
                if (minPlacementScale >= this$1.maxScale) {
                    return minPlacementScale;
                }
            }
        }
    }
    return minPlacementScale;
};
CollisionTile.prototype.queryRenderedSymbols = function queryRenderedSymbols(queryGeometry, scale) {
    var sourceLayerFeatures = {};
    var result = [];
    if (queryGeometry.length === 0 || this.grid.length === 0 && this.ignoredGrid.length === 0) {
        return result;
    }
    var collisionBoxArray = this.collisionBoxArray;
    var rotationMatrix = this.rotationMatrix;
    var yStretch = this.yStretch;
    var rotatedQuery = [];
    var minX = Infinity;
    var minY = Infinity;
    var maxX = -Infinity;
    var maxY = -Infinity;
    for (var i = 0; i < queryGeometry.length; i++) {
        var ring = queryGeometry[i];
        for (var k = 0; k < ring.length; k++) {
            var p = ring[k].matMult(rotationMatrix);
            minX = Math.min(minX, p.x);
            minY = Math.min(minY, p.y);
            maxX = Math.max(maxX, p.x);
            maxY = Math.max(maxY, p.y);
            rotatedQuery.push(p);
        }
    }
    var features = this.grid.query(minX, minY, maxX, maxY);
    var ignoredFeatures = this.ignoredGrid.query(minX, minY, maxX, maxY);
    for (var i$1 = 0; i$1 < ignoredFeatures.length; i$1++) {
        features.push(ignoredFeatures[i$1]);
    }
    var roundedScale = Math.pow(2, Math.ceil(Math.log(scale) / Math.LN2 * 10) / 10);
    for (var i$2 = 0; i$2 < features.length; i$2++) {
        var blocking = collisionBoxArray.get(features[i$2]);
        var sourceLayer = blocking.sourceLayerIndex;
        var featureIndex = blocking.featureIndex;
        if (sourceLayerFeatures[sourceLayer] === undefined) {
            sourceLayerFeatures[sourceLayer] = {};
        }
        if (sourceLayerFeatures[sourceLayer][featureIndex]) {
            continue;
        }
        if (roundedScale < blocking.placementScale || roundedScale > blocking.maxScale) {
            continue;
        }
        var anchor = blocking.anchorPoint.matMult(rotationMatrix);
        var x1 = anchor.x + blocking.x1 / scale;
        var y1 = anchor.y + blocking.y1 / scale * yStretch;
        var x2 = anchor.x + blocking.x2 / scale;
        var y2 = anchor.y + blocking.y2 / scale * yStretch;
        var bbox = [
            new Point(x1, y1),
            new Point(x2, y1),
            new Point(x2, y2),
            new Point(x1, y2)
        ];
        if (!intersectionTests.polygonIntersectsPolygon(rotatedQuery, bbox)) {
            continue;
        }
        sourceLayerFeatures[sourceLayer][featureIndex] = true;
        result.push(features[i$2]);
    }
    return result;
};
CollisionTile.prototype.getPlacementScale = function getPlacementScale(minPlacementScale, anchorPoint, box, blockingAnchorPoint, blocking) {
    var anchorDiffX = anchorPoint.x - blockingAnchorPoint.x;
    var anchorDiffY = anchorPoint.y - blockingAnchorPoint.y;
    var s1 = (blocking.x1 - box.x2) / anchorDiffX;
    var s2 = (blocking.x2 - box.x1) / anchorDiffX;
    var s3 = (blocking.y1 - box.y2) * this.yStretch / anchorDiffY;
    var s4 = (blocking.y2 - box.y1) * this.yStretch / anchorDiffY;
    if (isNaN(s1) || isNaN(s2)) {
        s1 = s2 = 1;
    }
    if (isNaN(s3) || isNaN(s4)) {
        s3 = s4 = 1;
    }
    var collisionFreeScale = Math.min(Math.max(s1, s2), Math.max(s3, s4));
    var blockingMaxScale = blocking.maxScale;
    var boxMaxScale = box.maxScale;
    if (collisionFreeScale > blockingMaxScale) {
        collisionFreeScale = blockingMaxScale;
    }
    if (collisionFreeScale > boxMaxScale) {
        collisionFreeScale = boxMaxScale;
    }
    if (collisionFreeScale > minPlacementScale && collisionFreeScale >= blocking.placementScale) {
        minPlacementScale = collisionFreeScale;
    }
    return minPlacementScale;
};
CollisionTile.prototype.insertCollisionFeature = function insertCollisionFeature(collisionFeature, minPlacementScale, ignorePlacement) {
    var this$1 = this;
    var grid = ignorePlacement ? this.ignoredGrid : this.grid;
    var collisionBoxArray = this.collisionBoxArray;
    for (var k = collisionFeature.boxStartIndex; k < collisionFeature.boxEndIndex; k++) {
        var box = collisionBoxArray.get(k);
        box.placementScale = minPlacementScale;
        if (minPlacementScale < this$1.maxScale) {
            grid.insert(k, box.bbox0, box.bbox1, box.bbox2, box.bbox3);
        }
    }
};
module.exports = CollisionTile;
},{"../data/extent":20,"../util/intersection_tests":36,"grid-index":14,"point-geometry":42}],32:[function(require,module,exports){
'use strict';
var window = require('./window');
var AJAXError = function (Error) {
    function AJAXError(message, status) {
        Error.call(this, message);
        this.status = status;
    }
    if (Error)
        AJAXError.__proto__ = Error;
    AJAXError.prototype = Object.create(Error && Error.prototype);
    AJAXError.prototype.constructor = AJAXError;
    return AJAXError;
}(Error);
exports.getJSON = function (url, callback) {
    var xhr = new window.XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.setRequestHeader('Accept', 'application/json');
    xhr.onerror = function (e) {
        callback(e);
    };
    xhr.onload = function () {
        if (xhr.status >= 200 && xhr.status < 300 && xhr.response) {
            var data;
            try {
                data = JSON.parse(xhr.response);
            } catch (err) {
                return callback(err);
            }
            callback(null, data);
        } else {
            callback(new AJAXError(xhr.statusText, xhr.status));
        }
    };
    xhr.send();
    return xhr;
};
exports.getArrayBuffer = function (url, callback) {
    var xhr = new window.XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.responseType = 'arraybuffer';
    xhr.onerror = function (e) {
        callback(e);
    };
    xhr.onload = function () {
        if (xhr.response.byteLength === 0 && xhr.status === 200) {
            return callback(new Error('http status 200 returned without content.'));
        }
        if (xhr.status >= 200 && xhr.status < 300 && xhr.response) {
            callback(null, {
                data: xhr.response,
                cacheControl: xhr.getResponseHeader('Cache-Control'),
                expires: xhr.getResponseHeader('Expires')
            });
        } else {
            callback(new AJAXError(xhr.statusText, xhr.status));
        }
    };
    xhr.send();
    return xhr;
};
function sameOrigin(url) {
    var a = window.document.createElement('a');
    a.href = url;
    return a.protocol === window.document.location.protocol && a.host === window.document.location.host;
}
var transparentPngUrl = 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAC0lEQVQYV2NgAAIAAAUAAarVyFEAAAAASUVORK5CYII=';
exports.getImage = function (url, callback) {
    return exports.getArrayBuffer(url, function (err, imgData) {
        if (err) {
            return callback(err);
        }
        var img = new window.Image();
        var URL = window.URL || window.webkitURL;
        img.onload = function () {
            callback(null, img);
            URL.revokeObjectURL(img.src);
        };
        var blob = new window.Blob([new Uint8Array(imgData.data)], { type: 'image/png' });
        img.cacheControl = imgData.cacheControl;
        img.expires = imgData.expires;
        img.src = imgData.data.byteLength ? URL.createObjectURL(blob) : transparentPngUrl;
    });
};
exports.getVideo = function (urls, callback) {
    var video = window.document.createElement('video');
    video.onloadstart = function () {
        callback(null, video);
    };
    for (var i = 0; i < urls.length; i++) {
        var s = window.document.createElement('source');
        if (!sameOrigin(urls[i])) {
            video.crossOrigin = 'Anonymous';
        }
        s.src = urls[i];
        video.appendChild(s);
    }
    return video;
};
},{"./window":33}],33:[function(require,module,exports){
'use strict';
module.exports = self;
},{}],34:[function(require,module,exports){
'use strict';
var DictionaryCoder = function DictionaryCoder(strings) {
    var this$1 = this;
    this._stringToNumber = {};
    this._numberToString = [];
    for (var i = 0; i < strings.length; i++) {
        var string = strings[i];
        this$1._stringToNumber[string] = i;
        this$1._numberToString[i] = string;
    }
};
DictionaryCoder.prototype.encode = function encode(string) {
    return this._stringToNumber[string];
};
DictionaryCoder.prototype.decode = function decode(n) {
    return this._numberToString[n];
};
module.exports = DictionaryCoder;
},{}],35:[function(require,module,exports){
'use strict';
var util = require('./util');
function _addEventListener(type, listener, listenerList) {
    listenerList[type] = listenerList[type] || [];
    listenerList[type].push(listener);
}
function _removeEventListener(type, listener, listenerList) {
    if (listenerList && listenerList[type]) {
        var index = listenerList[type].indexOf(listener);
        if (index !== -1) {
            listenerList[type].splice(index, 1);
        }
    }
}
var Evented = function Evented() {
};
Evented.prototype.on = function on(type, listener) {
    this._listeners = this._listeners || {};
    _addEventListener(type, listener, this._listeners);
    return this;
};
Evented.prototype.off = function off(type, listener) {
    _removeEventListener(type, listener, this._listeners);
    _removeEventListener(type, listener, this._oneTimeListeners);
    return this;
};
Evented.prototype.once = function once(type, listener) {
    this._oneTimeListeners = this._oneTimeListeners || {};
    _addEventListener(type, listener, this._oneTimeListeners);
    return this;
};
Evented.prototype.fire = function fire(type, data) {
    var this$1 = this;
    if (this.listens(type)) {
        data = util.extend({}, data, {
            type: type,
            target: this
        });
        var listeners = this._listeners && this._listeners[type] ? this._listeners[type].slice() : [];
        for (var i = 0; i < listeners.length; i++) {
            listeners[i].call(this$1, data);
        }
        var oneTimeListeners = this._oneTimeListeners && this._oneTimeListeners[type] ? this._oneTimeListeners[type].slice() : [];
        for (var i$1 = 0; i$1 < oneTimeListeners.length; i$1++) {
            oneTimeListeners[i$1].call(this$1, data);
            _removeEventListener(type, oneTimeListeners[i$1], this$1._oneTimeListeners);
        }
        if (this._eventedParent) {
            this._eventedParent.fire(type, util.extend({}, data, typeof this._eventedParentData === 'function' ? this._eventedParentData() : this._eventedParentData));
        }
    } else if (util.endsWith(type, 'error')) {
        console.error(data && data.error || data || 'Empty error event');
    }
    return this;
};
Evented.prototype.listens = function listens(type) {
    return this._listeners && this._listeners[type] && this._listeners[type].length > 0 || this._oneTimeListeners && this._oneTimeListeners[type] && this._oneTimeListeners[type].length > 0 || this._eventedParent && this._eventedParent.listens(type);
};
Evented.prototype.setEventedParent = function setEventedParent(parent, data) {
    this._eventedParent = parent;
    this._eventedParentData = data;
    return this;
};
module.exports = Evented;
},{"./util":38}],36:[function(require,module,exports){
'use strict';
var isCounterClockwise = require('./util').isCounterClockwise;
module.exports = {
    multiPolygonIntersectsBufferedMultiPoint: multiPolygonIntersectsBufferedMultiPoint,
    multiPolygonIntersectsMultiPolygon: multiPolygonIntersectsMultiPolygon,
    multiPolygonIntersectsBufferedMultiLine: multiPolygonIntersectsBufferedMultiLine,
    polygonIntersectsPolygon: polygonIntersectsPolygon,
    distToSegmentSquared: distToSegmentSquared
};
function polygonIntersectsPolygon(polygonA, polygonB) {
    for (var i = 0; i < polygonA.length; i++) {
        if (polygonContainsPoint(polygonB, polygonA[i])) {
            return true;
        }
    }
    for (var i$1 = 0; i$1 < polygonB.length; i$1++) {
        if (polygonContainsPoint(polygonA, polygonB[i$1])) {
            return true;
        }
    }
    if (lineIntersectsLine(polygonA, polygonB)) {
        return true;
    }
    return false;
}
function multiPolygonIntersectsBufferedMultiPoint(multiPolygon, rings, radius) {
    for (var j = 0; j < multiPolygon.length; j++) {
        var polygon = multiPolygon[j];
        for (var i = 0; i < rings.length; i++) {
            var ring = rings[i];
            for (var k = 0; k < ring.length; k++) {
                var point = ring[k];
                if (polygonContainsPoint(polygon, point)) {
                    return true;
                }
                if (pointIntersectsBufferedLine(point, polygon, radius)) {
                    return true;
                }
            }
        }
    }
    return false;
}
function multiPolygonIntersectsMultiPolygon(multiPolygonA, multiPolygonB) {
    if (multiPolygonA.length === 1 && multiPolygonA[0].length === 1) {
        return multiPolygonContainsPoint(multiPolygonB, multiPolygonA[0][0]);
    }
    for (var m = 0; m < multiPolygonB.length; m++) {
        var ring = multiPolygonB[m];
        for (var n = 0; n < ring.length; n++) {
            if (multiPolygonContainsPoint(multiPolygonA, ring[n])) {
                return true;
            }
        }
    }
    for (var j = 0; j < multiPolygonA.length; j++) {
        var polygon = multiPolygonA[j];
        for (var i = 0; i < polygon.length; i++) {
            if (multiPolygonContainsPoint(multiPolygonB, polygon[i])) {
                return true;
            }
        }
        for (var k = 0; k < multiPolygonB.length; k++) {
            if (lineIntersectsLine(polygon, multiPolygonB[k])) {
                return true;
            }
        }
    }
    return false;
}
function multiPolygonIntersectsBufferedMultiLine(multiPolygon, multiLine, radius) {
    for (var i = 0; i < multiLine.length; i++) {
        var line = multiLine[i];
        for (var j = 0; j < multiPolygon.length; j++) {
            var polygon = multiPolygon[j];
            if (polygon.length >= 3) {
                for (var k = 0; k < line.length; k++) {
                    if (polygonContainsPoint(polygon, line[k])) {
                        return true;
                    }
                }
            }
            if (lineIntersectsBufferedLine(polygon, line, radius)) {
                return true;
            }
        }
    }
    return false;
}
function lineIntersectsBufferedLine(lineA, lineB, radius) {
    if (lineA.length > 1) {
        if (lineIntersectsLine(lineA, lineB)) {
            return true;
        }
        for (var j = 0; j < lineB.length; j++) {
            if (pointIntersectsBufferedLine(lineB[j], lineA, radius)) {
                return true;
            }
        }
    }
    for (var k = 0; k < lineA.length; k++) {
        if (pointIntersectsBufferedLine(lineA[k], lineB, radius)) {
            return true;
        }
    }
    return false;
}
function lineIntersectsLine(lineA, lineB) {
    if (lineA.length === 0 || lineB.length === 0) {
        return false;
    }
    for (var i = 0; i < lineA.length - 1; i++) {
        var a0 = lineA[i];
        var a1 = lineA[i + 1];
        for (var j = 0; j < lineB.length - 1; j++) {
            var b0 = lineB[j];
            var b1 = lineB[j + 1];
            if (lineSegmentIntersectsLineSegment(a0, a1, b0, b1)) {
                return true;
            }
        }
    }
    return false;
}
function lineSegmentIntersectsLineSegment(a0, a1, b0, b1) {
    return isCounterClockwise(a0, b0, b1) !== isCounterClockwise(a1, b0, b1) && isCounterClockwise(a0, a1, b0) !== isCounterClockwise(a0, a1, b1);
}
function pointIntersectsBufferedLine(p, line, radius) {
    var radiusSquared = radius * radius;
    if (line.length === 1) {
        return p.distSqr(line[0]) < radiusSquared;
    }
    for (var i = 1; i < line.length; i++) {
        var v = line[i - 1], w = line[i];
        if (distToSegmentSquared(p, v, w) < radiusSquared) {
            return true;
        }
    }
    return false;
}
function distToSegmentSquared(p, v, w) {
    var l2 = v.distSqr(w);
    if (l2 === 0) {
        return p.distSqr(v);
    }
    var t = ((p.x - v.x) * (w.x - v.x) + (p.y - v.y) * (w.y - v.y)) / l2;
    if (t < 0) {
        return p.distSqr(v);
    }
    if (t > 1) {
        return p.distSqr(w);
    }
    return p.distSqr(w.sub(v)._mult(t)._add(v));
}
function multiPolygonContainsPoint(rings, p) {
    var c = false, ring, p1, p2;
    for (var k = 0; k < rings.length; k++) {
        ring = rings[k];
        for (var i = 0, j = ring.length - 1; i < ring.length; j = i++) {
            p1 = ring[i];
            p2 = ring[j];
            if (p1.y > p.y !== p2.y > p.y && p.x < (p2.x - p1.x) * (p.y - p1.y) / (p2.y - p1.y) + p1.x) {
                c = !c;
            }
        }
    }
    return c;
}
function polygonContainsPoint(ring, p) {
    var c = false;
    for (var i = 0, j = ring.length - 1; i < ring.length; j = i++) {
        var p1 = ring[i];
        var p2 = ring[j];
        if (p1.y > p.y !== p2.y > p.y && p.x < (p2.x - p1.x) * (p.y - p1.y) / (p2.y - p1.y) + p1.x) {
            c = !c;
        }
    }
    return c;
}
},{"./util":38}],37:[function(require,module,exports){
'use strict';
module.exports = createStructArrayType;
var viewTypes = {
    'Int8': Int8Array,
    'Uint8': Uint8Array,
    'Uint8Clamped': Uint8ClampedArray,
    'Int16': Int16Array,
    'Uint16': Uint16Array,
    'Int32': Int32Array,
    'Uint32': Uint32Array,
    'Float32': Float32Array,
    'Float64': Float64Array
};
var Struct = function Struct(structArray, index) {
    this._structArray = structArray;
    this._pos1 = index * this.size;
    this._pos2 = this._pos1 / 2;
    this._pos4 = this._pos1 / 4;
    this._pos8 = this._pos1 / 8;
};
var DEFAULT_CAPACITY = 128;
var RESIZE_MULTIPLIER = 5;
var StructArray = function StructArray(serialized) {
    this.isTransferred = false;
    if (serialized !== undefined) {
        this.arrayBuffer = serialized.arrayBuffer;
        this.length = serialized.length;
        this.capacity = this.arrayBuffer.byteLength / this.bytesPerElement;
        this._refreshViews();
    } else {
        this.capacity = -1;
        this.resize(0);
    }
};
StructArray.serialize = function serialize() {
    return {
        members: this.prototype.members,
        alignment: this.prototype.StructType.prototype.alignment,
        bytesPerElement: this.prototype.bytesPerElement
    };
};
StructArray.prototype.serialize = function serialize(transferables) {
    this._trim();
    if (transferables) {
        this.isTransferred = true;
        transferables.push(this.arrayBuffer);
    }
    return {
        length: this.length,
        arrayBuffer: this.arrayBuffer
    };
};
StructArray.prototype.get = function get(index) {
    return new this.StructType(this, index);
};
StructArray.prototype._trim = function _trim() {
    if (this.length !== this.capacity) {
        this.capacity = this.length;
        this.arrayBuffer = this.arrayBuffer.slice(0, this.length * this.bytesPerElement);
        this._refreshViews();
    }
};
StructArray.prototype.resize = function resize(n) {
    this.length = n;
    if (n > this.capacity) {
        this.capacity = Math.max(n, Math.floor(this.capacity * RESIZE_MULTIPLIER), DEFAULT_CAPACITY);
        this.arrayBuffer = new ArrayBuffer(this.capacity * this.bytesPerElement);
        var oldUint8Array = this.uint8;
        this._refreshViews();
        if (oldUint8Array) {
            this.uint8.set(oldUint8Array);
        }
    }
};
StructArray.prototype._refreshViews = function _refreshViews() {
    var this$1 = this;
    for (var i = 0, list = this$1._usedTypes; i < list.length; i += 1) {
        var type = list[i];
        this$1[getArrayViewName(type)] = new viewTypes[type](this$1.arrayBuffer);
    }
};
StructArray.prototype.toArray = function toArray(startIndex, endIndex) {
    var this$1 = this;
    var array = [];
    for (var i = startIndex; i < endIndex; i++) {
        var struct = this$1.get(i);
        array.push(struct);
    }
    return array;
};
var structArrayTypeCache = {};
function createStructArrayType(options) {
    var key = JSON.stringify(options);
    if (structArrayTypeCache[key]) {
        return structArrayTypeCache[key];
    }
    var alignment = options.alignment === undefined ? 1 : options.alignment;
    var offset = 0;
    var maxSize = 0;
    var usedTypes = ['Uint8'];
    var members = options.members.map(function (member) {
        if (usedTypes.indexOf(member.type) < 0) {
            usedTypes.push(member.type);
        }
        var typeSize = sizeOf(member.type);
        var memberOffset = offset = align(offset, Math.max(alignment, typeSize));
        var components = member.components || 1;
        maxSize = Math.max(maxSize, typeSize);
        offset += typeSize * components;
        return {
            name: member.name,
            type: member.type,
            components: components,
            offset: memberOffset
        };
    });
    var size = align(offset, Math.max(maxSize, alignment));
    var StructType = function (Struct) {
        function StructType() {
            Struct.apply(this, arguments);
        }
        if (Struct)
            StructType.__proto__ = Struct;
        StructType.prototype = Object.create(Struct && Struct.prototype);
        StructType.prototype.constructor = StructType;
        return StructType;
    }(Struct);
    StructType.prototype.alignment = alignment;
    StructType.prototype.size = size;
    for (var i = 0, list = members; i < list.length; i += 1) {
        var member = list[i];
        for (var c = 0; c < member.components; c++) {
            var name = member.name + (member.components === 1 ? '' : c);
            Object.defineProperty(StructType.prototype, name, {
                get: createGetter(member, c),
                set: createSetter(member, c)
            });
        }
    }
    var StructArrayType = function (StructArray) {
        function StructArrayType() {
            StructArray.apply(this, arguments);
        }
        if (StructArray)
            StructArrayType.__proto__ = StructArray;
        StructArrayType.prototype = Object.create(StructArray && StructArray.prototype);
        StructArrayType.prototype.constructor = StructArrayType;
        return StructArrayType;
    }(StructArray);
    StructArrayType.prototype.members = members;
    StructArrayType.prototype.StructType = StructType;
    StructArrayType.prototype.bytesPerElement = size;
    StructArrayType.prototype.emplaceBack = createEmplaceBack(members, size);
    StructArrayType.prototype._usedTypes = usedTypes;
    structArrayTypeCache[key] = StructArrayType;
    return StructArrayType;
}
function align(offset, size) {
    return Math.ceil(offset / size) * size;
}
function sizeOf(type) {
    return viewTypes[type].BYTES_PER_ELEMENT;
}
function getArrayViewName(type) {
    return type.toLowerCase();
}
function createEmplaceBack(members, bytesPerElement) {
    var usedTypeSizes = [];
    var argNames = [];
    var body = 'var i = this.length;\n' + 'this.resize(this.length + 1);\n';
    for (var i = 0, list = members; i < list.length; i += 1) {
        var member = list[i];
        var size = sizeOf(member.type);
        if (usedTypeSizes.indexOf(size) < 0) {
            usedTypeSizes.push(size);
            body += 'var o' + size.toFixed(0) + ' = i * ' + (bytesPerElement / size).toFixed(0) + ';\n';
        }
        for (var c = 0; c < member.components; c++) {
            var argName = 'v' + argNames.length;
            var index = 'o' + size.toFixed(0) + ' + ' + (member.offset / size + c).toFixed(0);
            body += 'this.' + getArrayViewName(member.type) + '[' + index + '] = ' + argName + ';\n';
            argNames.push(argName);
        }
    }
    body += 'return i;';
    return new Function(argNames.toString(), body);
}
function createMemberComponentString(member, component) {
    var elementOffset = 'this._pos' + sizeOf(member.type).toFixed(0);
    var componentOffset = (member.offset / sizeOf(member.type) + component).toFixed(0);
    var index = elementOffset + ' + ' + componentOffset;
    return 'this._structArray.' + getArrayViewName(member.type) + '[' + index + ']';
}
function createGetter(member, c) {
    return new Function('return ' + createMemberComponentString(member, c) + ';');
}
function createSetter(member, c) {
    return new Function('x', createMemberComponentString(member, c) + ' = x;');
}
},{}],38:[function(require,module,exports){
'use strict';
var UnitBezier = require('@mapbox/unitbezier');
var Coordinate = require('../geo/coordinate');
var Point = require('point-geometry');
exports.easeCubicInOut = function (t) {
    if (t <= 0) {
        return 0;
    }
    if (t >= 1) {
        return 1;
    }
    var t2 = t * t, t3 = t2 * t;
    return 4 * (t < 0.5 ? t3 : 3 * (t - t2) + t3 - 0.75);
};
exports.bezier = function (p1x, p1y, p2x, p2y) {
    var bezier = new UnitBezier(p1x, p1y, p2x, p2y);
    return function (t) {
        return bezier.solve(t);
    };
};
exports.ease = exports.bezier(0.25, 0.1, 0.25, 1);
exports.clamp = function (n, min, max) {
    return Math.min(max, Math.max(min, n));
};
exports.wrap = function (n, min, max) {
    var d = max - min;
    var w = ((n - min) % d + d) % d + min;
    return w === min ? max : w;
};
exports.asyncAll = function (array, fn, callback) {
    if (!array.length) {
        return callback(null, []);
    }
    var remaining = array.length;
    var results = new Array(array.length);
    var error = null;
    array.forEach(function (item, i) {
        fn(item, function (err, result) {
            if (err) {
                error = err;
            }
            results[i] = result;
            if (--remaining === 0) {
                callback(error, results);
            }
        });
    });
};
exports.values = function (obj) {
    var result = [];
    for (var k in obj) {
        result.push(obj[k]);
    }
    return result;
};
exports.keysDifference = function (obj, other) {
    var difference = [];
    for (var i in obj) {
        if (!(i in other)) {
            difference.push(i);
        }
    }
    return difference;
};
exports.extend = function (dest, source0, source1, source2) {
    var arguments$1 = arguments;
    for (var i = 1; i < arguments.length; i++) {
        var src = arguments$1[i];
        for (var k in src) {
            dest[k] = src[k];
        }
    }
    return dest;
};
exports.pick = function (src, properties) {
    var result = {};
    for (var i = 0; i < properties.length; i++) {
        var k = properties[i];
        if (k in src) {
            result[k] = src[k];
        }
    }
    return result;
};
var id = 1;
exports.uniqueId = function () {
    return id++;
};
exports.bindAll = function (fns, context) {
    fns.forEach(function (fn) {
        if (!context[fn]) {
            return;
        }
        context[fn] = context[fn].bind(context);
    });
};
exports.getCoordinatesCenter = function (coords) {
    var minX = Infinity;
    var minY = Infinity;
    var maxX = -Infinity;
    var maxY = -Infinity;
    for (var i = 0; i < coords.length; i++) {
        minX = Math.min(minX, coords[i].column);
        minY = Math.min(minY, coords[i].row);
        maxX = Math.max(maxX, coords[i].column);
        maxY = Math.max(maxY, coords[i].row);
    }
    var dx = maxX - minX;
    var dy = maxY - minY;
    var dMax = Math.max(dx, dy);
    var zoom = Math.max(0, Math.floor(-Math.log(dMax) / Math.LN2));
    return new Coordinate((minX + maxX) / 2, (minY + maxY) / 2, 0).zoomTo(zoom);
};
exports.endsWith = function (string, suffix) {
    return string.indexOf(suffix, string.length - suffix.length) !== -1;
};
exports.mapObject = function (input, iterator, context) {
    var this$1 = this;
    var output = {};
    for (var key in input) {
        output[key] = iterator.call(context || this$1, input[key], key, input);
    }
    return output;
};
exports.filterObject = function (input, iterator, context) {
    var this$1 = this;
    var output = {};
    for (var key in input) {
        if (iterator.call(context || this$1, input[key], key, input)) {
            output[key] = input[key];
        }
    }
    return output;
};
exports.deepEqual = function (a, b) {
    if (Array.isArray(a)) {
        if (!Array.isArray(b) || a.length !== b.length) {
            return false;
        }
        for (var i = 0; i < a.length; i++) {
            if (!exports.deepEqual(a[i], b[i])) {
                return false;
            }
        }
        return true;
    }
    if (typeof a === 'object' && a !== null && b !== null) {
        if (!(typeof b === 'object')) {
            return false;
        }
        var keys = Object.keys(a);
        if (keys.length !== Object.keys(b).length) {
            return false;
        }
        for (var key in a) {
            if (!exports.deepEqual(a[key], b[key])) {
                return false;
            }
        }
        return true;
    }
    return a === b;
};
exports.clone = function (input) {
    if (Array.isArray(input)) {
        return input.map(exports.clone);
    } else if (typeof input === 'object' && input) {
        return exports.mapObject(input, exports.clone);
    } else {
        return input;
    }
};
exports.arraysIntersect = function (a, b) {
    for (var l = 0; l < a.length; l++) {
        if (b.indexOf(a[l]) >= 0) {
            return true;
        }
    }
    return false;
};
var warnOnceHistory = {};
exports.warnOnce = function (message) {
    if (!warnOnceHistory[message]) {
        if (typeof console !== 'undefined') {
            console.warn(message);
        }
        warnOnceHistory[message] = true;
    }
};
exports.isCounterClockwise = function (a, b, c) {
    return (c.y - a.y) * (b.x - a.x) > (b.y - a.y) * (c.x - a.x);
};
exports.calculateSignedArea = function (ring) {
    var sum = 0;
    for (var i = 0, len = ring.length, j = len - 1, p1 = void 0, p2 = void 0; i < len; j = i++) {
        p1 = ring[i];
        p2 = ring[j];
        sum += (p2.x - p1.x) * (p1.y + p2.y);
    }
    return sum;
};
exports.isClosedPolygon = function (points) {
    if (points.length < 4) {
        return false;
    }
    var p1 = points[0];
    var p2 = points[points.length - 1];
    if (Math.abs(p1.x - p2.x) > 0 || Math.abs(p1.y - p2.y) > 0) {
        return false;
    }
    return Math.abs(exports.calculateSignedArea(points)) > 0.01;
};
exports.sphericalToCartesian = function (spherical) {
    var r = spherical[0];
    var azimuthal = spherical[1], polar = spherical[2];
    azimuthal += 90;
    azimuthal *= Math.PI / 180;
    polar *= Math.PI / 180;
    return [
        r * Math.cos(azimuthal) * Math.sin(polar),
        r * Math.sin(azimuthal) * Math.sin(polar),
        r * Math.cos(polar)
    ];
};
exports.parseCacheControl = function (cacheControl) {
    var re = /(?:^|(?:\s*\,\s*))([^\x00-\x20\(\)<>@\,;\:\\"\/\[\]\?\=\{\}\x7F]+)(?:\=(?:([^\x00-\x20\(\)<>@\,;\:\\"\/\[\]\?\=\{\}\x7F]+)|(?:\"((?:[^"\\]|\\.)*)\")))?/g;
    var header = {};
    cacheControl.replace(re, function ($0, $1, $2, $3) {
        var value = $2 || $3;
        header[$1] = value ? value.toLowerCase() : true;
        return '';
    });
    if (header['max-age']) {
        var maxAge = parseInt(header['max-age'], 10);
        if (isNaN(maxAge)) {
            delete header['max-age'];
        } else {
            header['max-age'] = maxAge;
        }
    }
    return header;
};
},{"../geo/coordinate":23,"@mapbox/unitbezier":3,"point-geometry":42}],39:[function(require,module,exports){
'use strict';
var Feature = function Feature(vectorTileFeature, z, x, y) {
    this.type = 'Feature';
    this._vectorTileFeature = vectorTileFeature;
    vectorTileFeature._z = z;
    vectorTileFeature._x = x;
    vectorTileFeature._y = y;
    this.properties = vectorTileFeature.properties;
    if (vectorTileFeature.id != null) {
        this.id = vectorTileFeature.id;
    }
};
var prototypeAccessors = { geometry: {} };
prototypeAccessors.geometry.get = function () {
    if (this._geometry === undefined) {
        this._geometry = this._vectorTileFeature.toGeoJSON(this._vectorTileFeature._x, this._vectorTileFeature._y, this._vectorTileFeature._z).geometry;
    }
    return this._geometry;
};
prototypeAccessors.geometry.set = function (g) {
    this._geometry = g;
};
Feature.prototype.toJSON = function toJSON() {
    var this$1 = this;
    var json = { geometry: this.geometry };
    for (var i in this$1) {
        if (i === '_geometry' || i === '_vectorTileFeature') {
            continue;
        }
        json[i] = this$1[i];
    }
    return json;
};
Object.defineProperties(Feature.prototype, prototypeAccessors);
module.exports = Feature;
},{}],40:[function(require,module,exports){
'use strict';

// lightweight Buffer shim for pbf browser build
// based on code from github.com/feross/buffer (MIT-licensed)

module.exports = Buffer;

var ieee754 = require('ieee754');

var BufferMethods;

function Buffer(length) {
    var arr;
    if (length && length.length) {
        arr = length;
        length = arr.length;
    }
    var buf = new Uint8Array(length || 0);
    if (arr) buf.set(arr);

    buf.readUInt32LE = BufferMethods.readUInt32LE;
    buf.writeUInt32LE = BufferMethods.writeUInt32LE;
    buf.readInt32LE = BufferMethods.readInt32LE;
    buf.writeInt32LE = BufferMethods.writeInt32LE;
    buf.readFloatLE = BufferMethods.readFloatLE;
    buf.writeFloatLE = BufferMethods.writeFloatLE;
    buf.readDoubleLE = BufferMethods.readDoubleLE;
    buf.writeDoubleLE = BufferMethods.writeDoubleLE;
    buf.toString = BufferMethods.toString;
    buf.write = BufferMethods.write;
    buf.slice = BufferMethods.slice;
    buf.copy = BufferMethods.copy;

    buf._isBuffer = true;
    return buf;
}

var lastStr, lastStrEncoded;

BufferMethods = {
    readUInt32LE: function(pos) {
        return ((this[pos]) |
            (this[pos + 1] << 8) |
            (this[pos + 2] << 16)) +
            (this[pos + 3] * 0x1000000);
    },

    writeUInt32LE: function(val, pos) {
        this[pos] = val;
        this[pos + 1] = (val >>> 8);
        this[pos + 2] = (val >>> 16);
        this[pos + 3] = (val >>> 24);
    },

    readInt32LE: function(pos) {
        return ((this[pos]) |
            (this[pos + 1] << 8) |
            (this[pos + 2] << 16)) +
            (this[pos + 3] << 24);
    },

    readFloatLE:  function(pos) { return ieee754.read(this, pos, true, 23, 4); },
    readDoubleLE: function(pos) { return ieee754.read(this, pos, true, 52, 8); },

    writeFloatLE:  function(val, pos) { return ieee754.write(this, val, pos, true, 23, 4); },
    writeDoubleLE: function(val, pos) { return ieee754.write(this, val, pos, true, 52, 8); },

    toString: function(encoding, start, end) {
        var str = '',
            tmp = '';

        start = start || 0;
        end = Math.min(this.length, end || this.length);

        for (var i = start; i < end; i++) {
            var ch = this[i];
            if (ch <= 0x7F) {
                str += decodeURIComponent(tmp) + String.fromCharCode(ch);
                tmp = '';
            } else {
                tmp += '%' + ch.toString(16);
            }
        }

        str += decodeURIComponent(tmp);

        return str;
    },

    write: function(str, pos) {
        var bytes = str === lastStr ? lastStrEncoded : encodeString(str);
        for (var i = 0; i < bytes.length; i++) {
            this[pos + i] = bytes[i];
        }
    },

    slice: function(start, end) {
        return this.subarray(start, end);
    },

    copy: function(buf, pos) {
        pos = pos || 0;
        for (var i = 0; i < this.length; i++) {
            buf[pos + i] = this[i];
        }
    }
};

BufferMethods.writeInt32LE = BufferMethods.writeUInt32LE;

Buffer.byteLength = function(str) {
    lastStr = str;
    lastStrEncoded = encodeString(str);
    return lastStrEncoded.length;
};

Buffer.isBuffer = function(buf) {
    return !!(buf && buf._isBuffer);
};

function encodeString(str) {
    var length = str.length,
        bytes = [];

    for (var i = 0, c, lead; i < length; i++) {
        c = str.charCodeAt(i); // code point

        if (c > 0xD7FF && c < 0xE000) {

            if (lead) {
                if (c < 0xDC00) {
                    bytes.push(0xEF, 0xBF, 0xBD);
                    lead = c;
                    continue;

                } else {
                    c = lead - 0xD800 << 10 | c - 0xDC00 | 0x10000;
                    lead = null;
                }

            } else {
                if (c > 0xDBFF || (i + 1 === length)) bytes.push(0xEF, 0xBF, 0xBD);
                else lead = c;

                continue;
            }

        } else if (lead) {
            bytes.push(0xEF, 0xBF, 0xBD);
            lead = null;
        }

        if (c < 0x80) bytes.push(c);
        else if (c < 0x800) bytes.push(c >> 0x6 | 0xC0, c & 0x3F | 0x80);
        else if (c < 0x10000) bytes.push(c >> 0xC | 0xE0, c >> 0x6 & 0x3F | 0x80, c & 0x3F | 0x80);
        else bytes.push(c >> 0x12 | 0xF0, c >> 0xC & 0x3F | 0x80, c >> 0x6 & 0x3F | 0x80, c & 0x3F | 0x80);
    }
    return bytes;
}

},{"ieee754":15}],41:[function(require,module,exports){
(function (global){
'use strict';

module.exports = Pbf;

var Buffer = global.Buffer || require('./buffer');

function Pbf(buf) {
    this.buf = !Buffer.isBuffer(buf) ? new Buffer(buf || 0) : buf;
    this.pos = 0;
    this.length = this.buf.length;
}

Pbf.Varint  = 0; // varint: int32, int64, uint32, uint64, sint32, sint64, bool, enum
Pbf.Fixed64 = 1; // 64-bit: double, fixed64, sfixed64
Pbf.Bytes   = 2; // length-delimited: string, bytes, embedded messages, packed repeated fields
Pbf.Fixed32 = 5; // 32-bit: float, fixed32, sfixed32

var SHIFT_LEFT_32 = (1 << 16) * (1 << 16),
    SHIFT_RIGHT_32 = 1 / SHIFT_LEFT_32,
    POW_2_63 = Math.pow(2, 63);

Pbf.prototype = {

    destroy: function() {
        this.buf = null;
    },

    // === READING =================================================================

    readFields: function(readField, result, end) {
        end = end || this.length;

        while (this.pos < end) {
            var val = this.readVarint(),
                tag = val >> 3,
                startPos = this.pos;

            readField(tag, result, this);

            if (this.pos === startPos) this.skip(val);
        }
        return result;
    },

    readMessage: function(readField, result) {
        return this.readFields(readField, result, this.readVarint() + this.pos);
    },

    readFixed32: function() {
        var val = this.buf.readUInt32LE(this.pos);
        this.pos += 4;
        return val;
    },

    readSFixed32: function() {
        var val = this.buf.readInt32LE(this.pos);
        this.pos += 4;
        return val;
    },

    // 64-bit int handling is based on github.com/dpw/node-buffer-more-ints (MIT-licensed)

    readFixed64: function() {
        var val = this.buf.readUInt32LE(this.pos) + this.buf.readUInt32LE(this.pos + 4) * SHIFT_LEFT_32;
        this.pos += 8;
        return val;
    },

    readSFixed64: function() {
        var val = this.buf.readUInt32LE(this.pos) + this.buf.readInt32LE(this.pos + 4) * SHIFT_LEFT_32;
        this.pos += 8;
        return val;
    },

    readFloat: function() {
        var val = this.buf.readFloatLE(this.pos);
        this.pos += 4;
        return val;
    },

    readDouble: function() {
        var val = this.buf.readDoubleLE(this.pos);
        this.pos += 8;
        return val;
    },

    readVarint: function() {
        var buf = this.buf,
            val, b;

        b = buf[this.pos++]; val  =  b & 0x7f;        if (b < 0x80) return val;
        b = buf[this.pos++]; val |= (b & 0x7f) << 7;  if (b < 0x80) return val;
        b = buf[this.pos++]; val |= (b & 0x7f) << 14; if (b < 0x80) return val;
        b = buf[this.pos++]; val |= (b & 0x7f) << 21; if (b < 0x80) return val;

        return readVarintRemainder(val, this);
    },

    readVarint64: function() {
        var startPos = this.pos,
            val = this.readVarint();

        if (val < POW_2_63) return val;

        var pos = this.pos - 2;
        while (this.buf[pos] === 0xff) pos--;
        if (pos < startPos) pos = startPos;

        val = 0;
        for (var i = 0; i < pos - startPos + 1; i++) {
            var b = ~this.buf[startPos + i] & 0x7f;
            val += i < 4 ? b << i * 7 : b * Math.pow(2, i * 7);
        }

        return -val - 1;
    },

    readSVarint: function() {
        var num = this.readVarint();
        return num % 2 === 1 ? (num + 1) / -2 : num / 2; // zigzag encoding
    },

    readBoolean: function() {
        return Boolean(this.readVarint());
    },

    readString: function() {
        var end = this.readVarint() + this.pos,
            str = this.buf.toString('utf8', this.pos, end);
        this.pos = end;
        return str;
    },

    readBytes: function() {
        var end = this.readVarint() + this.pos,
            buffer = this.buf.slice(this.pos, end);
        this.pos = end;
        return buffer;
    },

    // verbose for performance reasons; doesn't affect gzipped size

    readPackedVarint: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readVarint());
        return arr;
    },
    readPackedSVarint: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readSVarint());
        return arr;
    },
    readPackedBoolean: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readBoolean());
        return arr;
    },
    readPackedFloat: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readFloat());
        return arr;
    },
    readPackedDouble: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readDouble());
        return arr;
    },
    readPackedFixed32: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readFixed32());
        return arr;
    },
    readPackedSFixed32: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readSFixed32());
        return arr;
    },
    readPackedFixed64: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readFixed64());
        return arr;
    },
    readPackedSFixed64: function() {
        var end = this.readVarint() + this.pos, arr = [];
        while (this.pos < end) arr.push(this.readSFixed64());
        return arr;
    },

    skip: function(val) {
        var type = val & 0x7;
        if (type === Pbf.Varint) while (this.buf[this.pos++] > 0x7f) {}
        else if (type === Pbf.Bytes) this.pos = this.readVarint() + this.pos;
        else if (type === Pbf.Fixed32) this.pos += 4;
        else if (type === Pbf.Fixed64) this.pos += 8;
        else throw new Error('Unimplemented type: ' + type);
    },

    // === WRITING =================================================================

    writeTag: function(tag, type) {
        this.writeVarint((tag << 3) | type);
    },

    realloc: function(min) {
        var length = this.length || 16;

        while (length < this.pos + min) length *= 2;

        if (length !== this.length) {
            var buf = new Buffer(length);
            this.buf.copy(buf);
            this.buf = buf;
            this.length = length;
        }
    },

    finish: function() {
        this.length = this.pos;
        this.pos = 0;
        return this.buf.slice(0, this.length);
    },

    writeFixed32: function(val) {
        this.realloc(4);
        this.buf.writeUInt32LE(val, this.pos);
        this.pos += 4;
    },

    writeSFixed32: function(val) {
        this.realloc(4);
        this.buf.writeInt32LE(val, this.pos);
        this.pos += 4;
    },

    writeFixed64: function(val) {
        this.realloc(8);
        this.buf.writeInt32LE(val & -1, this.pos);
        this.buf.writeUInt32LE(Math.floor(val * SHIFT_RIGHT_32), this.pos + 4);
        this.pos += 8;
    },

    writeSFixed64: function(val) {
        this.realloc(8);
        this.buf.writeInt32LE(val & -1, this.pos);
        this.buf.writeInt32LE(Math.floor(val * SHIFT_RIGHT_32), this.pos + 4);
        this.pos += 8;
    },

    writeVarint: function(val) {
        val = +val;

        if (val > 0xfffffff) {
            writeBigVarint(val, this);
            return;
        }

        this.realloc(4);

        this.buf[this.pos++] =           val & 0x7f  | (val > 0x7f ? 0x80 : 0); if (val <= 0x7f) return;
        this.buf[this.pos++] = ((val >>>= 7) & 0x7f) | (val > 0x7f ? 0x80 : 0); if (val <= 0x7f) return;
        this.buf[this.pos++] = ((val >>>= 7) & 0x7f) | (val > 0x7f ? 0x80 : 0); if (val <= 0x7f) return;
        this.buf[this.pos++] =   (val >>> 7) & 0x7f;
    },

    writeSVarint: function(val) {
        this.writeVarint(val < 0 ? -val * 2 - 1 : val * 2);
    },

    writeBoolean: function(val) {
        this.writeVarint(Boolean(val));
    },

    writeString: function(str) {
        str = String(str);
        var bytes = Buffer.byteLength(str);
        this.writeVarint(bytes);
        this.realloc(bytes);
        this.buf.write(str, this.pos);
        this.pos += bytes;
    },

    writeFloat: function(val) {
        this.realloc(4);
        this.buf.writeFloatLE(val, this.pos);
        this.pos += 4;
    },

    writeDouble: function(val) {
        this.realloc(8);
        this.buf.writeDoubleLE(val, this.pos);
        this.pos += 8;
    },

    writeBytes: function(buffer) {
        var len = buffer.length;
        this.writeVarint(len);
        this.realloc(len);
        for (var i = 0; i < len; i++) this.buf[this.pos++] = buffer[i];
    },

    writeRawMessage: function(fn, obj) {
        this.pos++; // reserve 1 byte for short message length

        // write the message directly to the buffer and see how much was written
        var startPos = this.pos;
        fn(obj, this);
        var len = this.pos - startPos;

        if (len >= 0x80) reallocForRawMessage(startPos, len, this);

        // finally, write the message length in the reserved place and restore the position
        this.pos = startPos - 1;
        this.writeVarint(len);
        this.pos += len;
    },

    writeMessage: function(tag, fn, obj) {
        this.writeTag(tag, Pbf.Bytes);
        this.writeRawMessage(fn, obj);
    },

    writePackedVarint:   function(tag, arr) { this.writeMessage(tag, writePackedVarint, arr);   },
    writePackedSVarint:  function(tag, arr) { this.writeMessage(tag, writePackedSVarint, arr);  },
    writePackedBoolean:  function(tag, arr) { this.writeMessage(tag, writePackedBoolean, arr);  },
    writePackedFloat:    function(tag, arr) { this.writeMessage(tag, writePackedFloat, arr);    },
    writePackedDouble:   function(tag, arr) { this.writeMessage(tag, writePackedDouble, arr);   },
    writePackedFixed32:  function(tag, arr) { this.writeMessage(tag, writePackedFixed32, arr);  },
    writePackedSFixed32: function(tag, arr) { this.writeMessage(tag, writePackedSFixed32, arr); },
    writePackedFixed64:  function(tag, arr) { this.writeMessage(tag, writePackedFixed64, arr);  },
    writePackedSFixed64: function(tag, arr) { this.writeMessage(tag, writePackedSFixed64, arr); },

    writeBytesField: function(tag, buffer) {
        this.writeTag(tag, Pbf.Bytes);
        this.writeBytes(buffer);
    },
    writeFixed32Field: function(tag, val) {
        this.writeTag(tag, Pbf.Fixed32);
        this.writeFixed32(val);
    },
    writeSFixed32Field: function(tag, val) {
        this.writeTag(tag, Pbf.Fixed32);
        this.writeSFixed32(val);
    },
    writeFixed64Field: function(tag, val) {
        this.writeTag(tag, Pbf.Fixed64);
        this.writeFixed64(val);
    },
    writeSFixed64Field: function(tag, val) {
        this.writeTag(tag, Pbf.Fixed64);
        this.writeSFixed64(val);
    },
    writeVarintField: function(tag, val) {
        this.writeTag(tag, Pbf.Varint);
        this.writeVarint(val);
    },
    writeSVarintField: function(tag, val) {
        this.writeTag(tag, Pbf.Varint);
        this.writeSVarint(val);
    },
    writeStringField: function(tag, str) {
        this.writeTag(tag, Pbf.Bytes);
        this.writeString(str);
    },
    writeFloatField: function(tag, val) {
        this.writeTag(tag, Pbf.Fixed32);
        this.writeFloat(val);
    },
    writeDoubleField: function(tag, val) {
        this.writeTag(tag, Pbf.Fixed64);
        this.writeDouble(val);
    },
    writeBooleanField: function(tag, val) {
        this.writeVarintField(tag, Boolean(val));
    }
};

function readVarintRemainder(val, pbf) {
    var buf = pbf.buf, b;

    b = buf[pbf.pos++]; val += (b & 0x7f) * 0x10000000;         if (b < 0x80) return val;
    b = buf[pbf.pos++]; val += (b & 0x7f) * 0x800000000;        if (b < 0x80) return val;
    b = buf[pbf.pos++]; val += (b & 0x7f) * 0x40000000000;      if (b < 0x80) return val;
    b = buf[pbf.pos++]; val += (b & 0x7f) * 0x2000000000000;    if (b < 0x80) return val;
    b = buf[pbf.pos++]; val += (b & 0x7f) * 0x100000000000000;  if (b < 0x80) return val;
    b = buf[pbf.pos++]; val += (b & 0x7f) * 0x8000000000000000; if (b < 0x80) return val;

    throw new Error('Expected varint not more than 10 bytes');
}

function writeBigVarint(val, pbf) {
    pbf.realloc(10);

    var maxPos = pbf.pos + 10;

    while (val >= 1) {
        if (pbf.pos >= maxPos) throw new Error('Given varint doesn\'t fit into 10 bytes');
        var b = val & 0xff;
        pbf.buf[pbf.pos++] = b | (val >= 0x80 ? 0x80 : 0);
        val /= 0x80;
    }
}

function reallocForRawMessage(startPos, len, pbf) {
    var extraLen =
        len <= 0x3fff ? 1 :
        len <= 0x1fffff ? 2 :
        len <= 0xfffffff ? 3 : Math.ceil(Math.log(len) / (Math.LN2 * 7));

    // if 1 byte isn't enough for encoding message length, shift the data to the right
    pbf.realloc(extraLen);
    for (var i = pbf.pos - 1; i >= startPos; i--) pbf.buf[i + extraLen] = pbf.buf[i];
}

function writePackedVarint(arr, pbf)   { for (var i = 0; i < arr.length; i++) pbf.writeVarint(arr[i]);   }
function writePackedSVarint(arr, pbf)  { for (var i = 0; i < arr.length; i++) pbf.writeSVarint(arr[i]);  }
function writePackedFloat(arr, pbf)    { for (var i = 0; i < arr.length; i++) pbf.writeFloat(arr[i]);    }
function writePackedDouble(arr, pbf)   { for (var i = 0; i < arr.length; i++) pbf.writeDouble(arr[i]);   }
function writePackedBoolean(arr, pbf)  { for (var i = 0; i < arr.length; i++) pbf.writeBoolean(arr[i]);  }
function writePackedFixed32(arr, pbf)  { for (var i = 0; i < arr.length; i++) pbf.writeFixed32(arr[i]);  }
function writePackedSFixed32(arr, pbf) { for (var i = 0; i < arr.length; i++) pbf.writeSFixed32(arr[i]); }
function writePackedFixed64(arr, pbf)  { for (var i = 0; i < arr.length; i++) pbf.writeFixed64(arr[i]);  }
function writePackedSFixed64(arr, pbf) { for (var i = 0; i < arr.length; i++) pbf.writeSFixed64(arr[i]); }

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./buffer":40}],42:[function(require,module,exports){
'use strict';

module.exports = Point;

function Point(x, y) {
    this.x = x;
    this.y = y;
}

Point.prototype = {
    clone: function() { return new Point(this.x, this.y); },

    add:     function(p) { return this.clone()._add(p);     },
    sub:     function(p) { return this.clone()._sub(p);     },
    mult:    function(k) { return this.clone()._mult(k);    },
    div:     function(k) { return this.clone()._div(k);     },
    rotate:  function(a) { return this.clone()._rotate(a);  },
    matMult: function(m) { return this.clone()._matMult(m); },
    unit:    function() { return this.clone()._unit(); },
    perp:    function() { return this.clone()._perp(); },
    round:   function() { return this.clone()._round(); },

    mag: function() {
        return Math.sqrt(this.x * this.x + this.y * this.y);
    },

    equals: function(p) {
        return this.x === p.x &&
               this.y === p.y;
    },

    dist: function(p) {
        return Math.sqrt(this.distSqr(p));
    },

    distSqr: function(p) {
        var dx = p.x - this.x,
            dy = p.y - this.y;
        return dx * dx + dy * dy;
    },

    angle: function() {
        return Math.atan2(this.y, this.x);
    },

    angleTo: function(b) {
        return Math.atan2(this.y - b.y, this.x - b.x);
    },

    angleWith: function(b) {
        return this.angleWithSep(b.x, b.y);
    },

    // Find the angle of the two vectors, solving the formula for the cross product a x b = |a||b|sin(θ) for θ.
    angleWithSep: function(x, y) {
        return Math.atan2(
            this.x * y - this.y * x,
            this.x * x + this.y * y);
    },

    _matMult: function(m) {
        var x = m[0] * this.x + m[1] * this.y,
            y = m[2] * this.x + m[3] * this.y;
        this.x = x;
        this.y = y;
        return this;
    },

    _add: function(p) {
        this.x += p.x;
        this.y += p.y;
        return this;
    },

    _sub: function(p) {
        this.x -= p.x;
        this.y -= p.y;
        return this;
    },

    _mult: function(k) {
        this.x *= k;
        this.y *= k;
        return this;
    },

    _div: function(k) {
        this.x /= k;
        this.y /= k;
        return this;
    },

    _unit: function() {
        this._div(this.mag());
        return this;
    },

    _perp: function() {
        var y = this.y;
        this.y = this.x;
        this.x = -y;
        return this;
    },

    _rotate: function(angle) {
        var cos = Math.cos(angle),
            sin = Math.sin(angle),
            x = cos * this.x - sin * this.y,
            y = sin * this.x + cos * this.y;
        this.x = x;
        this.y = y;
        return this;
    },

    _round: function() {
        this.x = Math.round(this.x);
        this.y = Math.round(this.y);
        return this;
    }
};

// constructs Point from an array if necessary
Point.convert = function (a) {
    if (a instanceof Point) {
        return a;
    }
    if (Array.isArray(a)) {
        return new Point(a[0], a[1]);
    }
    return a;
};

},{}],43:[function(require,module,exports){
'use strict';

var kdbush = require('kdbush');

module.exports = supercluster;

function supercluster(options) {
    return new SuperCluster(options);
}

function SuperCluster(options) {
    this.options = extend(Object.create(this.options), options);
    this.trees = new Array(this.options.maxZoom + 1);
}

SuperCluster.prototype = {
    options: {
        minZoom: 0,   // min zoom to generate clusters on
        maxZoom: 16,  // max zoom level to cluster the points on
        radius: 40,   // cluster radius in pixels
        extent: 512,  // tile extent (radius is calculated relative to it)
        nodeSize: 64, // size of the KD-tree leaf node, affects performance
        log: false,   // whether to log timing info

        // a reduce function for calculating custom cluster properties
        reduce: null, // function (accumulated, props) { accumulated.sum += props.sum; }

        // initial properties of a cluster (before running the reducer)
        initial: function () { return {}; }, // function () { return {sum: 0}; },

        // properties to use for individual points when running the reducer
        map: function (props) { return props; } // function (props) { return {sum: props.my_value}; },
    },

    load: function (points) {
        var log = this.options.log;

        if (log) console.time('total time');

        var timerId = 'prepare ' + points.length + ' points';
        if (log) console.time(timerId);

        this.points = points;

        // generate a cluster object for each point
        var clusters = points.map(createPointCluster);
        if (log) console.timeEnd(timerId);

        // cluster points on max zoom, then cluster the results on previous zoom, etc.;
        // results in a cluster hierarchy across zoom levels
        for (var z = this.options.maxZoom; z >= this.options.minZoom; z--) {
            var now = +Date.now();

            // index input points into a KD-tree
            this.trees[z + 1] = kdbush(clusters, getX, getY, this.options.nodeSize, Float32Array);

            clusters = this._cluster(clusters, z); // create a new set of clusters for the zoom

            if (log) console.log('z%d: %d clusters in %dms', z, clusters.length, +Date.now() - now);
        }

        // index top-level clusters
        this.trees[this.options.minZoom] = kdbush(clusters, getX, getY, this.options.nodeSize, Float32Array);

        if (log) console.timeEnd('total time');

        return this;
    },

    getClusters: function (bbox, zoom) {
        var tree = this.trees[this._limitZoom(zoom)];
        var ids = tree.range(lngX(bbox[0]), latY(bbox[3]), lngX(bbox[2]), latY(bbox[1]));
        var clusters = [];
        for (var i = 0; i < ids.length; i++) {
            var c = tree.points[ids[i]];
            clusters.push(c.numPoints ? getClusterJSON(c) : this.points[c.id]);
        }
        return clusters;
    },

    getChildren: function (clusterId, clusterZoom) {
        var origin = this.trees[clusterZoom + 1].points[clusterId];
        var r = this.options.radius / (this.options.extent * Math.pow(2, clusterZoom));
        var points = this.trees[clusterZoom + 1].within(origin.x, origin.y, r);
        var children = [];
        for (var i = 0; i < points.length; i++) {
            var c = this.trees[clusterZoom + 1].points[points[i]];
            if (c.parentId === clusterId) {
                children.push(c.numPoints ? getClusterJSON(c) : this.points[c.id]);
            }
        }
        return children;
    },

    getLeaves: function (clusterId, clusterZoom, limit, offset) {
        limit = limit || 10;
        offset = offset || 0;

        var leaves = [];
        this._appendLeaves(leaves, clusterId, clusterZoom, limit, offset, 0);

        return leaves;
    },

    getTile: function (z, x, y) {
        var tree = this.trees[this._limitZoom(z)];
        var z2 = Math.pow(2, z);
        var extent = this.options.extent;
        var r = this.options.radius;
        var p = r / extent;
        var top = (y - p) / z2;
        var bottom = (y + 1 + p) / z2;

        var tile = {
            features: []
        };

        this._addTileFeatures(
            tree.range((x - p) / z2, top, (x + 1 + p) / z2, bottom),
            tree.points, x, y, z2, tile);

        if (x === 0) {
            this._addTileFeatures(
                tree.range(1 - p / z2, top, 1, bottom),
                tree.points, z2, y, z2, tile);
        }
        if (x === z2 - 1) {
            this._addTileFeatures(
                tree.range(0, top, p / z2, bottom),
                tree.points, -1, y, z2, tile);
        }

        return tile.features.length ? tile : null;
    },

    getClusterExpansionZoom: function (clusterId, clusterZoom) {
        while (clusterZoom < this.options.maxZoom) {
            var children = this.getChildren(clusterId, clusterZoom);
            clusterZoom++;
            if (children.length !== 1) break;
            clusterId = children[0].properties.cluster_id;
        }
        return clusterZoom;
    },

    _appendLeaves: function (result, clusterId, clusterZoom, limit, offset, skipped) {
        var children = this.getChildren(clusterId, clusterZoom);

        for (var i = 0; i < children.length; i++) {
            var props = children[i].properties;

            if (props.cluster) {
                if (skipped + props.point_count <= offset) {
                    // skip the whole cluster
                    skipped += props.point_count;
                } else {
                    // enter the cluster
                    skipped = this._appendLeaves(
                        result, props.cluster_id, clusterZoom + 1, limit, offset, skipped);
                    // exit the cluster
                }
            } else if (skipped < offset) {
                // skip a single point
                skipped++;
            } else {
                // add a single point
                result.push(children[i]);
            }
            if (result.length === limit) break;
        }

        return skipped;
    },

    _addTileFeatures: function (ids, points, x, y, z2, tile) {
        for (var i = 0; i < ids.length; i++) {
            var c = points[ids[i]];
            tile.features.push({
                type: 1,
                geometry: [[
                    Math.round(this.options.extent * (c.x * z2 - x)),
                    Math.round(this.options.extent * (c.y * z2 - y))
                ]],
                tags: c.numPoints ? getClusterProperties(c) : this.points[c.id].properties
            });
        }
    },

    _limitZoom: function (z) {
        return Math.max(this.options.minZoom, Math.min(z, this.options.maxZoom + 1));
    },

    _cluster: function (points, zoom) {
        var clusters = [];
        var r = this.options.radius / (this.options.extent * Math.pow(2, zoom));

        // loop through each point
        for (var i = 0; i < points.length; i++) {
            var p = points[i];
            // if we've already visited the point at this zoom level, skip it
            if (p.zoom <= zoom) continue;
            p.zoom = zoom;

            // find all nearby points
            var tree = this.trees[zoom + 1];
            var neighborIds = tree.within(p.x, p.y, r);

            var numPoints = p.numPoints || 1;
            var wx = p.x * numPoints;
            var wy = p.y * numPoints;

            var clusterProperties = null;

            if (this.options.reduce) {
                clusterProperties = this.options.initial();
                this._accumulate(clusterProperties, p);
            }

            for (var j = 0; j < neighborIds.length; j++) {
                var b = tree.points[neighborIds[j]];
                // filter out neighbors that are too far or already processed
                if (zoom < b.zoom) {
                    var numPoints2 = b.numPoints || 1;
                    b.zoom = zoom; // save the zoom (so it doesn't get processed twice)
                    wx += b.x * numPoints2; // accumulate coordinates for calculating weighted center
                    wy += b.y * numPoints2;
                    numPoints += numPoints2;
                    b.parentId = i;

                    if (this.options.reduce) {
                        this._accumulate(clusterProperties, b);
                    }
                }
            }

            if (numPoints === 1) {
                clusters.push(p);
            } else {
                p.parentId = i;
                clusters.push(createCluster(wx / numPoints, wy / numPoints, numPoints, i, clusterProperties));
            }
        }

        return clusters;
    },

    _accumulate: function (clusterProperties, point) {
        var properties = point.numPoints ?
            point.properties :
            this.options.map(this.points[point.id].properties);

        this.options.reduce(clusterProperties, properties);
    }
};

function createCluster(x, y, numPoints, id, properties) {
    return {
        x: x, // weighted cluster center
        y: y,
        zoom: Infinity, // the last zoom the cluster was processed at
        id: id, // index of the first child of the cluster in the zoom level tree
        properties: properties,
        parentId: -1, // parent cluster id
        numPoints: numPoints
    };
}

function createPointCluster(p, id) {
    var coords = p.geometry.coordinates;
    return {
        x: lngX(coords[0]), // projected point coordinates
        y: latY(coords[1]),
        zoom: Infinity, // the last zoom the point was processed at
        id: id, // index of the source feature in the original input array
        parentId: -1 // parent cluster id
    };
}

function getClusterJSON(cluster) {
    return {
        type: 'Feature',
        properties: getClusterProperties(cluster),
        geometry: {
            type: 'Point',
            coordinates: [xLng(cluster.x), yLat(cluster.y)]
        }
    };
}

function getClusterProperties(cluster) {
    var count = cluster.numPoints;
    var abbrev = count >= 10000 ? Math.round(count / 1000) + 'k' :
                 count >= 1000 ? (Math.round(count / 100) / 10) + 'k' : count;
    return extend(extend({}, cluster.properties), {
        cluster: true,
        cluster_id: cluster.id,
        point_count: count,
        point_count_abbreviated: abbrev
    });
}

// longitude/latitude to spherical mercator in [0..1] range
function lngX(lng) {
    return lng / 360 + 0.5;
}
function latY(lat) {
    var sin = Math.sin(lat * Math.PI / 180),
        y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);
    return y < 0 ? 0 :
           y > 1 ? 1 : y;
}

// spherical mercator to longitude/latitude
function xLng(x) {
    return (x - 0.5) * 360;
}
function yLat(y) {
    var y2 = (180 - y * 360) * Math.PI / 180;
    return 360 * Math.atan(Math.exp(y2)) / Math.PI - 90;
}

function extend(dest, src) {
    for (var id in src) dest[id] = src[id];
    return dest;
}

function getX(p) {
    return p.x;
}
function getY(p) {
    return p.y;
}

},{"kdbush":16}],44:[function(require,module,exports){
module.exports.VectorTile = require('./lib/vectortile.js');
module.exports.VectorTileFeature = require('./lib/vectortilefeature.js');
module.exports.VectorTileLayer = require('./lib/vectortilelayer.js');

},{"./lib/vectortile.js":45,"./lib/vectortilefeature.js":46,"./lib/vectortilelayer.js":47}],45:[function(require,module,exports){
'use strict';

var VectorTileLayer = require('./vectortilelayer');

module.exports = VectorTile;

function VectorTile(pbf, end) {
    this.layers = pbf.readFields(readTile, {}, end);
}

function readTile(tag, layers, pbf) {
    if (tag === 3) {
        var layer = new VectorTileLayer(pbf, pbf.readVarint() + pbf.pos);
        if (layer.length) layers[layer.name] = layer;
    }
}


},{"./vectortilelayer":47}],46:[function(require,module,exports){
'use strict';

var Point = require('point-geometry');

module.exports = VectorTileFeature;

function VectorTileFeature(pbf, end, extent, keys, values) {
    // Public
    this.properties = {};
    this.extent = extent;
    this.type = 0;

    // Private
    this._pbf = pbf;
    this._geometry = -1;
    this._keys = keys;
    this._values = values;

    pbf.readFields(readFeature, this, end);
}

function readFeature(tag, feature, pbf) {
    if (tag == 1) feature.id = pbf.readVarint();
    else if (tag == 2) readTag(pbf, feature);
    else if (tag == 3) feature.type = pbf.readVarint();
    else if (tag == 4) feature._geometry = pbf.pos;
}

function readTag(pbf, feature) {
    var end = pbf.readVarint() + pbf.pos;

    while (pbf.pos < end) {
        var key = feature._keys[pbf.readVarint()],
            value = feature._values[pbf.readVarint()];
        feature.properties[key] = value;
    }
}

VectorTileFeature.types = ['Unknown', 'Point', 'LineString', 'Polygon'];

VectorTileFeature.prototype.loadGeometry = function() {
    var pbf = this._pbf;
    pbf.pos = this._geometry;

    var end = pbf.readVarint() + pbf.pos,
        cmd = 1,
        length = 0,
        x = 0,
        y = 0,
        lines = [],
        line;

    while (pbf.pos < end) {
        if (!length) {
            var cmdLen = pbf.readVarint();
            cmd = cmdLen & 0x7;
            length = cmdLen >> 3;
        }

        length--;

        if (cmd === 1 || cmd === 2) {
            x += pbf.readSVarint();
            y += pbf.readSVarint();

            if (cmd === 1) { // moveTo
                if (line) lines.push(line);
                line = [];
            }

            line.push(new Point(x, y));

        } else if (cmd === 7) {

            // Workaround for https://github.com/mapbox/mapnik-vector-tile/issues/90
            if (line) {
                line.push(line[0].clone()); // closePolygon
            }

        } else {
            throw new Error('unknown command ' + cmd);
        }
    }

    if (line) lines.push(line);

    return lines;
};

VectorTileFeature.prototype.bbox = function() {
    var pbf = this._pbf;
    pbf.pos = this._geometry;

    var end = pbf.readVarint() + pbf.pos,
        cmd = 1,
        length = 0,
        x = 0,
        y = 0,
        x1 = Infinity,
        x2 = -Infinity,
        y1 = Infinity,
        y2 = -Infinity;

    while (pbf.pos < end) {
        if (!length) {
            var cmdLen = pbf.readVarint();
            cmd = cmdLen & 0x7;
            length = cmdLen >> 3;
        }

        length--;

        if (cmd === 1 || cmd === 2) {
            x += pbf.readSVarint();
            y += pbf.readSVarint();
            if (x < x1) x1 = x;
            if (x > x2) x2 = x;
            if (y < y1) y1 = y;
            if (y > y2) y2 = y;

        } else if (cmd !== 7) {
            throw new Error('unknown command ' + cmd);
        }
    }

    return [x1, y1, x2, y2];
};

VectorTileFeature.prototype.toGeoJSON = function(x, y, z) {
    var size = this.extent * Math.pow(2, z),
        x0 = this.extent * x,
        y0 = this.extent * y,
        coords = this.loadGeometry(),
        type = VectorTileFeature.types[this.type],
        i, j;

    function project(line) {
        for (var j = 0; j < line.length; j++) {
            var p = line[j], y2 = 180 - (p.y + y0) * 360 / size;
            line[j] = [
                (p.x + x0) * 360 / size - 180,
                360 / Math.PI * Math.atan(Math.exp(y2 * Math.PI / 180)) - 90
            ];
        }
    }

    switch (this.type) {
    case 1:
        var points = [];
        for (i = 0; i < coords.length; i++) {
            points[i] = coords[i][0];
        }
        coords = points;
        project(coords);
        break;

    case 2:
        for (i = 0; i < coords.length; i++) {
            project(coords[i]);
        }
        break;

    case 3:
        coords = classifyRings(coords);
        for (i = 0; i < coords.length; i++) {
            for (j = 0; j < coords[i].length; j++) {
                project(coords[i][j]);
            }
        }
        break;
    }

    if (coords.length === 1) {
        coords = coords[0];
    } else {
        type = 'Multi' + type;
    }

    var result = {
        type: "Feature",
        geometry: {
            type: type,
            coordinates: coords
        },
        properties: this.properties
    };

    if ('id' in this) {
        result.id = this.id;
    }

    return result;
};

// classifies an array of rings into polygons with outer rings and holes

function classifyRings(rings) {
    var len = rings.length;

    if (len <= 1) return [rings];

    var polygons = [],
        polygon,
        ccw;

    for (var i = 0; i < len; i++) {
        var area = signedArea(rings[i]);
        if (area === 0) continue;

        if (ccw === undefined) ccw = area < 0;

        if (ccw === area < 0) {
            if (polygon) polygons.push(polygon);
            polygon = [rings[i]];

        } else {
            polygon.push(rings[i]);
        }
    }
    if (polygon) polygons.push(polygon);

    return polygons;
}

function signedArea(ring) {
    var sum = 0;
    for (var i = 0, len = ring.length, j = len - 1, p1, p2; i < len; j = i++) {
        p1 = ring[i];
        p2 = ring[j];
        sum += (p2.x - p1.x) * (p1.y + p2.y);
    }
    return sum;
}

},{"point-geometry":42}],47:[function(require,module,exports){
'use strict';

var VectorTileFeature = require('./vectortilefeature.js');

module.exports = VectorTileLayer;

function VectorTileLayer(pbf, end) {
    // Public
    this.version = 1;
    this.name = null;
    this.extent = 4096;
    this.length = 0;

    // Private
    this._pbf = pbf;
    this._keys = [];
    this._values = [];
    this._features = [];

    pbf.readFields(readLayer, this, end);

    this.length = this._features.length;
}

function readLayer(tag, layer, pbf) {
    if (tag === 15) layer.version = pbf.readVarint();
    else if (tag === 1) layer.name = pbf.readString();
    else if (tag === 5) layer.extent = pbf.readVarint();
    else if (tag === 2) layer._features.push(pbf.pos);
    else if (tag === 3) layer._keys.push(pbf.readString());
    else if (tag === 4) layer._values.push(readValueMessage(pbf));
}

function readValueMessage(pbf) {
    var value = null,
        end = pbf.readVarint() + pbf.pos;

    while (pbf.pos < end) {
        var tag = pbf.readVarint() >> 3;

        value = tag === 1 ? pbf.readString() :
            tag === 2 ? pbf.readFloat() :
            tag === 3 ? pbf.readDouble() :
            tag === 4 ? pbf.readVarint64() :
            tag === 5 ? pbf.readVarint() :
            tag === 6 ? pbf.readSVarint() :
            tag === 7 ? pbf.readBoolean() : null;
    }

    return value;
}

// return feature `i` from this layer as a `VectorTileFeature`
VectorTileLayer.prototype.feature = function(i) {
    if (i < 0 || i >= this._features.length) throw new Error('feature index out of bounds');

    this._pbf.pos = this._features[i];

    var end = this._pbf.readVarint() + this._pbf.pos;
    return new VectorTileFeature(this._pbf, end, this.extent, this._keys, this._values);
};

},{"./vectortilefeature.js":46}],48:[function(require,module,exports){
var Pbf = require('pbf')
var vtpb = require('./vector-tile-pb')
var GeoJSONWrapper = require('./lib/geojson_wrapper')

module.exports = fromVectorTileJs
module.exports.fromVectorTileJs = fromVectorTileJs
module.exports.fromGeojsonVt = fromGeojsonVt
module.exports.GeoJSONWrapper = GeoJSONWrapper

/**
 * Serialize a vector-tile-js-created tile to pbf
 *
 * @param {Object} tile
 * @return {Buffer} uncompressed, pbf-serialized tile data
 */
function fromVectorTileJs (tile) {
  var layers = []
  for (var l in tile.layers) {
    layers.push(prepareLayer(tile.layers[l]))
  }

  var out = new Pbf()
  vtpb.tile.write({ layers: layers }, out)
  return out.finish()
}

/**
 * Serialized a geojson-vt-created tile to pbf.
 *
 * @param {Object} layers - An object mapping layer names to geojson-vt-created vector tile objects
 * @return {Buffer} uncompressed, pbf-serialized tile data
 */
function fromGeojsonVt (layers) {
  var l = {}
  for (var k in layers) {
    l[k] = new GeoJSONWrapper(layers[k].features)
    l[k].name = k
  }
  return fromVectorTileJs({layers: l})
}

/**
 * Prepare the given layer to be serialized by the auto-generated pbf
 * serializer by encoding the feature geometry and properties.
 */
function prepareLayer (layer) {
  var preparedLayer = {
    name: layer.name || '',
    version: layer.version || 1,
    extent: layer.extent || 4096,
    keys: [],
    values: [],
    features: []
  }

  var keycache = {}
  var valuecache = {}

  for (var i = 0; i < layer.length; i++) {
    var feature = layer.feature(i)
    feature.geometry = encodeGeometry(feature.loadGeometry())

    var tags = []
    for (var key in feature.properties) {
      var keyIndex = keycache[key]
      if (typeof keyIndex === 'undefined') {
        preparedLayer.keys.push(key)
        keyIndex = preparedLayer.keys.length - 1
        keycache[key] = keyIndex
      }
      var value = wrapValue(feature.properties[key])
      var valueIndex = valuecache[value.key]
      if (typeof valueIndex === 'undefined') {
        preparedLayer.values.push(value)
        valueIndex = preparedLayer.values.length - 1
        valuecache[value.key] = valueIndex
      }
      tags.push(keyIndex)
      tags.push(valueIndex)
    }

    feature.tags = tags
    preparedLayer.features.push(feature)
  }

  return preparedLayer
}

function command (cmd, length) {
  return (length << 3) + (cmd & 0x7)
}

function zigzag (num) {
  return (num << 1) ^ (num >> 31)
}

/**
 * Encode a polygon's geometry into an array ready to be serialized
 * to mapbox vector tile specified geometry data.
 *
 * @param {Array} Rings, each being an array of [x, y] tile-space coordinates
 * @return {Array} encoded geometry
 */
function encodeGeometry (geometry) {
  var encoded = []
  var x = 0
  var y = 0
  var rings = geometry.length
  for (var r = 0; r < rings; r++) {
    var ring = geometry[r]
    encoded.push(command(1, 1)) // moveto
    for (var i = 0; i < ring.length; i++) {
      if (i === 1) {
        encoded.push(command(2, ring.length - 1)) // lineto
      }
      var dx = ring[i].x - x
      var dy = ring[i].y - y
      encoded.push(zigzag(dx), zigzag(dy))
      x += dx
      y += dy
    }
  }

  return encoded
}

/**
 * Wrap a property value according to its type. The returned object
 * is of the form { xxxx_value: primitiveValue }, which is what the generated
 * protobuf serializer expects.
 */
function wrapValue (value) {
  var result
  var type = typeof value
  if (type === 'string') {
    result = { string_value: value }
  } else if (type === 'boolean') {
    result = { bool_value: value }
  } else if (type === 'number') {
    if (value % 1 !== 0) {
      result = { double_value: value }
    } else if (value < 0) {
      result = { sint_value: value }
    } else {
      result = { uint_value: value }
    }
  } else {
    value = JSON.stringify(value)
    result = { string_value: value }
  }

  result.key = type + ':' + value
  return result
}

},{"./lib/geojson_wrapper":49,"./vector-tile-pb":50,"pbf":41}],49:[function(require,module,exports){
'use strict'

var Point = require('point-geometry')
var VectorTileFeature = require('vector-tile').VectorTileFeature

module.exports = GeoJSONWrapper

// conform to vectortile api
function GeoJSONWrapper (features) {
  this.features = features
  this.length = features.length
}

GeoJSONWrapper.prototype.feature = function (i) {
  return new FeatureWrapper(this.features[i])
}

function FeatureWrapper (feature) {
  this.id = typeof feature.id === 'number' ? feature.id : undefined
  this.type = feature.type
  this.rawGeometry = feature.type === 1 ? [feature.geometry] : feature.geometry
  this.properties = feature.tags
  this.extent = 4096
}

FeatureWrapper.prototype.loadGeometry = function () {
  var rings = this.rawGeometry
  this.geometry = []

  for (var i = 0; i < rings.length; i++) {
    var ring = rings[i]
    var newRing = []
    for (var j = 0; j < ring.length; j++) {
      newRing.push(new Point(ring[j][0], ring[j][1]))
    }
    this.geometry.push(newRing)
  }
  return this.geometry
}

FeatureWrapper.prototype.bbox = function () {
  if (!this.geometry) this.loadGeometry()

  var rings = this.geometry
  var x1 = Infinity
  var x2 = -Infinity
  var y1 = Infinity
  var y2 = -Infinity

  for (var i = 0; i < rings.length; i++) {
    var ring = rings[i]

    for (var j = 0; j < ring.length; j++) {
      var coord = ring[j]

      x1 = Math.min(x1, coord.x)
      x2 = Math.max(x2, coord.x)
      y1 = Math.min(y1, coord.y)
      y2 = Math.max(y2, coord.y)
    }
  }

  return [x1, y1, x2, y2]
}

FeatureWrapper.prototype.toGeoJSON = VectorTileFeature.prototype.toGeoJSON

},{"point-geometry":42,"vector-tile":44}],50:[function(require,module,exports){
'use strict';

// tile ========================================

var tile = exports.tile = {read: readTile, write: writeTile};

tile.GeomType = {
    "Unknown": 0,
    "Point": 1,
    "LineString": 2,
    "Polygon": 3
};

function readTile(pbf, end) {
    return pbf.readFields(readTileField, {"layers": []}, end);
}

function readTileField(tag, tile, pbf) {
    if (tag === 3) tile.layers.push(readLayer(pbf, pbf.readVarint() + pbf.pos));
}

function writeTile(tile, pbf) {
    var i;
    if (tile.layers !== undefined) for (i = 0; i < tile.layers.length; i++) pbf.writeMessage(3, writeLayer, tile.layers[i]);
}

// value ========================================

tile.value = {read: readValue, write: writeValue};

function readValue(pbf, end) {
    return pbf.readFields(readValueField, {}, end);
}

function readValueField(tag, value, pbf) {
    if (tag === 1) value.string_value = pbf.readString();
    else if (tag === 2) value.float_value = pbf.readFloat();
    else if (tag === 3) value.double_value = pbf.readDouble();
    else if (tag === 4) value.int_value = pbf.readVarint();
    else if (tag === 5) value.uint_value = pbf.readVarint();
    else if (tag === 6) value.sint_value = pbf.readSVarint();
    else if (tag === 7) value.bool_value = pbf.readBoolean();
}

function writeValue(value, pbf) {
    if (value.string_value !== undefined) pbf.writeStringField(1, value.string_value);
    if (value.float_value !== undefined) pbf.writeFloatField(2, value.float_value);
    if (value.double_value !== undefined) pbf.writeDoubleField(3, value.double_value);
    if (value.int_value !== undefined) pbf.writeVarintField(4, value.int_value);
    if (value.uint_value !== undefined) pbf.writeVarintField(5, value.uint_value);
    if (value.sint_value !== undefined) pbf.writeSVarintField(6, value.sint_value);
    if (value.bool_value !== undefined) pbf.writeBooleanField(7, value.bool_value);
}

// feature ========================================

tile.feature = {read: readFeature, write: writeFeature};

function readFeature(pbf, end) {
    var feature = pbf.readFields(readFeatureField, {}, end);
    if (feature.type === undefined) feature.type = "Unknown";
    return feature;
}

function readFeatureField(tag, feature, pbf) {
    if (tag === 1) feature.id = pbf.readVarint();
    else if (tag === 2) feature.tags = pbf.readPackedVarint();
    else if (tag === 3) feature.type = pbf.readVarint();
    else if (tag === 4) feature.geometry = pbf.readPackedVarint();
}

function writeFeature(feature, pbf) {
    if (feature.id !== undefined) pbf.writeVarintField(1, feature.id);
    if (feature.tags !== undefined) pbf.writePackedVarint(2, feature.tags);
    if (feature.type !== undefined) pbf.writeVarintField(3, feature.type);
    if (feature.geometry !== undefined) pbf.writePackedVarint(4, feature.geometry);
}

// layer ========================================

tile.layer = {read: readLayer, write: writeLayer};

function readLayer(pbf, end) {
    return pbf.readFields(readLayerField, {"features": [], "keys": [], "values": []}, end);
}

function readLayerField(tag, layer, pbf) {
    if (tag === 15) layer.version = pbf.readVarint();
    else if (tag === 1) layer.name = pbf.readString();
    else if (tag === 2) layer.features.push(readFeature(pbf, pbf.readVarint() + pbf.pos));
    else if (tag === 3) layer.keys.push(pbf.readString());
    else if (tag === 4) layer.values.push(readValue(pbf, pbf.readVarint() + pbf.pos));
    else if (tag === 5) layer.extent = pbf.readVarint();
}

function writeLayer(layer, pbf) {
    if (layer.version !== undefined) pbf.writeVarintField(15, layer.version);
    if (layer.name !== undefined) pbf.writeStringField(1, layer.name);
    var i;
    if (layer.features !== undefined) for (i = 0; i < layer.features.length; i++) pbf.writeMessage(2, writeFeature, layer.features[i]);
    if (layer.keys !== undefined) for (i = 0; i < layer.keys.length; i++) pbf.writeStringField(3, layer.keys[i]);
    if (layer.values !== undefined) for (i = 0; i < layer.values.length; i++) pbf.writeMessage(4, writeValue, layer.values[i]);
    if (layer.extent !== undefined) pbf.writeVarintField(5, layer.extent);
}

},{}],51:[function(require,module,exports){
var bundleFn = arguments[3];
var sources = arguments[4];
var cache = arguments[5];

var stringify = JSON.stringify;

module.exports = function (fn, options) {
    var wkey;
    var cacheKeys = Object.keys(cache);

    for (var i = 0, l = cacheKeys.length; i < l; i++) {
        var key = cacheKeys[i];
        var exp = cache[key].exports;
        // Using babel as a transpiler to use esmodule, the export will always
        // be an object with the default export as a property of it. To ensure
        // the existing api and babel esmodule exports are both supported we
        // check for both
        if (exp === fn || exp && exp.default === fn) {
            wkey = key;
            break;
        }
    }

    if (!wkey) {
        wkey = Math.floor(Math.pow(16, 8) * Math.random()).toString(16);
        var wcache = {};
        for (var i = 0, l = cacheKeys.length; i < l; i++) {
            var key = cacheKeys[i];
            wcache[key] = key;
        }
        sources[wkey] = [
            Function(['require','module','exports'], '(' + fn + ')(self)'),
            wcache
        ];
    }
    var skey = Math.floor(Math.pow(16, 8) * Math.random()).toString(16);

    var scache = {}; scache[wkey] = wkey;
    sources[skey] = [
        Function(['require'], (
            // try to call default if defined to also support babel esmodule
            // exports
            'var f = require(' + stringify(wkey) + ');' +
            '(f.default ? f.default : f)(self);'
        )),
        scache
    ];

    var workerSources = {};
    resolveSources(skey);

    function resolveSources(key) {
        workerSources[key] = true;

        for (var depPath in sources[key][1]) {
            var depKey = sources[key][1][depPath];
            if (!workerSources[depKey]) {
                resolveSources(depKey);
            }
        }
    }

    var src = '(' + bundleFn + ')({'
        + Object.keys(workerSources).map(function (key) {
            return stringify(key) + ':['
                + sources[key][0]
                + ',' + stringify(sources[key][1]) + ']'
            ;
        }).join(',')
        + '},{},[' + stringify(skey) + '])'
    ;

    var URL = window.URL || window.webkitURL || window.mozURL || window.msURL;

    var blob = new Blob([src], { type: 'text/javascript' });
    if (options && options.bare) { return blob; }
    var workerUrl = URL.createObjectURL(blob);
    var worker = new Worker(workerUrl);
    worker.objectURL = workerUrl;
    return worker;
};

},{}],52:[function(require,module,exports){
module.exports.RADIUS = 6378137;
module.exports.FLATTENING = 1/298.257223563;
module.exports.POLAR_RADIUS = 6356752.3142;

},{}]},{},[1])(1)
});