var GeoJSONSource = require('mapbox-gl/src/source/geojson_source');
var webworkify = require('webworkify');

class CustomGeoJSONSource extends GeoJSONSource {
}

CustomGeoJSONSource.workerSourceURL = URL.createObjectURL(webworkify(require('./custom_geojson_worker.js'), {bare: true}));

module.exports = CustomGeoJSONSource;
