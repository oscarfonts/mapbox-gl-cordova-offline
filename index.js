var mapboxgl = require('mapbox-gl/src/index');
var CustomGeoJSONSource = require('./custom_geojson_source');

mapboxgl.SourceTypes = {
    CustomGeoJSON: CustomGeoJSONSource
}

module.exports = mapboxgl;
