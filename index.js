const mapboxgl = require('mapbox-gl/src/index');
const CustomGeoJSONSource = require('./custom_geojson_source');
const CustomVectorTileSource = require('./custom_vector_tile_source');

mapboxgl.SourceTypes = {
    CustomGeoJSON: CustomGeoJSONSource,
    CustomVectorTile: CustomVectorTileSource
};

module.exports = mapboxgl;
