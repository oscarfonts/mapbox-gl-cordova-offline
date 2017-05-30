const GeoJSONSource = require('mapbox-gl/src/source/geojson_source');
const webworkify = require('webworkify');

class CustomGeoJSONSource extends GeoJSONSource {
    constructor(id, options, dispatcher, eventedParent) {
        super(id, options, dispatcher, eventedParent);
        this.type="custom_geojson";
    }
}

CustomGeoJSONSource.workerSourceURL = URL.createObjectURL(webworkify(require('./custom_geojson_worker.js'), {bare: true}));

module.exports = CustomGeoJSONSource;
