const GeoJSONWorkerSource = require('mapbox-gl/src/source/geojson_worker_source');

function customLoadGeoJSON (options, callback) {
  GeoJSONWorkerSource.prototype.loadGeoJSON.call(this, options, function (err, data) {
    // Read the data normally, then call this custom callback, which shows options and data in console
    if (err) { return callback(err) }
    console.log(options, data);
    callback(null, data);
  })
}

class CustomGeoJSONWorkerSource extends GeoJSONWorkerSource {
  constructor(actor, style) {
      super(actor, style, customLoadGeoJSON);
  }
}

module.exports = function(self) {
  self.registerWorkerSource('custom_geojson', CustomGeoJSONWorkerSource)
};
