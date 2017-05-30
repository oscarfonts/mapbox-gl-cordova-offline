const VectorTileWorkerSource = require('mapbox-gl/src/source/vector_tile_worker_source');

function customLoadVectorData (options, callback) {
  console.log("Called customLoadVectorData with x, y, z =", options.coord.x, options.coord.y, options.coord.z);
  VectorTileWorkerSource.prototype.loadVectorData.call(this, options, function (err, data) {
    // Read the data normally, then call this custom callback, which shows options and data in console
    if (err) { return callback(err) }
    console.log(data);
    callback(null, data);
  })
}

class CustomVectorTileWorkerSource extends VectorTileWorkerSource {
  constructor(actor, layerIndex) {
      super(actor, layerIndex, customLoadVectorData);
  }
}

module.exports = function(self) {
  self.registerWorkerSource('custom_vector', CustomVectorTileWorkerSource)
};
