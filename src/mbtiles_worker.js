const VectorTileWorkerSource = require('mapbox-gl/src/source/vector_tile_worker_source');

module.exports = function (self) {
    self.registerWorkerSource('mbtiles', VectorTileWorkerSource)
};
