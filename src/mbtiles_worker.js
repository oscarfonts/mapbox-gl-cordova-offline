const vt = require('vector-tile');
const Protobuf = require('pbf');
const VectorTileWorkerSource = require('mapbox-gl/src/source/vector_tile_worker_source');

class MBTilesWorkerSource extends VectorTileWorkerSource {
    loadVectorData(params, callback) {
        const vectorTile = new vt.VectorTile(new Protobuf(params.vectorTile.data));
        vectorTile.rawData = params.vectorTile.data;
        vectorTile.cacheControl = null;
        vectorTile.expires = null;
        callback(undefined, vectorTile);
    }
}

module.exports = function (self) {
    self.registerWorkerSource('mbtiles', MBTilesWorkerSource)
};
