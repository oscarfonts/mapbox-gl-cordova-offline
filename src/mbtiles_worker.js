const pako = require('pako/lib/inflate');
const Protobuf = require('pbf');
const VectorTileWorkerSource = require('mapbox-gl/src/source/vector_tile_worker_source');
const vt = require('vector-tile');
const base64js = require('base64-js');


class MBTilesWorkerSource extends VectorTileWorkerSource {
    loadVectorData(params, callback) {
        let rawData = params.base64Data ? pako.inflate(this.convertBase64ToUint8Array(params.base64Data)) : undefined;
        const vectorTile = new vt.VectorTile(new Protobuf(rawData));
        vectorTile.rawData = rawData;
        vectorTile.cacheControl = null;
        vectorTile.expires = null;
        callback(undefined, vectorTile);
    }

    convertBase64ToUint8Array(base64) {
        return base64js.toByteArray(base64)
    }
}

module.exports = function (self) {
    self.registerWorkerSource('mbtiles', MBTilesWorkerSource)
};
