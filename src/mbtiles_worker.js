const pako = require('pako/lib/inflate');
const Protobuf = require('pbf');
const VectorTileWorkerSource = require('mapbox-gl/src/source/vector_tile_worker_source');
const vt = require('vector-tile');


class MBTilesWorkerSource extends VectorTileWorkerSource {
    loadVectorData(params, callback) {
        var rawData = params.base64Data ? pako.inflate(this.convertBase64ToUint8Array(params.base64Data)) : undefined;
        const vectorTile = new vt.VectorTile(new Protobuf(rawData));
        vectorTile.rawData = rawData;
        vectorTile.cacheControl = null;
        vectorTile.expires = null;
        callback(undefined, vectorTile);
    }

    convertBase64ToUint8Array(base64) {
        var raw = self.atob(base64);
        var rawLength = raw.length;
        var array = new Uint8Array(new ArrayBuffer(rawLength));

        for (var i = 0; i < rawLength; i++) {
            array[i] = raw.charCodeAt(i);
        }

        return array;
    }
}

module.exports = function (self) {
    self.registerWorkerSource('mbtiles', MBTilesWorkerSource)
};
