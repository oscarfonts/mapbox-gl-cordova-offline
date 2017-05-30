const VectorTileSource = require('mapbox-gl/src/source/vector_tile_source');
const webworkify = require('webworkify');

class CustomVectorTileSource extends VectorTileSource {
    constructor(id, options, dispatcher, eventedParent) {
        super(id, options, dispatcher, eventedParent);
        this.type="custom_vector";
    }
}

CustomVectorTileSource.workerSourceURL = URL.createObjectURL(webworkify(require('./custom_vector_tile_worker.js'), {bare: true}));

module.exports = CustomVectorTileSource;
