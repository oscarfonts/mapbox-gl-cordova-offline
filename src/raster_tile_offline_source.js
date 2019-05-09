import RasterTileSource from 'mapbox-gl/src/source/raster_tile_source';
import { extend, pick } from 'mapbox-gl/src/util/util'
import Texture from 'mapbox-gl/src/render/texture'
import Database from './database'

class RasterTileSourceOffline extends RasterTileSource {

    constructor(id, options, dispatcher, eventedParent) {
        super(id, options, dispatcher, eventedParent);
        this.id = id;
        this.dispatcher = dispatcher;
        this.setEventedParent(eventedParent);

        this.type = 'rasteroffline';
        this.minzoom = 0;
        this.maxzoom = 22;
        this.roundZoom = true;
        this.scheme = 'xyz';
        this.tileSize = 512;
        this.imageFormat = 'png';
        this._loaded = false;
        this._options = extend({}, options);
        extend(this, pick(options, ['scheme', 'tileSize', 'imageFormat']));

        this._transparentPngUrl = 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAC0lEQVQYV2NgAAIAAAUAAarVyFEAAAAASUVORK5CYII=';

        this.db = this.openDatabase(options.path)
    }

    openDatabase(dbLocation) {
        return Database.openDatabase(dbLocation)
    }

    copyDatabaseFile(dbLocation, dbName, targetDir) {
        return Database.copyDatabaseFile(dbLocation, dbName, targetDir)
    }

    loadTile(tile, callback) {

        tile.request = this._getImage(tile.tileID.canonical, done.bind(this));

        function done(err, img) {
            delete tile.request;

            if (tile.aborted) {
                tile.state = 'unloaded';
                callback(null);
            } else if (err) {
                tile.state = 'errored';
                callback(err);
            } else if (img) {
                if (this.map._refreshExpiredTiles) tile.setExpiryData(img);
                delete img.cacheControl;
                delete img.expires;

                const context = this.map.painter.context;
                const gl = context.gl;
                tile.texture = this.map.painter.getTileTexture(img.width);
                if (tile.texture) {
                    tile.texture.update(img, { useMipmap: true });
                } else {
                    tile.texture = new Texture(context, img, gl.RGBA, { useMipmap: true });
                    tile.texture.bind(gl.LINEAR, gl.CLAMP_TO_EDGE, gl.LINEAR_MIPMAP_NEAREST);

                    if (context.extTextureFilterAnisotropic) {
                        gl.texParameterf(gl.TEXTURE_2D, context.extTextureFilterAnisotropic.TEXTURE_MAX_ANISOTROPY_EXT,
                            context.extTextureFilterAnisotropicMax);
                    }
                }

                tile.state = 'loaded';

                callback(null);
            }
        }
    }

    _getBlob(coord, callback){

        const coordY = Math.pow(2, coord.z) -1 - coord.y;

        const query = 'SELECT BASE64(tile_data) AS base64_tile_data FROM tiles WHERE zoom_level = ? AND tile_column = ? AND tile_row = ?';
        const params = [coord.z, coord.x, coordY];

        const base64Prefix = 'data:image/' + this.imageFormat + ';base64,';

        this.db.then((db) => {
            db.transaction((txn) => {
                txn.executeSql(query, params, (tx, res) => {
                    if (res.rows.length) {
                        callback(undefined,
                            {
                                data: base64Prefix + res.rows.item(0).base64_tile_data,
                                cacheControl: null,
                                expires: null
                            });

                    } else {
                        console.error('tile ' + params.join(',') + ' not found');
                        callback(undefined,
                            {
                                data: this._transparentPngUrl,
                                cacheControl: null,
                                expires: null
                            });
                    }
                });
            }, (error) => {
                callback(error); // Error executing SQL
            });
        }).catch((err) => {
            callback(err);
        });
    }


    _getImage(coord, callback) {

        return this._getBlob(coord, (err, imgData) => {
            if (err) return callback(err);

            const img = new window.Image();
            const URL = window.URL || window.webkitURL;
            img.onload = () => {
                callback(null, img);
                URL.revokeObjectURL(img.src);
            };
            img.cacheControl = imgData.cacheControl;
            img.expires = imgData.expires;
            img.src = imgData.data;
        });

    }
}

export default RasterTileSourceOffline;
