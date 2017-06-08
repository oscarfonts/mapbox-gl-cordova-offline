const VectorTileSource = require('mapbox-gl/src/source/vector_tile_source');
const normalizeURL = require('mapbox-gl/src/util/mapbox').normalizeTileURL;
const webworkify = require('webworkify');
const pako = require('pako/lib/inflate');

class MBTilesSource extends VectorTileSource {
    constructor(id, options, dispatcher, eventedParent) {
        super(id, options, dispatcher, eventedParent);
        this.type = "mbtiles";
    }

    readTile(z, x, y, callback) {
        if ('sqlitePlugin' in self) {
            copyDatabaseFile('barcelona.mbtiles').then(function () {
                var db = sqlitePlugin.openDatabase({name: 'barcelona.mbtiles', location: 'default'});
                db.transaction(function (txn) {
                    txn.executeSql('SELECT BASE64(tile_data) AS data FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?', [z, x, Math.pow(2,z)-y-1], function (tx, res) {
                        var pbf_gz = convertBase64ToUint8Array(res.rows.item(0).data);
                        var pbf = pako.inflate(pbf_gz);
                        callback(null, {
                            data: pbf,
                            cacheControl: null,
                            expires: null
                        }); // Tile contents read
                    });
                }, function (error) {
                    callback(error); // Error executing SQL
                });
            }).catch(function (err) {
                callback(err); // Error copying database file
            });
        } else {
            callback("cordova-sqlite-ext plugin not installed"); // Required plugin not installed
        }

        function convertBase64ToUint8Array(base64) {
            var raw = self.atob(base64);
            var rawLength = raw.length;
            var array = new Uint8Array(new ArrayBuffer(rawLength));

            for (var i = 0; i < rawLength; i++) {
                array[i] = raw.charCodeAt(i);
            }

            return array;
        }

        function copyDatabaseFile(dbName) {
            var sourceFileName = cordova.file.applicationDirectory + 'www/data/' + dbName;

            return Promise.all([
                new Promise(function (resolve, reject) {
                    resolveLocalFileSystemURL(sourceFileName, resolve, reject);
                }),
                new Promise(function (resolve, reject) {
                    // If android
                    resolveLocalFileSystemURL(cordova.file.applicationStorageDirectory, function (dir) {
                        dir.getDirectory('databases', {create: true}, function (subdir) {
                            resolve(subdir);
                        });
                    }, reject);
                    // TODO else if ios
                    //resolveLocalFileSystemURL(cordova.file.documentsDirectory, resolve, reject);
                })
            ]).then(function (files) {
                var sourceFile = files[0];
                var targetDir = files[1];
                return new Promise(function (resolve, reject) {
                    targetDir.getFile(dbName, {}, resolve, reject);
                }).then(function () {
                    console.log("file already copied");
                }).catch(function () {
                    console.log("file doesn't exist, copying it");
                    return new Promise(function (resolve, reject) {
                        sourceFile.copyTo(targetDir, dbName, resolve, reject);
                    }).then(function () {
                        console.log("database file copied");
                    });
                });
            });
        }

    }


    loadTile(tile, callback) {
        const overscaling = tile.coord.z > this.maxzoom ? Math.pow(2, tile.coord.z - this.maxzoom) : 1;
        this.readTile(tile.coord.z, tile.coord.x, tile.coord.y, dispatch.bind(this));

        function dispatch(err, vectorTile) {
            if (err) {
                done(err);
            }

            const params = {
                url: normalizeURL(tile.coord.url(this.tiles, this.maxzoom, this.scheme), this.url),
                uid: tile.uid,
                coord: tile.coord,
                zoom: tile.coord.z,
                tileSize: this.tileSize * overscaling,
                type: this.type,
                source: this.id,
                overscaling: overscaling,
                angle: this.map.transform.angle,
                pitch: this.map.transform.pitch,
                showCollisionBoxes: this.map.showCollisionBoxes,
                vectorTile: vectorTile
            };

            if (!tile.workerID || tile.state === 'expired') {
                tile.workerID = this.dispatcher.send('loadTile', params, done.bind(this));
            } else if (tile.state === 'loading') {
                // schedule tile reloading after it has been loaded
                tile.reloadCallback = callback;
            } else {
                this.dispatcher.send('reloadTile', params, done.bind(this), tile.workerID);
            }

            function done(err, data) {
                if (tile.aborted)
                    return;

                if (err) {
                    return callback(err);
                }

                if (this.map._refreshExpiredTiles) tile.setExpiryData(data);
                tile.loadVectorData(data, this.map.painter);

                if (tile.redoWhenDone) {
                    tile.redoWhenDone = false;
                    tile.redoPlacement(this);
                }

                callback(null);

                if (tile.reloadCallback) {
                    this.loadTile(tile, tile.reloadCallback);
                    tile.reloadCallback = null;
                }
            }
        }

    }

}

MBTilesSource.workerSourceURL = URL.createObjectURL(webworkify(require('./mbtiles_worker.js'), {bare: true}));

module.exports = MBTilesSource;
