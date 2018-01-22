const VectorTileSource = require('mapbox-gl/src/source/vector_tile_source');
const webworkify = require('webworkify');

class MBTilesSource extends VectorTileSource {
    constructor(id, options, dispatcher, eventedParent) {
        super(id, options, dispatcher, eventedParent);
        this.type = "mbtiles";
        this.db = this.openDatabase(options.name);
    }

    openDatabase(name) {
        if ('sqlitePlugin' in self) {
            return this.copyDatabaseFile(name).then(function () {
                var params = {name: name};
                if(device.platform === 'iOS') {
                    params.iosDatabaseLocation = 'Documents';
                } else {
                    params.location = 'default';
                }
                return sqlitePlugin.openDatabase(params);
            });
        } else {
            return Promise.reject(new Error("cordova-sqlite-ext plugin not available. " +
                "Please install the plugin and make sure this code is run after onDeviceReady event"));
        }
    }

    copyDatabaseFile(dbName) {
        const sourceFileName = cordova.file.applicationDirectory + 'www/data/' + dbName;

        if(!('device' in self)) {
            return Promise.reject(new Error("cordova-plugin-device not available. " +
                "Please install the plugin and make sure this code is run after onDeviceReady event"));
        }

        return Promise.all([
            new Promise(function (resolve, reject) {
                resolveLocalFileSystemURL(sourceFileName, resolve, reject);
            }),
            new Promise(function (resolve, reject) {
                if(device.platform === 'Android') {
                    resolveLocalFileSystemURL(cordova.file.applicationStorageDirectory, function (dir) {
                        dir.getDirectory('databases', {create: true}, function (subdir) {
                            resolve(subdir);
                        });
                    }, reject);
                } else if(device.platform === 'iOS') {
                    resolveLocalFileSystemURL(cordova.file.documentsDirectory, resolve, reject);
                } else {
                    reject("Platform not supported");
                }
            })
        ]).then(function (files) {
            const sourceFile = files[0];
            const targetDir = files[1];
            return new Promise(function (resolve, reject) {
                targetDir.getFile(dbName, {}, resolve, reject);
            }).catch(function () {
                console.log("Copying database to application storage directory");
                return new Promise(function (resolve, reject) {
                    sourceFile.copyTo(targetDir, dbName, resolve, reject);
                }).then(function () {
                    console.log("Database copied");
                });
            });
        });
    }

    readTile(z, x, y, callback) {
        const query = 'SELECT BASE64(tile_data) AS base64_tile_data FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?';
        const params = [z, x, y];
        this.db.then(function(db) {
            db.transaction(function (txn) {
                txn.executeSql(query, params, function (tx, res) {
                    const base64Data = res.rows.length ? res.rows.item(0).base64_tile_data : undefined;
                    callback(undefined, base64Data); // Tile contents read, callback success.
                });
            }, function (error) {
                callback(error); // Error executing SQL
            });
        }).catch(function(err) {
            callback(err);
        });
    }

    loadTile(tile, callback) {
        const overscaling = tile.coord.z > this.maxzoom ? Math.pow(2, tile.coord.z - this.maxzoom) : 1;

        const z = Math.min(tile.coord.z, this.maxzoom || tile.coord.z); // Don't try to get data over maxzoom
        const x = tile.coord.x;
        const y = Math.pow(2,z)-tile.coord.y-1; // Tiles on database are tms (inverted y axis)

        this.readTile(z, x, y, dispatch.bind(this));

        function dispatch(err, base64Data) {
            if (err) {
                return callback(err);
            }

            const params = {
                base64Data: base64Data,
                uid: tile.uid,
                coord: tile.coord,
                zoom: tile.coord.z,
                tileSize: this.tileSize * overscaling,
                type: this.type,
                source: this.id,
                overscaling: overscaling,
                angle: this.map.transform.angle,
                pitch: this.map.transform.pitch,
                showCollisionBoxes: this.map.showCollisionBoxes
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
