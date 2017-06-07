const vt = require('vector-tile');
const Protobuf = require('pbf');
const pako = require('pako/dist/pako_inflate');
const VectorTileWorkerSource = require('mapbox-gl/src/source/vector_tile_worker_source');

class MBTilesWorkerSource extends VectorTileWorkerSource {
    loadVectorData(params, callback) {
        const xhr = this.readTile(params.coord.z, params.coord.x, params.coord.y, done.bind(this));
        return function abort() {
            xhr.abort();
        };
        function done(err, response) {
            if (err) {
                return callback(err);
            }
            const vectorTile = new vt.VectorTile(new Protobuf(response.data));
            vectorTile.rawData = response.data;
            vectorTile.cacheControl = response.cacheControl;
            vectorTile.expires = response.expires;
            callback(err, vectorTile);
        }
    }

     readTile(z, x, y, callback) {
        if ('sqlitePlugin' in self) {
            this.copyDatabaseFile('barcelona.mbtiles').then(function () {
                var db = sqlitePlugin.openDatabase({name: 'barcelona.mbtiles', location: 'default'});
                db.transaction(function (txn) {
                    txn.executeSql('SELECT BASE64(tile_data) AS data FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?', [z, x, y], function (tx, res) {
                        var pbf_gz = this.convertBase64ToUint8Array(res.rows.item(0).data);
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

    copyDatabaseFile(dbName) {
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

module.exports = function (self) {
    self.registerWorkerSource('mbtiles', MBTilesWorkerSource)
};
