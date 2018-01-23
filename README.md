# mapbox-gl + cordova + mbtiles

A mapbox-gl build capable of reading offline maps in mbtiles format in cordova.
Tested on Android, should work on iOS as well.


## Run sample application

```
cordova platform add android
cordova run android
```

```
cordova platform add ios
cordova run ios
```

Will use `www/data/barcelona.mbtiles` as sample data source, and `www/data/klokantech-basic.json` as sample style, both
coming from the OpenMapTiles project: https://openmaptiles.org/


## Develop with live reload

1. Get the IP address from your development computer (`ifconfig`).
2. Edit `www/index.html` and put your IP address in the script tag that loads the `mapbox-gl-cordova-mbtiles.js` resource:
   `<script src='http://xxx.xxx.xxx.xxx:8080/www/mapbox-gl-cordova-mbtiles.js'></script>`. For live reload to work,
   change also the IP_ADDRESS_AND_PORT var, and uncomment the code block at the end of the document.
3. Make sure to have your device attached (`adb devices`).
4. Run `npm start`.
5. For remote debugging, open `chrome://inspect/#devices` in your chrome browser.

Every time the contents in `src/` are changed, the file `www/mapbox-gl-cordova-mbtiles.js` will be rebuilt, and the
web view will be reloaded.


## Use in your application

Use the bundled library from `www/mapbox-gl-cordova-mbtiles.js` which is based in mapbox-gl-js v.0.37.0, or install it
as npm dependency (`npm install oscarfonts/mapbox-gl-cordova-mbtiles`).

Add the following cordova plugins via "cordova plugin add" command:

    * "cordova-plugin-device"
    * "cordova-plugin-file"
    * "cordova-sqlite-ext"


Mapbox-gl is not able to load a style with an mbtiles source directy (at instantiation), so use a base style to create
the map, then register the mbtiles source type, and finally load the source.
  
```javascript
var map = new mapboxgl.Map({
    container: 'map',
    style: 'data/empty.json'
});

map.on('load', function () {
    map.addSourceType('mbtiles', mapboxgl.MBTilesSource, function () {
        map.addSource('openmaptiles', {
            type: 'mbtiles',
            name: 'barcelona.mbtiles' /* The mbtiles file should be located in www/data/ */ 
        });
    });
});
```

See `www/index.html` in this repo for a complete example.
