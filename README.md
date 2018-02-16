# mapbox-gl + cordova + mbtiles

A mapbox-gl build capable of reading offline maps in mbtiles format in cordova.
Tested on Android, and (to a lesser extent) on iOS.


## Run exsample application

```
cordova platform add android
cordova run android
```

```
cordova platform add ios
cordova run ios
```

Will use `www/data/2017-07-03_spain_barcelona.mbtiles` as sample data source, and `www/styles/osm-bright/` as sample
style, both coming from the OpenMapTiles project: https://openmaptiles.org/


## Integrate in your application

Use the bundled library from `www/mapbox-gl-cordova-mbtiles.js` which is based in mapbox-gl-js v.0.44.0, or install it
as npm dependency (`npm install oscarfonts/mapbox-gl-cordova-mbtiles`).

Add the following cordova plugins via "cordova plugin add" command:

    * "cordova-plugin-device"
    * "cordova-plugin-file"
    * "cordova-sqlite-ext"


Added a new Map instantiator called `mapboxgl.OfflineMap`. It returns a **promise** instead of a map directly, as the
instantiation is asynchronous. Access the `map` as a `then` callback. For instance:
  
```javascript
       new mapboxgl.OfflineMap({
            container: 'map',
            style: 'styles/osm-bright/style-offline.json'
       }).then((map) => {
           map.addControl(new mapboxgl.NavigationControl());
       });
```

See `www/index.html` in this repo for a working example.

In your style, you can specify the `mbtiles` source and the location to the file, relative to application's base URL:

```json
"sources": {
    "openmaptiles": {
        "type": "mbtiles",
        "path": "data/2017-07-03_spain_barcelona.mbtiles"
    }
}
```

See styles in `styles` directory for working examples.


### Where to download data and styles

* Vector tiles for other regions can be downloaded here: https://openmaptiles.com/downloads/planet/
* Fonts are available in https://github.com/openmaptiles/fonts (see gh-pages branch).
* Additional styles can be found in OpenMapTiles repos (see gh-pages branches): https://github.com/openmaptiles


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
