# Offline vector maps in Cordova using Mapbox GL JS

A mapbox-gl-js build capable of reading local mbtiles in cordova.
Tested on Android, and (to a lesser extent) on iOS.


## Run example application

```
npm install
cordova platform add android
cordova run android
```

**The Android emulator browser won't display WebGL content, please run on a physical Android device**

```
npm install
cordova platform add ios
cordova run ios
```


Will use `www/data/2017-07-03_spain_barcelona.mbtiles` as sample data source, and `www/styles/osm-bright/style-offline.json`
as style, both coming from the OpenMapTiles project: https://openmaptiles.org/


## Integrate in your application

Use the bundled library from `www/mapbox-gl-cordova-offline.js` which is based in mapbox-gl-js v.0.53.1, or install it
as npm dependency (`npm install oscarfonts/mapbox-gl-cordova-offline`).

Add the following cordova plugins via "cordova plugin add" command:

    * "cordova-plugin-device"
    * "cordova-plugin-file"
    * "cordova-sqlite-ext"

Use the OfflineMap constructor. It returns a **promise** instead of a map, as the
offline map creation process is asynchronous:
  
```javascript
       new mapboxgl.OfflineMap({
            container: 'map',
            style: 'styles/osm-bright/style-offline.json'
       }).then(function(map) {
           map.addControl(new mapboxgl.NavigationControl());
       });
```

See `www/index.html` in this repo for a working example.


### Offline data sources (mbtiles)

In your style, you can specify offline tile sources specifying `mbtiles` as the source type,
and the location to the mbtiles file as a relative path:

```json
"sources": {
    "openmaptiles": {
        "type": "mbtiles",
        "path": "data/2017-07-03_spain_barcelona.mbtiles"
    }
}
```

Additional styles can be found in OpenMapTiles repos (see gh-pages branches): https://github.com/openmaptiles
Vector tiles for other regions can be downloaded here: https://openmaptiles.com/downloads/planet/

#### Raster tiles

You can add raster tiles in offline mode. In the style you have to add a new source with type `rasteroffline` with the 
relative path:

```json
"sources": {
    "hillshading": {
        "type": "rasteroffline",
        "path": "data/2016-11-28-hillshade-spain_barcelona.mbtiles"
    }
}
```

The style related to this source will be the same that a common Mapbox raster layer. You have to relate the layer with
the source using the source name, and the type of the layer have to be `raster` as a common raster Mapbox layer.

```json
{
      "id": "hillshading",
      "type": "raster",
      "source": "hillshading",
      "layout": {
        "visibility": "visible"
      },
      "paint": {
        "raster-opacity": {
          "base": 0.5,
          "stops": [
            [
              3,
              0
            ],
            [
              5,
              0.5
            ],
            [
              12,
              0.5
            ]
          ]
        },
        "raster-fade-duration": 300,
        "raster-contrast": 0
      }
    }
```

### Offline sprites (icon set) 

Copy the files `sprite.json`, `sprite.png`, `sprite@2x.json` and `sprite@2x.png` as local resources and
reference them as a relative path in your style:

```json
"sprite": "styles/osm-bright/sprite"
```


### Offline glyphs (fonts) 

Search "text-font" attributes in your style. Download the needed fonts from https://github.com/openmaptiles/fonts
(see gh-pages branch) and copy them locally. Set the relative path in the "glyphs" property of the
style:

```json
"glyphs": "fonts/{fontstack}/{range}.pbf"
```


## Enable live reload for development

1. Get the your development computer's IP address (`ifconfig`).
2. Edit `www/index.html` and put your IP address in the script tag that loads the `mapbox-gl-cordova-offline.js` resource:
   `<script src='http://xxx.xxx.xxx.xxx:8080/www/mapbox-gl-cordova-offline.js'></script>`. For live reload to work,
   change also the IP_ADDRESS_AND_PORT var, and uncomment the code block at the end of the document.
3. Run `npm start`.

Every time the contents in `src/` are changed, the file `www/mapbox-gl-cordova-offline.js` will be rebuilt, and the
web view will be reloaded.
