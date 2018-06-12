import mapboxgl from 'mapbox-gl/src/index'

import offlineMap from './offline_map'
import language from './openmaptiles-language.js'

mapboxgl.OfflineMap = offlineMap;
language(mapboxgl);

export default mapboxgl
