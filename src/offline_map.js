// @flow

'use strict';

import MBTilesSource from './mbtiles_source'
import Map from 'mapbox-gl/src/ui/map'
import {extend} from 'mapbox-gl/src/util/util'
import window from 'mapbox-gl/src/util/window'

const readJSON = (url) => new Promise((resolve, reject) => {
    const xhr = new window.XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.setRequestHeader('Accept', 'application/json');
    xhr.onerror = (e) => reject(e);
    xhr.onload = () => {
        const isFile = xhr.responseURL.indexOf('file://') === 0;
        if (((xhr.status >= 200 && xhr.status < 300) || isFile) && xhr.response) {
            try {
                resolve(JSON.parse(xhr.response));
            } catch (err) {
                reject(err);
            }
        } else {
            reject(new Error(xhr.statusText, xhr.status));
        }
    };
    xhr.send();
    return xhr;
});

const dereferenceStyle = (options) => {
    if (typeof options.style === 'string' || options.style instanceof String) {
        return readJSON(options.style).then((style) => extend({}, options, {style: style}));
    } else {
        return Promise.resolve(options);
    }
};

const absoluteSpriteUrl = (options) => {
    const style = options.style;
    const hasProtocol = /^.+:\/\//;
    const path = window.location.origin + window.location.pathname.split('/').slice(0, -1).join('/');

    if (('sprite' in style) && !style.sprite.match(hasProtocol)) {
        style.sprite = path + '/' +  style.sprite; // eslint-disable-line prefer-template
    }
    return options;
};

const createEmptyMap = (options) => new Promise((resolve) => {
    const emptyMapStyle = extend({}, options.style, {
        sources: {},
        layers: []
    });
    const emptyMapOptions = extend({}, options, {style: emptyMapStyle});
    const map = new Map(emptyMapOptions);
    map.once('load', () => map.addSourceType('mbtiles', MBTilesSource, () => resolve(map)));
});

const loadSources = (style) => (map) => {
    Object.keys(style.sources).map((sourceName) => map.addSource(sourceName, style.sources[sourceName]));
    return map;
};

const loadLayers = (style) => (map) => {
    style.layers.map((layer) => map.addLayer(layer));
    return map;
};

const OfflineMap = (options) =>
    dereferenceStyle(options).then(absoluteSpriteUrl).then((newOptions) =>
        createEmptyMap(newOptions)
            .then(loadSources(newOptions.style))
            .then(loadLayers(newOptions.style))
    );

export default OfflineMap
