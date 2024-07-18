var POI = /* color: #d63000 */ee.Geometry.Point([144.71840415665153, -67.51034505589391]);

/*
Grounding line
*/
var region = 'ZoneWilkes1';

var coverageCriteria = 25; // %

var countNum = 3; // years

var iceShelfName = 'Totten';

// var latCri = 0.9;
// var lonCri = 0.9;
var latCri = 0;
var lonCri = 0;

var saveToAsset = true;

var AOI = ee.FeatureCollection('users/aszinck/GlaciologicalData/regionsCollection').filterMetadata('region','equals',region).geometry();

var IS = ee.FeatureCollection('users/azcryosphere/GlaciologicalData/FilteredIS').filterMetadata('IS','equals',1).select('IS').filterBounds(POI);


Map.setOptions('SATELLITE');
Map.centerObject(POI,7.5);

print('/BMB_covCri_' + coverageCriteria.toString() + '_countNum_' + countNum.toString() + '/' + iceShelfName)

/*
Visualization parameters
*/

var palettes = require('users/gena/packages:palettes');
var elevVis = {min:0, max:30, palette:'450256,3B1C8C,21908D,5AC865,F9E721,white'};
var vmin = 10;
var trendVis = {min: -vmin, max: vmin, palette:palettes.cmocean.Balance[7]};
var meltVis = {min:0, max: 150, palette:'black,darkmagenta,red,orange,yellow,white'};
var meltVisGour = {min:-100, max: 100, palette:palettes.cmocean.Balance[7].reverse()};


/*
REMA collection to use
*/

var startYear = 2010;
var endYear = 2017;



if (region == 'ZoneVictoria1' || region == 'ZoneVictoria2' || region == 'ZoneWilkes1') {
  var user = 'users/azcryosphere';
} else if (region == 'ZoneAmundsen1' || region == 'ZoneAmundsen2' || region == 'ZoneAmundsen3') {
  var user = 'users/aszinck';
} else if (region == 'ZoneBellinghausen1' || region == 'ZoneBellinghausen2' || region == 'ZoneLarsen') {
  var user = 'users/azburgee';
} else if (region == 'ZoneWilkes2') {
  var user = 'projects/ee-zinckhoogendoorn/assets';
} else {
  print('Region not recognized!');
}






var facFunc = function(mask){
  return(
  function(image){
  return image.resample('bicubic').reproject({crs: 'EPSG:3031',scale: 50});
})};




var bicubicFunc = function(image){
  var mask = image.select('elevation').mask();
  var facIm = ee.ImageCollection('users/azcryosphere/IMAU-FDM/IMAU-FDM_fac_ERA5_extrapolated').filterDate(image.date().advance(-10,'day'),image.date().advance(10,'day')).select('b1');
  facIm = facIm.map(facFunc(mask));
  facIm = facIm.mean().updateMask(mask).rename('fac');
  image = image.addBands(facIm);
return image};

var addK3band = function(image) {
  var k3old = image.select(['divergence']).multiply(image.select(['elevation']).subtract(image.select(['fac']))).rename('k3');
  return image.addBands(k3old);
};

// REMA = REMA.map(bicubicFunc);
// print(REMA)
var REMAraw = ee.ImageCollection(user+"/"+region+"/DisplacedREMA").filterBounds(IS);
print('REMAraw',REMAraw);
var REMA = ee.ImageCollection([]);

for (var i= startYear; i < endYear; i = i +1) {
  var REMAic = ee.ImageCollection(user+"/"+region+"/DisplacedREMA").filterDate(i+'-07-01',(i+1)+'-07-01').filterBounds(IS);
  REMAic = REMAic.select('elevation','divergence','timeband');
  REMAic = REMAic.filterMetadata('latFrac','greater_than',latCri).filterMetadata('lonFrac','greater_than',lonCri);
  REMAic = REMAic.map(bicubicFunc);
  REMAic = REMAic.map(addK3band);
  REMAic = REMAic.median();
  REMAic = REMAic.set('system:time_start',ee.Date((i+1)+'-01-01').millis());
  REMAic = REMAic.set('Year',i+'/'+(i+1));
  REMA = REMA.merge(ee.ImageCollection.fromImages([REMAic]));
}


var listOfImages = REMA.toList(REMA.size());
// print("old list", listOfImages);

var newList = listOfImages.map(function comprobeBandsNumber(ele){
  var new_list = ee.List([]); 
  var count = ee.Image(ele).bandNames().size();
  var comp = ee.Algorithms.If(count.eq(5), ele, 0);
  new_list = new_list.add(comp);
  return new_list;
}).flatten();

//removing zeroes in new list
newList = newList.removeAll([0]);

// print("new list", newList);

//creating new collection
var REMA = ee.ImageCollection(newList);





var oceanmask = ee.Image("users/aszinck/GlaciologicalData/BMoceanmask");
var REMAMosaic = ee.Image('users/aszinck/REMAmosaics/'+region+'REMAdemMosaic').updateMask(oceanmask.mask()).clip(IS);
// Map.addLayer(REMAMosaic,elevVis,'REMAMosaic')



var yearlyCoverage = function(image){

  var firstREMA = image.select('elevation');

  var maskREMA = firstREMA.mask();
  var maskMosaic = REMAMosaic.mask();
  
  var reg = IS.geometry();
  var scale = 1000;
  
  
  var fullArea = ee.Image.constant(1).updateMask(maskMosaic).reduceRegion({
    reducer:ee.Reducer.sum(), 
    geometry:reg,
    scale:scale,
    crs:'EPSG:3031',
    tileScale:16})
    .get('constant');

  
  var areaREMA = ee.Image.constant(1).updateMask(maskREMA).reduceRegion({
    reducer:ee.Reducer.sum(), 
    geometry:reg,
    scale:scale,
    crs:'EPSG:3031',
    tileScale:16})
    .get('constant');

  
  var Perc = ee.Number(areaREMA).divide(ee.Number(fullArea)).multiply(ee.Number(100))

  return image.set('OverlapPercent',Perc)
}

REMA = REMA.map(yearlyCoverage);




// var REMA = ee.ImageCollection("users/azcryosphere/"+region+"/DisplacedREMA").filterBounds(IS);

print('REMA',REMA);

REMA = REMA.filterMetadata('OverlapPercent','not_less_than',coverageCriteria);

print('REMA coverage filtered',REMA);


//////////// Same as Eulerian filtering ///////////

var countREMA = REMA.select('elevation').count().rename('count');
var countMask = countREMA.mask().updateMask(countREMA.select('count').gte(countNum));
Map.addLayer(countMask,{min:0,max:10},'countMask');

REMA = REMA.map(function(image){return image.updateMask(countMask)});

///////////////////////////////////////////////////



var lagrangian = REMA.select(['timeband','elevation']).reduce(ee.Reducer.linearFit());
lagrangian = lagrangian.select(['scale']).rename('Lagrangian');


var timdiff = (REMA.select('timeband').max()).subtract(REMA.select('timeband').min());
// Map.addLayer(timdiff,{},'timdiff');

/////////////// SMB ////////////////

var startYear = ee.Date(REMA.first().get('system:time_start')).get('year');
print('First DEM year (XXXX-01-01):', startYear);
var listOfImages = REMA.toList(REMA.size());
var endYear = ee.Date(REMA.limit(1, 'system:time_start', false).first().get('system:time_start')).get('year');
print('Last DEM year (XXXX-01-01):', endYear);

var emptySMB = ee.ImageCollection([]);

var smbFunc = function(mask){
  return(
  function(image){
  return image.updateMask(mask).resample('bicubic').reproject({crs: 'EPSG:3031',scale: 1000});
})};


for (var i= startYear.getInfo(); i < endYear.getInfo(); i = i +1) {
  // print('i',i)
  var oceanmask = ee.Image("users/aszinck/GlaciologicalData/BMoceanmask");
  var REMAMosaic = ee.Image('users/aszinck/REMAmosaics/'+region+'REMAdemMosaic').updateMask(oceanmask.mask());
  var REMAmask = REMAMosaic.mask();
  var smbic = ee.ImageCollection('users/azcryosphere/RACMO_IMAU/RACMO_ERA5_extrapolated').filterDate(i+'-07-01',(i+1)+'-07-01');
  smbic = smbic.map(smbFunc(REMAmask));
  smbic = smbic.mean().updateMask(REMAmask);
  emptySMB = emptySMB.merge(ee.ImageCollection.fromImages([smbic]));
}

var smb = emptySMB.mean();

var water2ice = ee.Number(0.917);
var mm2m = ee.Number(0.001); 
smb = smb.divide(water2ice).multiply(mm2m).multiply(ee.Number(12));

/////////////// FAC CHANGE ////////////////

var emptyFAC = ee.ImageCollection([]);


for (var i= startYear.getInfo(); i < endYear.getInfo(); i = i +1) {
  var oceanmask = ee.Image("users/aszinck/GlaciologicalData/BMoceanmask");
  var REMAMosaic = ee.Image('users/aszinck/REMAmosaics/'+region+'REMAdemMosaic').updateMask(oceanmask.mask());
  var REMAmask = REMAMosaic.mask();
  var FACic = ee.ImageCollection('users/azcryosphere/IMAU-FDM/IMAU-FDM_fac_ERA5_extrapolated').filterDate(i+'-07-01',(i+1)+'-07-01');
  FACic = FACic.map(smbFunc(REMAmask));
  FACic = FACic.mean().updateMask(REMAmask).rename('fac');
  FACic = FACic.set('timeband',ee.Date((i+1)+'-01-01').millis())
  FACic = FACic.addBands(FACic.metadata('timeband').divide(1000 * 60 * 60 * 24 * 365).subtract(40)); 
  emptyFAC = emptyFAC.merge(ee.ImageCollection.fromImages([FACic]));
}

var facChange = emptyFAC.select(['timeband','fac']).reduce(ee.Reducer.linearFit()).select(['scale']).rename('facChange');
// print('facChange',facChange)


// Map.addLayer(facChange.clip(IS),{min:-2,max:2,palette:'green,white,blue'},'fac change',false);



/////////////// BMB ////////////////

var rhoWater = ee.Image.constant(1026).toFloat().select('constant');
rhoWater = rhoWater.updateMask(lagrangian.mask());

// Constant ice density of 917 kg m^(-3) 
var rhoIce = ee.Image.constant(917).toFloat().select('constant');
rhoIce = rhoIce.updateMask(lagrangian.mask());

// Map.addLayer(REMA.select(['k3']),{palette:palettes.cmocean.Balance[7]},'k3 i')
// Map.addLayer(REMA.select(['elevation']),elevVis,'REMA')
var k3 = REMA.select(['k3']).mean().divide(120); // mean(div(v)*(h-d))
// Map.addLayer(k3.clip(IS),{min:-10,max:10,palette:palettes.cmocean.Balance[7]},'k3',false)
// (h-h_f)*div(V)
var k0 = lagrangian.subtract(facChange); // Dh/Dt - Dhf/Dt
// var k0 = lagrangian;
var k1 = k0.add(k3); // Dh/Dt - Dhf/Dt + div(v)*(h-d)
var k2 = rhoWater.divide(rhoWater.subtract(rhoIce)); // rho_w/(rho_w-rho_i)
var kk = k1.multiply(k2); // (Dh/Dt + div(v)*(h-d)) * (rho_w/(rho_w-rho_i))
var basalMeltRate = smb.subtract(kk.toFloat()).rename('BMB'); // a - (Dh/Dt + div(v)*(h-d)) * (rho_w/(rho_w-rho_i))

lagrangian = lagrangian.addBands(basalMeltRate).addBands(facChange);

// lagrangian = lagrangian.set('dt',dt,'startYear',startYear + '/' + (startYear+1), 'endYear', endYear + '/' + (endYear+1));






// print(lagrangian)

var BMB = lagrangian.select('BMB').clip(IS);
print('BMB',BMB)
var DHDt = lagrangian.select('Lagrangian').clip(IS);

if (saveToAsset === true){

  
  
  Export.image.toAsset({
    image: BMB,
    description: 'BMB_'+iceShelfName,
    assetId: 'BURGEEout' + '/BMB_clip_covCri_' + coverageCriteria.toString() + '_countNum_' + countNum.toString() + '/' + iceShelfName,
    scale: 50,
    region: IS.geometry(),
    maxPixels: 10000000000000
  });
  
  
  Export.image.toAsset({
    image: DHDt,
    description: 'Lagrangian_'+iceShelfName,
    assetId: 'BURGEEout' + '/Lagrangian_clip_covCri_' + coverageCriteria.toString() + '_countNum_' + countNum.toString() + '/' + iceShelfName,
    scale: 50,
    region: IS.geometry(),
    maxPixels: 10000000000000
  });
  

}



Map.addLayer(lagrangian.select('Lagrangian').clip(IS),trendVis,'Lagrangian');
Map.addLayer(lagrangian.select('BMB').clip(IS),meltVis,'BMB');


var sui  = require('users/earthmapps/public:ui.js');
var colorleg = ui.Panel({
    style: {
        width: '300px',
        position: 'bottom-right',
        padding: '8px 15px'
    }
});
// Add legend colorbars to panel
// colorleg.widgets().set(0, sui.colorbar(trendVis,'elevation trend [m/yr]'));
// colorleg.widgets().set(0, sui.colorbar({min:-2,max:2,palette:'green,white,blue'},'fac change [m/yr]'));
// colorleg.widgets().set(0, sui.colorbar(meltVis,'basal mass balance [m/yr]'));
// colorleg.widgets().set(0, sui.colorbar(meltVisGour,'basal mass balance [m/yr]'));
// Map.add(colorleg);





