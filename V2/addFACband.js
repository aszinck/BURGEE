
var fac = ee.ImageCollection("users/azcryosphere/IMAU-FDM/IMAU-FDM_fac_ERA5_extrapolated");


var region = 'ZoneWilkes1';

var saveFacToAsset = true;




Map.setOptions('SATELLITE');


fac = fac.filterDate('2010-07-01','2022-01-01');


var tools = require('users/fitoprincipe/geetools:tools'); 


var downloadAll = function(collection, assetFolder, options) {
  
  var defaults = {
      name: null,
      scale: 1000,
      maxPixels: 1e13,
      region: null
    };
    
  var opt = tools.get_options(defaults, options);
  var n = collection.size().getInfo();
    
  var colList = collection.toList(n);
  
  for (var i = 0; i < n; i++) {
    var img = ee.Image(colList.get(i));
    var id = img.id().getInfo() || 'image_'+i.toString();
    // var region = opt.region || img.geometry(1).bounds(1).getInfo()["coordinates"];
    var region = opt.region || img.geometry(1);
    var assetId = assetFolder+'/'+id;
    
    Export.image.toAsset({
      image: img,
      description: id,
      assetId: assetId,
      region: region,
      scale: opt.scale,
      maxPixels: opt.maxPixels});
  }
};




var facFunc = function(mask){
  return(
  function(image){
  // return image.updateMask(mask).resample('bicubic').reproject({crs: 'EPSG:3031',scale: 50});
  return image.resample('bicubic').reproject({crs: 'EPSG:3031',scale: 50});
})};

var REMA = ee.ImageCollection("users/azcryosphere/"+region+"/FullyCorrectedREMA").select(['elevation']);
print('REMA',REMA.size());


var bicubicFunc = function(image){
  var mask = image.mask();
  var facIm = fac.filterDate(image.date().advance(-15,'day'),image.date().advance(15,'day'));
  facIm = facIm.map(facFunc(mask));
  facIm = facIm.mean().updateMask(mask).rename('fac');
  image = image.addBands(facIm);
return image};

REMA = REMA.map(bicubicFunc);
print(REMA)

// Map.addLayer(REMA.first().select('fac'),{min:17,max:18});

// Map.addLayer(REMA.first().geometry())


if (saveFacToAsset === true) { 
downloadAll(REMA,
  region + '/REMAwFAC', 
  {name:'REMAwFAC',
  scale:50
})}



