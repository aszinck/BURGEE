var POI = /* color: #d63000 */ee.Geometry.Point([120.48462694850305, -67.07988982442774]);

var region = 'ZoneWilkes1';


var startEpoch = '2010-07-01';
var endEpoch = '2011-07-01';


var REMAstart = '2010-07-01';
var REMAend = '2018-07-01';

var saveToAsset = true;

print('This script performs the Lagrangian displacement \nin the '+region+' region.');




/*
Grounding line
*/
var AOI = ee.FeatureCollection('users/aszinck/GlaciologicalData/regionsCollection').filterMetadata('region','equals',region).geometry();

Map.setOptions('SATELLITE');
Map.centerObject(POI,6.5);


var IS = ee.FeatureCollection('users/azcryosphere/GlaciologicalData/FilteredIS').filterBounds(POI);
Map.addLayer(IS,{color:'red'},'IS')

// IS = ee.Feature(IS.first()).intersection(ee.Feature(shelf4))
// IS = IS.geometry()
// Map.addLayer(IS,{color:'yellow'},'IS used')

// var newG = ee.Feature(shelf5).difference(ee.Feature(shelf4))
// IS = ee.Feature(IS.first()).intersection(newG)
// IS = IS.geometry()
// Map.addLayer(IS,{color:'yellow'},'IS used')




/*
Constant velocity field used for initial lagrangian displacement
*/

var itsLiveRaw = ee.Image('users/earthmapps/Quantarctica/ANT_G0240_0000_vxy');
var b1 = itsLiveRaw.select(['b1']);
var b2 = itsLiveRaw.select(['b2']);
b2 = b2.multiply(ee.Image.constant(-1));
var itsLive = b1.addBands(b2);

var velMag = b1.multiply(b1).add(b2.multiply(b2)).sqrt();
Map.addLayer(velMag,{min:0, max:800, palette:'purple,blue,cyan,lime,yellow,orange,red,brown,maroon'},'Velocity',false);


/*
Visualization parameters
*/

var palettes = require('users/gena/packages:palettes');
var elevVis = {min:0, max:120, palette:'450256,3B1C8C,21908D,5AC865,F9E721,white'};
var vmin = 10;
var trendVis = {min: -vmin, max: vmin, palette:palettes.cmocean.Balance[7]};
var meltVis = {min:0, max: 35, palette:'black,darkmagenta,red,orange,yellow,white'};
var meltVisGour = {min:-100, max: 100, palette:palettes.cmocean.Balance[7].reverse()};

var tools = require('users/fitoprincipe/geetools:tools'); 


/*
Divergence field
*/

var getDivergences = function(image){
  var imageMask = image.select('elevation').mask();
  var xgradMin = ee.Image('users/aszinck/Divergences/'+region+'Regularized_xgradMin');//.resample('bicubic');
  var ygradPlus = ee.Image('users/aszinck/Divergences/'+region+'Regularized_ygradPlus');//.resample('bicubic');
  
  var divXMYP = xgradMin.add(ygradPlus).updateMask(imageMask).rename('divergence');

  return image.addBands(divXMYP);
};



/*
REMA collection to use
*/

var REMA = ee.ImageCollection("users/azcryosphere/"+region+"/REMAwFAC").filterDate(REMAstart,REMAend).filterBounds(IS);

print('REMA',REMA);


Map.addLayer(REMA.filterDate(startEpoch,endEpoch).select('elevation'),elevVis,'Before')

/*
The huge function which is used to calculate both the Lagrangian elevation change and the basal melt rate

Input to function:
startYear: the first year to consider (2010 --> 2010/11) for each dt 
dt: the trend period
*/

var startYear = 2010;
var dt = 6; // years, and I would like it to be able to be any number between 1 and max dt

// var hugeFunction = function(startYear,dt,maxdhdt){

var endYear = startYear + dt;

// This would mean the refDate would have to be this:
var refYear = endYear + 1;
var refDate = ee.Date.fromYMD(refYear,1,1); // if endYear = 2017 then refDate = 2018-01-01


/*
Velocity displacement to refDate of all strips
*/

var velocityDisplacement = function(image){
  var idate = ee.Date(ee.Image(image).get('system:time_start')) ;
  var ddate = idate.difference(refDate,'year');
  var dis = itsLive.multiply(ddate);
  var dimage = image.displace(dis);
  return dimage.set({'system:time_start':idate.millis()});
};


/*
Create system time start band
*/

var addTimeBand = function(image) {
  var timeBand = image.metadata('system:time_start').divide(1000 * 60 * 60 * 24 * 365).subtract(40).rename('timeband');
  return image.addBands(timeBand.updateMask(image.select(['elevation']).mask()));
};



/*
Make yearly qualityMosaics of all bands
*/


var refImage = REMA.select('elevation').filterDate('2015-07-01','2018-07-01');
refImage = refImage.map(velocityDisplacement);


// var refImage = REMA.filterDate(endYear+'-08-01',(endYear+1)+'-06-01');


// refImage = refImage.qualityMosaic('elevation');
refImage = refImage.median();

Map.addLayer(refImage,elevVis,'refImage');

/*
secondDisplacement function using refImage
refImage is the qualityMosaic of the endYear
*/

var secondDisplacement = function(refImage){
  return(
    function(image){
      var elevImage = image.select(['elevation']);
      elevImage = elevImage.setDefaultProjection('EPSG:3031').clip(AOI);
      refImage = refImage.select('elevation').setDefaultProjection('EPSG:3031').clip(AOI);
      var imageDisplace = elevImage.displacement({
                                referenceImage:refImage,
                                maxOffset:300,
                                patchWidth:100,
                                stiffness:3});
      image = image.displace(imageDisplace);
      return image.clip(IS);
})};



var addK3band = function(image) {
  var k3old = image.select(['divergence']).multiply(image.select(['elevation']).subtract(image.select(['fac']))).rename('k3');
  return image.addBands(k3old);
};



var REMAold = REMA.filterDate(startEpoch,endEpoch);

REMAold = REMAold.map(getDivergences);
REMAold = REMAold.map(velocityDisplacement);
Map.addLayer(REMAold.select('elevation'),elevVis,'Vel dis')
REMAold = REMAold.map(addTimeBand);
REMAold = REMAold.map(addK3band);
REMAold = REMAold.map(secondDisplacement(refImage));


Map.addLayer(REMAold.select('elevation'),elevVis,'After');

print(REMAold)

////////////////////////////////////////////////////////////////////
//////////////////////////// IS MASKING ////////////////////////////
////////////////////////////////////////////////////////////////////


function clipCollection(image){
  var mask = ee.Image.constant(1).clip(IS).mask();
  image = image.updateMask(mask);
  var areaImage = ee.Image.constant(1).updateMask(image.select('elevation').mask());
  var area = areaImage.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: IS.geometry(),
    scale: 500
    });

  return image.set('Area',ee.Number(area.get('constant')));
}



REMAold = REMAold.map(clipCollection);

REMAold = REMAold.filter(ee.Filter.metadata('Area','greater_than',1));

print('After IS masking',REMAold);


var downloadAll = function(collection, assetFolder, options) {
  
  var defaults = {
      name: null,
      scale: 1000,
      maxPixels: 1e13,
      region: null
    }
    
  var opt = tools.get_options(defaults, options)
  var n = collection.size().getInfo();
    
  var colList = collection.toList(n);
  
  for (var i = 0; i < n; i++) {
    var img = ee.Image(colList.get(i));
    var id = img.id().getInfo() || 'image_'+i.toString();
    // var region = opt.region || img.geometry(1).bounds(1).getInfo()["coordinates"];
    var region = opt.region || img.geometry(1);
    var assetId = assetFolder+'/'+id
    
    Export.image.toAsset({
      image: img,
      description: id,
      assetId: assetId,
      region: region,
      scale: opt.scale,
      maxPixels: opt.maxPixels})
  }
}



if (saveToAsset === true) { 
  downloadAll(REMAold,
    region+'/DisplacedREMA', 
    {name:'DisplacedREMA',
    scale:50,
    region: IS.geometry(1).bounds(1).getInfo()
})}


