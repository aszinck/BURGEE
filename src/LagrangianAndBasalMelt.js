/*
This script calculates the lagrangian elevation change and the basal melt rate

How to run the script:
1. Copy-paste to the code editor on Google Earth Engine
2. Optional: Change data path in line 58 in case you want to use your own co-registered DEMs
3. Change to own firn and surface mass balance data where marked in the code
3. Click run

Line 282-284 contains more changable parameters controlling the period in which one wants to calculate the trends over
*/


////////////////////////////////////////////////////////////////////
///////////////////////// LOAD DATA ETC. ///////////////////////////
////////////////////////////////////////////////////////////////////

// Initialize map

var AOI =
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-113.46433216556532, -74.12581799259875],
          [-114.02463489994032, -74.23065574862235],
          [-113.94223743900282, -74.3777749785877],
          [-113.86533314212782, -74.52795235449202],
          [-114.46958118900282, -74.70572718049199],
          [-114.40915638431532, -74.7780083376015],
          [-114.20041614994032, -74.77656598815501],
          [-114.48056751712782, -75.1203397642524],
          [-114.14548450931532, -75.20614288663648],
          [-111.31101185306532, -74.80106779060921],
          [-111.50327259525282, -74.58792014866033],
          [-111.17917591556532, -74.55138145708696],
          [-111.10227161869032, -74.45009560839453],
          [-111.81088978275282, -74.19329137944453]]])


Map.setOptions('SATELLITE');
Map.centerObject(AOI,7.5);

var IS = ee.FeatureCollection('users/earthmapps/IceShelves').filterMetadata('SURFACE','equals','ice shelf').filterBounds(AOI);

var itsLiveRaw = ee.Image('users/earthmapps/Quantarctica/ANT_G0240_0000_vxy');
var b1 = itsLiveRaw.select(['b1']);
var b2 = itsLiveRaw.select(['b2']);
b2 = b2.multiply(ee.Image.constant(-1));
var itsLive = b1.addBands(b2);

var palettes = require('users/gena/packages:palettes');
var vmin = 10;
var trendVis = {min: -vmin, max: vmin, palette:palettes.cmocean.Balance[7]};
var meltVis = {min:0, max: 25, palette:'black,darkmagenta,red,orange,yellow,white'};


// Use this imageCollection of fully correted REMA strips or use your own DEMs
var REMA = ee.ImageCollection("users/azcryosphere/DotsonIceShelf/REMADotsonFullCorrection_lat60_1month").filterBounds(IS).select(['elevation']);


print('REMA',REMA);

/*
Get divergences function
*/

var getDivergences = function(image){
  var imageMask = image.mask();
  var xgradMin = ee.Image('users/azcryosphere/DotsonRegularized_xgradMin');//.resample('bicubic');
  var ygradPlus = ee.Image('users/azcryosphere/DotsonRegularized_ygradPlus');//.resample('bicubic');

  var divXMYP = xgradMin.add(ygradPlus).updateMask(imageMask).rename('divergence');

  return image.addBands(divXMYP);
};


/*
Get firn air content function
OBS! Here you have to use your own firn data!
*/

var getFAC = function(image){
  var firnAirContent = ee.ImageCollection('users/azcryosphere/IMAU-FDM/IMAU-FDM_ERA5_extrapolated_bicubic');
  firnAirContent = firnAirContent.filterDate(ee.Date(image.get('system:time_start')).advance(-10,'days'),ee.Date(image.get('system:time_start')).advance(10,'days')).mean().rename('fac');
  return image.addBands(firnAirContent.updateMask(image.select(['elevation']).mask()));
};




var hugeFunction = function(startYear,dt){

  var endYear = startYear + dt;



  // This would mean the refDate would have to be this:
  var refYear = endYear + 1;
  var refDate = ee.Date.fromYMD(refYear,1,1);


  /*
  Velocity displacement to refDate of all strips, both REMAold and REMAnew
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
    var timeBand = image.metadata('system:time_start').divide(1000 * 60 * 60 * 24 * 365).subtract(40);
    return image.addBands(timeBand.updateMask(image.select(['elevation']).mask()));
  };


  /*
  Make yearly DEMs
  */

  var refImage = REMA.filterDate(endYear+'-07-01',(endYear+1)+'-07-01');
  // print(refImage)
  // print('refImage size',refImage.size())
  if(refImage.size().getInfo() == 0){
    print('empty refImage');
    var lagrangian = ee.Image.constant(0);
  } else {

  refImage = refImage.qualityMosaic('elevation');
  // print(refImage)

  /*
  Feature tracking displacement
  */
  var secondDisplacement = function(image){
    var elevImage = image.select(['elevation']);
    elevImage = elevImage.setDefaultProjection('EPSG:3031').clip(AOI);
    // print(elevImage)
    refImage = refImage.select('elevation').setDefaultProjection('EPSG:3031').clip(AOI);
    var imageDisplace = elevImage.displacement({
                              referenceImage:refImage,
                              maxOffset:300,
                              patchWidth:100,
                              stiffness:3});
    image = image.displace(imageDisplace);
    // print(image)
    return image;
  };

  // Elevation change outlier criteria (changeable)
  var maxdhdt = 15; // m/yr

  var outlier_mask = function(image){
    return image.updateMask(image.lte(maxdhdt).and(image.gte(-maxdhdt)))};

  var emptyIC = ee.ImageCollection([]);
  var emptyLagr = ee.ImageCollection([]);

  for (var i= startYear; i < endYear; i = i +1) {
    // print(i)
    var REMAold = REMA.filterDate(i+'-08-01',(i+1)+'-06-01');
    // print(REMAold.size())

    if(REMAold.size().getInfo() == 0){
      print('empty collection');
      var lagrangian = ee.Image.constant(0)
    } else {

    REMAold = REMAold.map(getDivergences);
    REMAold = REMAold.map(getFAC);
    REMAold = REMAold.map(velocityDisplacement);
    REMAold = REMAold.map(addTimeBand);
    REMAold = REMAold.qualityMosaic('elevation');
    var k3old = REMAold.select(['divergence']).multiply(REMAold.select(['elevation']).subtract(REMAold.select(['fac']))).rename('k3');
    REMAold = REMAold.addBands(k3old)
    REMAold = REMAold.set('MosaicYear',i+'/'+(i+1));
    var REMAoldie = secondDisplacement(REMAold);

    emptyIC = emptyIC.merge(REMAold);
    emptyLagr = emptyLagr.merge(ee.ImageCollection.fromImages([REMAoldie]));
    // print(REMAold)
  }}



  /*
  Calculate the Lagrangian elevation change
  */



  var lastImage = REMA.filterDate(endYear+'-07-01',(endYear+1)+'-07-01');
  lastImage = lastImage.map(getDivergences);
  lastImage = lastImage.map(getFAC);
  lastImage = lastImage.map(velocityDisplacement);
  lastImage = lastImage.map(addTimeBand);
  lastImage = lastImage.qualityMosaic('elevation');
  var k3last = lastImage.select(['divergence']).multiply(lastImage.select(['elevation']).subtract(lastImage.select(['fac']))).rename('k3');
  lastImage = lastImage.addBands(k3last)
  emptyLagr = emptyLagr.merge(ee.ImageCollection.fromImages([lastImage.set('MosaicYear',endYear+'/'+(endYear+1))]));

  // print(emptyLagr)

  var lagrangian = emptyLagr.select(['system:time_start','elevation']).reduce(ee.Reducer.linearFit());
  lagrangian = outlier_mask(lagrangian.select(['scale']).rename('Lagrangian'));


  var facChange = emptyLagr.select(['system:time_start','fac']).reduce(ee.Reducer.linearFit()).select(['scale']).rename('facChange');


  /////////////// BMB ////////////////

  var startEpo = emptyLagr.first().select(['system:time_start']).reduceRegion({reducer:ee.Reducer.min(), geometry:IS.geometry(), scale:500}).get('system:time_start');
  startEpo = ee.Date(ee.Number(startEpo).add(40).multiply(1000 * 60 * 60 * 24 * 365));
  // print(startEpo)

  var endEpo = lastImage.select(['system:time_start']).reduceRegion({reducer:ee.Reducer.max(), geometry:IS.geometry(), scale:500}).get('system:time_start');
  endEpo = ee.Date(ee.Number(endEpo).add(40).multiply(1000 * 60 * 60 * 24 * 365));
  // print(endEpo)

  /*
  Get surface mass balance
  OBS! Here you have to use your own surface mass balance data!
  */

  // Surface mass balance
  var smb = ee.ImageCollection('users/azcryosphere/RACMO_IMAU/RACMO_ERA5_extrapolated_bicubic');
  smb = smb.filterDate(startEpo,endEpo).mean();
  // from mm weq to m ieq
  var water2ice = ee.Number(0.917);
  var mm2m = ee.Number(0.001);
  smb = smb.divide(water2ice).multiply(mm2m).multiply(ee.Number(12));
  smb = smb.updateMask(lagrangian.mask());
  var smb = smb.reproject({
        crs: 'EPSG:3031',
        scale: 1000
      });


  var rhoWater = ee.Image.constant(1026).toFloat().select('constant');
  rhoWater = rhoWater.updateMask(lagrangian.mask());

  // Constant ice density of 917 kg m^(-3)
  var rhoIce = ee.Image.constant(917).toFloat().select('constant');
  rhoIce = rhoIce.updateMask(lagrangian.mask());




  var k3 = emptyLagr.select(['k3']).mean(); // mean(div(v)*(h-d))

  // (h-h_f)*div(V)
  var k0 = lagrangian.subtract(facChange); // Dh/Dt - Dhf/Dt
  var k1 = k0.add(k3); // Dh/Dt - Dhf/Dt + div(v)*(h-d)
  var k2 = rhoWater.divide(rhoWater.subtract(rhoIce)); // rho_w/(rho_w-rho_i)
  var kk = k1.multiply(k2); // (Dh/Dt + div(v)*(h-d)) * (rho_w/(rho_w-rho_i))
  var basalMeltRate = smb.subtract(kk.toFloat()).rename('BMB'); // a - (Dh/Dt + div(v)*(h-d)) * (rho_w/(rho_w-rho_i))

  lagrangian = lagrangian.addBands(basalMeltRate).addBands(facChange);

  // print(lagrangian)
  }
  return lagrangian.set('dt',dt,'startYear',startYear + '/' + (startYear+1), 'endYear', endYear + '/' + (endYear+1));
};


var emptyHugeLagrangian = ee.ImageCollection([]);

////////////////////////////////////////////////////////////////////
/////////////////////// CHANGABLE PARAMETERS ///////////////////////
////////////////////////////////////////////////////////////////////

var startEpoch = 2010; // First year to be considered in the trend period (min is 2010 if using the provided imageCollection which means 2010/11)
var endEpoch = 2017; // Last year to be considered in the trend period (the last data in the provided imageCollection is from 2017/18 --> endEpoch 2017)
var dt = 7; // The number of years to calculate the trend over. Max dt is endEpoch-startEpoch

////////////////////////////////////////////////////////////////////

for (var j = startEpoch; j<= endEpoch; j = j+1) {
  var output = hugeFunction(j,dt);
  emptyHugeLagrangian = emptyHugeLagrangian.merge(ee.ImageCollection.fromImages([output]));
}

print(emptyHugeLagrangian);


var year = '2010/2011';


Map.addLayer(emptyHugeLagrangian.select('Lagrangian').filterMetadata('startYear','equals',year).first().clip(IS),trendVis,'Lagrangian startEpoch ' + year);
Map.addLayer(emptyHugeLagrangian.select('BMB').filterMetadata('startYear','equals',year).first().clip(IS),meltVis,'BMB startEpoch ' + year);


var sui  = require('users/earthmapps/public:ui.js');
var colorleg = ui.Panel({
    style: {
        width: '300px',
        position: 'bottom-right',
        padding: '8px 15px'
    }
});
// Add legend colorbars to panel
colorleg.widgets().set(0, sui.colorbar(trendVis,'elevation trend [m/yr]'));
colorleg.widgets().set(1, sui.colorbar(meltVis,'basal mass balance [m/yr]'));

Map.add(colorleg);
