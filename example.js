var SGP4 = require('./');

// Sample ISS TLE Data from 14 September 2015
var issLine1 = "1 25544U 98067A   15256.76500793  .00042673  00000-0  62367-3 0  9995";
var issLine2 = "2 25544  51.6468   1.7516 0001021 349.4970 110.1741 15.55194220961827";

// Create a satellite record
var issSatRec = SGP4.twoline2rv(issLine1, issLine2, SGP4.wgs84());

// This will print some info every 3/4th second
function printPosition() {
    // Current time
    var now = new Date();
    
    // This will contain ECI (http://en.wikipedia.org/wiki/Earth-centered_inertial) coordinates of position and velocity of the satellite
    var positionAndVelocity = SGP4.propogate(issSatRec, now.getUTCFullYear(), now.getUTCMonth() + 1, now.getUTCDate(), now.getUTCHours(), now.getUTCMinutes(), now.getUTCSeconds());
    
    // Prints ECI coordinates
    console.log(positionAndVelocity);
    
    // GMST required to get Lat/Long
    var gmst = SGP4.gstimeFromDate(now.getUTCFullYear(), now.getUTCMonth() + 1, now.getUTCDate(), now.getUTCHours(), now.getUTCMinutes(), now.getUTCSeconds());
    
    // Geodetic coordinates
    var geodeticCoordinates = SGP4.eciToGeodetic(positionAndVelocity.position, gmst);
    
    // Coordinates in degrees
    var longitude = SGP4.degreesLong(geodeticCoordinates.longitude);
    var latitude = SGP4.degreesLat(geodeticCoordinates.latitude);
    
    // Prints latitude of longitude of ISS
    console.log('Lat/Long: ' + latitude + ' ' + longitude);
    
    // Prints current speed of satellite in km/s
    console.log('Velocity: ' + geodeticCoordinates.velocity + ' km/s');
    
    // Prints orbital period of satellite in minutes
    // 2pi * sqrt(Relative Height / Gravity of Earth * Mass of Earth)
    console.log('Oribital Period: ' + ((2 * Math.PI) * (geodeticCoordinates.height + 6378.135)) * (Math.sqrt((geodeticCoordinates.height + 6378.135)/398600.8)) / 60);
    
    var observerPos = {
        longitude: -117.61199249999999 * SGP4.deg2rad,
        latitude: 33.4269728 * SGP4.deg2rad,
        height: 1
    };

    var satEcf = SGP4.eciToEcf(positionAndVelocity.position, gmst);
    var lookAngles = SGP4.topocentricToLookAngles(SGP4.topocentric(observerPos, satEcf));

    console.log("Azimuth: " + lookAngles.azimuth * SGP4.rad2deg);
    console.log("Elevation: " + lookAngles.elevation * SGP4.rad2deg);

    // Call printPosition in 750 ms
    setTimeout(printPosition, 750);
}
printPosition();