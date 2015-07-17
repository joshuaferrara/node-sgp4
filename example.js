var SGP4 = require('./');

// Sample ISS TLE Data
var issLine1 = "1 25544U 98067A   15073.52509288  .00016392  00000-0  24500-3 0  9990";
var issLine2 = "2 25544  51.6454 196.8789 0008842 100.5316 358.5831 15.55127532933300";

// Create a satellite record
var issSatRec = SGP4.twoline2rv(issLine1, issLine2, SGP4.wgs84());

// This will print some info every 3/4th second
function printPosition() {
    // Current time
    var now = new Date();
    
    // This will contain ECI (http://en.wikipedia.org/wiki/Earth-centered_inertial) coordinates of position and velocity of the satellite
    var positionAndVelocity = SGP4.propogate(issSatRec, now.getUTCFullYear(), now.getUTCMonth() + 1, now.getUTCDate(), now.getUTCHours(), now.getUTCMinutes(), now.getUTCSeconds());
    
    // Prints ECI coordinates
    //console.log(positionAndVelocity);
    
    // GMST required to get Lat/Long
    var gmst = SGP4.gstimeFromDate(now.getUTCFullYear(), now.getUTCMonth() + 1, now.getUTCDate(), now.getUTCHours(), now.getUTCMinutes(), now.getUTCSeconds());
    //console.log(gmst);

    // Geodetic coordinates
    var geodeticCoordinates = SGP4.eciToGeodetic(positionAndVelocity.position, gmst);
    //console.log(geodeticCoordinates);

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
    
    // Call printPosition in 750 ms
    setTimeout(printPosition, 750);
}
printPosition();
