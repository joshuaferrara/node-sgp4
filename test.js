require('mocha');
var SGP4 = require('./');
var assert = require('assert');
var issLine1 = "1 25544U 98067A   15073.52509288  .00016392  00000-0  24500-3 0  9990";
var issLine2 = "2 25544  51.6454 196.8789 0008842 100.5316 358.5831 15.55127532933300";

var firstOfYear = new Date(2015, 0, 1);

describe("SGP4", function() {
	describe(".twoline2rv(line1, line2, gravconst)", function() {
		it("should return an object with satnum of 25544", function() {
			var issSatRec = SGP4.twoline2rv(issLine1, issLine2, SGP4.wgs84());
			assert.equal(issSatRec.satnum, "25544");
		});
	});

	describe(".propogate(satrec, year, month, day, hour, minute, second)", function() {
		it("should return correct ECI coordinates and velocity of satellite at date Jan 1, 2015 00:00:00", function() {
			var issSatRec = SGP4.twoline2rv(issLine1, issLine2, SGP4.wgs84());
			var positionAndVelocity = SGP4.propogate(issSatRec, firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds());

			assert.equal(positionAndVelocity.position.x, 370.25446753024704);
			assert.equal(positionAndVelocity.position.y, 4385.686016393587);
			assert.equal(positionAndVelocity.position.z, -5168.413191512966);

			assert.equal(positionAndVelocity.velocity.x, -7.4483937043258495);
			assert.equal(positionAndVelocity.velocity.y, -1.055132910776903);
			assert.equal(positionAndVelocity.velocity.z, -1.4380797484368044);
		});
	});

	describe(".gstimeFromDate(year, mon, day, hr, minute, sec)", function() {
		it("should return correct sidereel time for Jan 1, 2015 00:00:00", function() {
			var gmst = SGP4.gstimeFromDate(firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds())
			assert.equal(gmst, 3.85121337656868);
		});
	});

	describe(".eciToGeodetic(eci_coords, gmst)", function() {
		it("should return correct geodetic coordinates for time Jan 1, 2015 00:00:00", function() {
			var issSatRec = SGP4.twoline2rv(issLine1, issLine2, SGP4.wgs84());
			var positionAndVelocity = SGP4.propogate(issSatRec, firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds());
			var gmst = SGP4.gstimeFromDate(firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds());
			var geodeticCoordinates = SGP4.eciToGeodetic(positionAndVelocity.position, gmst);

			assert.equal(geodeticCoordinates.longitude, -2.3646407195802626);
			assert.equal(geodeticCoordinates.latitude, -0.8684966567825487);
			assert.equal(geodeticCoordinates.height, 422.7895363497473);
			assert.equal(geodeticCoordinates.velocity, 7.655703497788085);
		});
	});

	describe(".degreesLong(geodeticLong)", function() {
		it("should return correct longitude in degrees", function() {
			var issSatRec = SGP4.twoline2rv(issLine1, issLine2, SGP4.wgs84());
			var positionAndVelocity = SGP4.propogate(issSatRec, firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds());
			var gmst = SGP4.gstimeFromDate(firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds());
			var geodeticCoordinates = SGP4.eciToGeodetic(positionAndVelocity.position, gmst);

			var longitude = SGP4.degreesLong(geodeticCoordinates.longitude);
			assert.equal(longitude, -135.48393329672706);
		});
	});

	describe(".degreesLat(geodeticLat)", function() {
		it("should return correct latitude in degrees", function() {
			var issSatRec = SGP4.twoline2rv(issLine1, issLine2, SGP4.wgs84());
			var positionAndVelocity = SGP4.propogate(issSatRec, firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds());
			var gmst = SGP4.gstimeFromDate(firstOfYear.getUTCFullYear(), firstOfYear.getUTCMonth() + 1, firstOfYear.getUTCDate(), firstOfYear.getUTCHours(), firstOfYear.getUTCMinutes(), firstOfYear.getUTCSeconds());
			var geodeticCoordinates = SGP4.eciToGeodetic(positionAndVelocity.position, gmst);

			var latitude = SGP4.degreesLat(geodeticCoordinates.latitude);
			assert.equal(latitude, -49.76119295486205);
		});
	});
});