// node-sgp4 test cases

// Requirements
var expect      = require('expect.js');
var fs          = require('fs');
var testData    = require('./testData.json');
var SGP4        = require('../');

function calcErr(theoretical, experimental) {
	return Math.abs((experimental - theoretical)/theoretical) * 100;
}

describe('SGP4', function() {
	describe('should return values with less than 5% error', function() {
		Object.keys(testData).forEach(function(testName) {
			var caseData = testData[testName];
			var satrec = SGP4.twoline2rv(caseData.tle.line1, caseData.tle.line2, SGP4.wgs84());

			caseData.results.forEach(function(theoreticalResult) {
				var experimentalResult = SGP4.sgp4(satrec, theoreticalResult.time);

				it("[" + testName + "][Position] Time: " + theoreticalResult.time, function() {
					expect(calcErr(theoreticalResult.position.x, experimentalResult.position.x)).to.be.lessThan(5);
					expect(calcErr(theoreticalResult.position.y, experimentalResult.position.y)).to.be.lessThan(5);
					expect(calcErr(theoreticalResult.position.z, experimentalResult.position.z)).to.be.lessThan(5);
				});

				it("[" + testName + "][Velocity] Time: " + theoreticalResult.time, function() {
					expect(calcErr(theoreticalResult.velocity.x, experimentalResult.velocity.x)).to.be.lessThan(5);
					expect(calcErr(theoreticalResult.velocity.y, experimentalResult.velocity.y)).to.be.lessThan(5);
					expect(calcErr(theoreticalResult.velocity.z, experimentalResult.velocity.z)).to.be.lessThan(5);
				});
			});		
		});
	});
});