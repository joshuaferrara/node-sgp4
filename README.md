# node-sgp4

SGP4 is one of five simplified pertubation models used to calculate various orbital information about satellites. This library was ported from python-sgp4 and uses some functions from satellite-js to convert coordinates.

If you need any support, click the badge below. *Support is not guaranteed and I provide it only when I am not busy.*
[![Join the chat at https://gitter.im/joshuaferrara/node-sgp4](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/joshuaferrara/node-sgp4?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![NPM](https://nodei.co/npm/sgp4.png?downloads=true&stars=true)](https://nodei.co/npm/sgp4/)
### Installation
`npm install sgp4`
### Example | [Google Maps Demo](http://tracking.ferrara.space/)
An example is provided in the top directory named `example.js`. To run the example, simply use the command `node example`. This example program will print the latitude and longitude of the ISS at the current point in time based off of TLE data from 14 March 2015.
### Notes
Most of the methods are self explanatory. The example does a great job of showing how to get different types of coordinates. There are many functions at the bottom of sgp4.js that allow conversion between different coordinate models.
### License
MIT
