/*

Loose collection of examples to cut and paste

*/

// blank script 

return new CSG({center:[0,0,0]}, {}, function(coords, params, attrs){

	var x = coords[0] - params.center[0];
	var y = coords[1] - params.center[1];
	var z = coords[2] - params.center[2];
        

	var result = 

	return [result, attrs];
});

/************************************************************************************************************************/


// row of torus structures
var r = 2;
r1 = 0.3;

var struct = function(center, r, r1){
  return CSG.torusX({center:center, R:r, r0:r1}).union(CSG.torusZ({center:center, R:r, r0:r1}).union(CSG.torusY({center:center, R:r, r0:r1})))};

var group = struct([0,0,0],r,r1);

for (var i = 5; i < 25; i+=5){

  group = group.union(struct([0,0,i],r,r1));
} 

return group;


/************************************************************************************************************************/

// row of torus structures  -with color
var r = 2;
r1 = 0.3;

var struct = function(center, r, r1, attrs){
  return CSG.torusX({center:center, R:r, r0:r1},attrs).union(CSG.torusZ({center:center, R:r, r0:r1},attrs).union(CSG.torusY({center:center, R:r, r0:r1},attrs)))};

var colors = [[0,0,1],[1,1,0],[1,0,1]]

var group = struct([0,0,0],r,r1,{color:[1,0,0]});

for (var i = 1; i < 4; i++){
  var s = struct([0,0,i*5],r,r1, {color:colors[i-1]})
  group = group.union(s);
} 

return group;

/************************************************************************************************************************/

// screw

var height = 7;

var c = CSG.ellipticCylinderZ({center:[0,0,0], a:1,b:2}, {color:[0,0.9,0.9]});

var b = CSG.block({vertex:[-2,-2,0], dx:4,dy:4,dz:height});

var tube = c.intersect(b);

return tube.twistZ({z1:-1, z2:1, theta1:0, theta2:-3.14});


/************************************************************************************************************************/


// blobby ball
var bc = [[-7.0, -7, -7],[-6.0,-4,-6],[-4.0, -7, -4],[-2.0, -4, -2],[0, -2, 0],[7, 1,4]];
var a = [0.7, 1, 1, 0.3, 1, 1];
var b = [3, 1, 1, 1, 1, 0.5];

return CSG.blobbyball({bc:bc,a:a,b:b,T:0.05})



/************************************************************************************************************************/

// convline - with stretch and twist

var begin = [[-8.0, 0.0, 0.0],
             [0.0, -8.0, 0.0],
             [0.0, 0.0, -8.0]];
var end = [[8.0, 0.0, 0.0],
       	   [0.0, 8.0, 0.0],
           [0.0, 0.0, 8.0]];
var S = [0.65, 0.65, 1];

return CSG.convLine({T:0.5, S:S, begin:begin,end:end}).stretch({x0:[0,0,-8],sx:1,sy:1,sz:1.5}).twistZ({z1:-8,z2:10,theta1:1,theta2:3});


/************************************************************************************************************************/
//convarc

var centers = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]];
var radius =  [8.5, 8.5, 8.5];
var theta =   [360.0, 360.0, 360.0];
var rotate_axes = [[0.0, 0.0, 1.0],[0.0, 1.0, 0.0],[1.0, 0.0, 0.0]];
var rotate_angles = [0.0, 90.0, 90.0];
var S =       [0.65, 0.65, 0.65];

return CSG.convArc({centers:centers,radius:radius,theta:theta,rotate_axes:rotate_axes,rotate_angles:rotate_angles,S:S,T:0.22});

/************************************************************************************************************************/


// scaled shifted egg

var scaleFactor = 0.5;

var egg = CSG.sphere({center:[0,0,0], radius:3}, {color:[1,0,0]}).taperZ({z1:-2,z2:4,s1:1,s2:0.5}).scale({sx:scaleFactor,sy:scaleFactor,sz:scaleFactor});

return egg.shift({dx:3,dy:3,dz:3})