/*

Loose collection os examples

*/


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