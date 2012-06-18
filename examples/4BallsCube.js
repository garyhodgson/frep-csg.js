var a = CSG.sphere({center:[-2,-2,-2], radius:2});
var a1 = CSG.sphere({center:[2,2,2], radius:2});
var a2 = CSG.sphere({center:[-2,2,2], radius:2});
var a3 = CSG.sphere({center:[-2,-2,2], radius:2});

var a4 = CSG.sphere({center:[2,2,-2], radius:2});
var a5 = CSG.sphere({center:[2,-2,-2], radius:2});
var a6 = CSG.sphere({center:[-2,2,-2], radius:2});
var a7 = CSG.sphere({center:[2,-2,2], radius:2});

var b = CSG.block({vertex:[-2,-2,-2], dx:4,dy:4,dz:4});

return  b.union(a).union(a1).union(a2).union(a3).union(a4).union(a5).union(a6).union(a7);