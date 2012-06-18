var sphereradius = 10;

var lattice = new CSG({sphereradius:sphereradius}, {color:[0.7,0.7,0.9]}, function(coords, params, attrs) {	

	var sphere = CSG.sphere({center:[0,0,0],radius:params.sphereradius});

	var sphereval = sphere.call(coords)[0];

	var lx = 0.55;
	var ly = lx;
	var lz = 0;

	var minX = 0
	var maxX = 105
	var mapMinX = 6.25
	var mapMaxX = 1

	var qx = (sphereval-minX) / (maxX-minX) * (mapMaxX-mapMinX) + mapMinX;
	var qy = (sphereval-minX) / (maxX-minX) * (mapMaxX-mapMinX) + mapMinX;
	var qz = 2.5;

	var px = -1.5;
	var py = px;
	var pz = 0;	

	var sx = Math.sin((qx * coords[0]) + px) - lx;
	var sy = Math.sin((qy * coords[1]) + py) - ly;
	var sz = Math.sin((qz * coords[2]) + pz) - lz;

	var rx = Math.min(sy,sz);
	var ry = Math.min(sx,sz);
	var rz = Math.min(sx,sy);

	var result =  Math.max(Math.max(rx, ry), rz);

	return [result, attrs];
});

var sphere = CSG.sphere({center:[0,0,0],radius:sphereradius}, {color:[1,1,1]});

var spherehollow = CSG.sphere({center:[0,0,0],radius:sphereradius-1});

var shell = sphere.subtract(spherehollow);

var block = CSG.block({vertex:[-10,-10,0],dx:20,dy:20,dz:10});

var support = lattice.intersect(spherehollow);

var supportshell = shell.union(support)

return supportshell.subtract(block);