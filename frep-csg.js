var ALPHA = 0.0;

CSG = function(params, attrs, func) {
	this.params = params;
	this.attrs = attrs;
	this.func = func;
	this.vertices = new Array();
	this.normals = new Array();
	this.indices = new Array();
}

CSG.prototype = {

	polygonise: function(grid, boundingBox, isosurface, callback) {

		var grid = grid?grid:50;
		var boundingBox = boundingBox?boundingBox:{min:{x:-5.0,y:-5.0,z:-5.0},max:{x:5.0,y:5.0,z:5.0}};
		var isosurface = isosurface?isosurface:0.0;

		var polygoniserWorker = new Worker('PolygoniserWorker.js')
		var that = this;

		polygoniserWorker.onmessage = function(e){
			if (e.data.msg != undefined){
				notify(e.data.msg);
			}
			if (e.data.progress != undefined){
				incrementProgress(1);
			}
			if (e.data.results != undefined){
				that.vertices = e.data.results.vertices;
				that.normals = e.data.results.normals;
				that.indices = e.data.results.indices;

				var mesh = new GL.Mesh({ normals: true, colors: true });
				var indexer = new GL.Indexer();

				for (var i = 2; i < that.indices.length; i+=3) {
					mesh.triangles.push([that.indices[i-2], that.indices[i - 1], that.indices[i]]);
				}

				mesh.vertices = that.vertices;
				mesh.normals = that.normals;
				mesh.colors = that.vertices.map(function(v){return [0.85,0.85,0.85,0.85];});
				mesh.computeWireframe();

				callback(mesh);
			}
		}

		polygoniserWorker.postMessage({'boundingBox':boundingBox, 'grid':grid, 'isosurface':isosurface, 'func': this.func.toString(), 'params':this.params, 'attrs':this.attrs})

	},
	
	toStl: function(callback){

		if (this.vertices == undefined || this.vertices.length == 0){
			notify("No vertices found, Polygonise!")
			return;
		}

		resetProgress(this.indices.length/3);
		var stlOutputWorker = new Worker('StlOutputWorker.js')
		stlOutputWorker.onmessage = function(e){
			if (e.data.msg != undefined){
				notify(e.data.msg);
			}
			if (e.data.url){
				// navigate to file, will download
				location.href = e.data.url;
				callback();
			}
			if (e.data.progress != undefined){
				incrementProgress(e.data.progress);
			}
		};
		stlOutputWorker.postMessage({'vertices':this.vertices,'indices':this.indices});

	},
	union: function(otherCsg){
		var csg = new CSG();
		var that = this;
		csg.func = function(coords) {
			var result1 = that.func(coords, that.params);
			var result2 = otherCsg.func(coords, otherCsg.params);
			return (1/(1+ALPHA))*(result2 + result1 + Math.sqrt(Math.abs(result2*result2 + result1*result1-2*ALPHA*result1*result2)));
		};
		return csg;
	},
	intersect: function(otherCsg){
		var csg = new CSG();
		var that = this;

		csg.func = function(coords) {
			var result1 = that.func(coords, that.params);
			var result2 = otherCsg.func(coords, otherCsg.params);

			return (1/(1+ALPHA))*(result2 + result1 - 
				Math.sqrt(Math.abs(result2*result2 + result1*result1 - 2*ALPHA*result1*result2)));
		};
		return csg;
	},
	subtract: function(otherCsg){
		var csg = new CSG();
		var that = this;

		csg.func = function(coords) {
			var result1 = otherCsg.func(coords, otherCsg.params);
			var result2 = that.func(coords, that.params);

			result1 = -result1;
					
			return (1/(1+ALPHA))*(result2 + result1 - 
				Math.sqrt(Math.abs(result2*result2 + result1*result1 - 2*ALPHA*result1*result2)));
		};
		return csg;
	}
};

/**Shapes **/

CSG.block = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var vertex = params.vertex;
		var dx = params.dx;
		var dy = params.dy;
		var dz = params.dz;

		var x0 = -(coords[0] - vertex[0]) * (coords[0] - (vertex[0] + dx));
		var y0 = -(coords[1] - vertex[1]) * (coords[1] - (vertex[1] + dy));
		var z0 = -(coords[2] - vertex[2]) * (coords[2] - (vertex[2] + dz));
		var i0 = x0 + y0 - Math.sqrt(x0 * x0 + y0 * y0);
		return i0 + z0 - Math.sqrt(i0 * i0 + z0 * z0);
	});
};


CSG.sphere = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var R = params.radius;
		var x = coords[0] - params.center[0];
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2];
		return (R * R) - (x * x) - (y * y) - (z * z);
	});
};

CSG.torusX = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var R = params.R;
		var r0 = params.r0;

		var x = coords[0] - params.center[0];
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2];

		return (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((y * y) + (z * z));
	});
};

CSG.gyroid = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var x = coords[0] - params.center[0];
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2];
        var r = x * x + y * y + z * z;
        var	ti = Math.abs(Math.sin(params.t / 100000)) / 10 + 0.06;
        var v = ti * r;
        return (Math.cos(x / v) * Math.sin(y / v) + Math.cos(y / v) * Math.sin(z / v) 
        	+ Math.cos(z / v) * Math.sin(x / v) + 1.0) - 0.1 * (1 - 0.016 * (r - 10 / r));
	});
};

CSG.ellipsoid = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var x = (coords[0] - params.center[0]) / params.a;
		var y = (coords[1] - params.center[1]) / params.b;
		var z = (coords[2] - params.center[2]) / params.c;
        
        return 1 - (x * x) - (y * y) - (z * z);
	});
};

CSG.cylinderY = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var x = coords[0] - params.center[0];
		var z = coords[2] - params.center[2];
        return (params.R * params.R) - (x * x) - (z * z);
	});
};

CSG.cylinderX = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2];
        return (params.R * params.R) - (y * y) - (z * z);
	});
};

CSG.cylinderZ = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var x = coords[0] - params.center[0];
		var y = coords[1] - params.center[1];
        return (params.R * params.R) - (x * x) - (y * y);
	});
};



CSG.heart = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var x = (coords[0] - params.center[0]);
		var y = (coords[1] - params.center[1]);
		var z = (coords[2] - params.center[2]);

		var pow = Math.pow
		return pow(pow(x,2)+(9/4)*pow(y,2)+pow(z,2)-1,3) - pow(x,2)*pow(z,3)-(9/80)*pow(y,2)*pow(z,3);
	});
};
