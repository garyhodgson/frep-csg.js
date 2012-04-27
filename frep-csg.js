var ALPHA = 0.0;

var MAX_BOUNDING_BOX = {min:{x:-200.0,y:-200.0,z:-200.0},max:{x:200.0,y:200.0,z:200.0}};
var MIN_BOUNDING_BOX = {min:{x:-5.0,y:-5.0,z:-5.0},max:{x:5.0,y:5.0,z:5.0}};

CSG = function(params, attrs, func) {
	this.params = params;
	this.attrs = attrs;
	this.func = func;
	this.funcDef = func.toString();	
	this.vertices = new Array();
	this.normals = new Array();
	this.indices = new Array();
}

CSG.prototype = {
	call: function(coords){
		return this.func(coords, this.params, this.attrs);
	},
	polygonise: function(gridSize, boundingBox, isosurface, numWorkers, callback) {

		var gridSize = gridSize?gridSize:50;
		var grid = {x:gridSize,y:gridSize,z:gridSize}
		var boundingBox = boundingBox?boundingBox:{min:{x:-5.0,y:-5.0,z:-5.0},max:{x:5.0,y:5.0,z:5.0}};
		var isosurface = isosurface?isosurface:0.0;

		var numWorkers = numWorkers?numWorkers:4;
		var polygoniserWorkers = new Array(numWorkers);
		var resultsCount = 0;
		
		this.verticesArray = new Array(numWorkers);
		this.normalsArray = new Array(numWorkers);
		this.indicesArray = new Array(numWorkers);
		var that = this;

		var division = (Math.abs(boundingBox.min.z)+boundingBox.max.z)/numWorkers;

		for (var i = 0; i <numWorkers; i++){
		
			polygoniserWorkers[i] = new Worker('PolygoniserWorker.js')			

			polygoniserWorkers[i].onmessage = function(e){
				
				if (e.data.msg != undefined) notify(e.data.msg);
				if (e.data.progress != undefined) incrementProgress(1);

				if (e.data.results != undefined) {

					var workerId = parseInt(e.data.worker)

					that.verticesArray[workerId] = e.data.results.vertices;
					that.normalsArray[workerId] = e.data.results.normals;
					that.indicesArray[workerId] = e.data.results.indices;

					resultsCount++;
					
					if (resultsCount == numWorkers){
						var offset = 0

						for(var h=0; h < numWorkers; h++){
							that.vertices = that.vertices.concat(that.verticesArray[h]);
							that.normals = that.normals.concat(that.normalsArray[h]);

							var indexArray = that.indicesArray[h];
							for (var k in indexArray){
								var newIndex = parseInt(indexArray[k]) + offset
								that.indices.push(newIndex)
							}
							offset += that.verticesArray[h].length

						}

						var mesh = new GL.Mesh({ normals: true, colors: true });

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
			}
			var bb = $.extend(true, {}, boundingBox);
			bb.min.z = boundingBox.min.z + (division*i);
			bb.max.z = bb.min.z + division;
			grid.z = gridSize/numWorkers

			polygoniserWorkers[i].postMessage({'worker':i, 'boundingBox':bb, 'grid':grid, 'isosurface':isosurface, 'funcDef': this.funcDef, 'params':this.params, 'attrs':this.attrs})
		}
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
		var params = [this.params, otherCsg.params];
		var attrs =  [this.attrs, otherCsg.attrs];
		var funcDef = "	var result1 = "+this.func.toString()+"(coords,params[0],attrs[0]);\
						var result2 = "+otherCsg.func.toString()+"(coords,params[1],attrs[1]); \
						return (1/(1+ALPHA))*(result2 + result1 + Math.sqrt(Math.abs(result2*result2 + result1*result1-2*ALPHA*result1*result2)));";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	add: function(otherCsg){
		var params = [this.params, otherCsg.params];
		var attrs =  [this.attrs, otherCsg.attrs];
		var funcDef = "	var result1 = "+this.func.toString()+"(coords,params[0],attrs[0]);\
						var result2 = "+otherCsg.func.toString()+"(coords,params[1],attrs[1]); \
						return result1 + result2";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	intersect: function(otherCsg){
		var params = [this.params, otherCsg.params];
		var attrs =  [this.attrs, otherCsg.attrs];
		var funcDef = "	var result1 = "+this.func.toString()+"(coords,params[0],attrs[0]);\
						var result2 = "+otherCsg.func.toString()+"(coords,params[1],attrs[1]); \
						return (1/(1+ALPHA))*(result2 + result1 - Math.sqrt(Math.abs(result2*result2 + result1*result1 - 2*ALPHA*result1*result2)));";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	subtract: function(otherCsg){
		var params = [this.params, otherCsg.params];
		var attrs =  [this.attrs, otherCsg.attrs];
		var funcDef = "	var result1 = "+this.func.toString()+"(coords,params[0],attrs[0]);\
						var result2 = "+otherCsg.func.toString()+"(coords,params[1],attrs[1]); \
						result1 = -result1; \
						return (1/(1+ALPHA))*(result2 + result1 - Math.sqrt(Math.abs(result2*result2 + result1*result1 - 2*ALPHA*result1*result2)));";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	rotateZ: function(theta){
		var params = [this.params,{'theta':theta}];
		var attrs =  [this.attrs];
		var funcDef = "	var theta = params[1].theta; \
						var ct = Math.cos(theta); \
						var st = Math.sin(theta); \
						var xr = coords[0] * ct + coords[1] * st; \
						var yr = -coords[0] * st + coords[1] * ct; \
						coords[0] = xr; \
						coords[1] = yr; \
						return "+this.func.toString()+"(coords,params[0],attrs[0]);\
						";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	//x1, x2 - end points of x-interval
 	//theta1, theta2 - rotation angles in radians for end points
	twistX: function(theta1, theta2, x1, x2){
		var params = [this.params, {'theta1':theta1,'theta2':theta2,'x1':x1,'x2':x2}];
		var attrs =  [this.attrs];
		var funcDef = "	var theta1 = params[1].theta1; \
						var theta2 = params[1].theta2; \
						var x1 = params[1].x1; \
						var x2 = params[1].x2; \
						var t = (coords[0]-x1)/(x2-x1); \
						var theta = (1-t)*theta1 + t*theta2; \
						var ct = Math.cos(theta); \
						var st = Math.sin(theta); \
						var yr = coords[1] * ct + coords[2] * st; \
						var zr = -coords[1] * st + coords[2] * ct; \
						coords[1] = yr; \
						coords[2] = zr; \
						return "+this.func.toString()+"(coords,params[0],attrs[0]);\
						";
		return new CSG(params, attrs, funcDef);
	}
};

/**Shapes **/
CSG.block = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
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
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.radius;
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
		return (R * R) - (x * x) - (y * y) - (z * z); 
		});
};

CSG.sphereOld = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.radius; 
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
		return (R * R) - (x * x) - (y * y) - (z * z); 
	});
};

CSG.torusX = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.R; 
		var r0 = params.r0; 
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
		return (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((y * y) + (z * z)); 
	});
};

CSG.torusY = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.R; 
		var r0 = params.r0; 
		var x = coords[0] - params.center[0];
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2]; 
		return (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((x * x) + (z * z)); 
	});
};

CSG.torusZ = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.R; 
		var r0 = params.r0; 
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
		return (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((x * x) + (y * y)); 
	});
};

CSG.gyroid = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
        var r = x * x + y * y + z * z; 
        var	ti = Math.abs(Math.sin(params.t / 100000)) / 10 + 0.06; 
        var v = ti * r; 
        return (Math.cos(x / v) * Math.sin(y / v) + Math.cos(y / v) * Math.sin(z / v)  + Math.cos(z / v) * Math.sin(x / v) + 1.0) - 0.1 * (1 - 0.016 * (r - 10 / r)); 
	});
};

CSG.ellipsoid = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = (coords[0] - params.center[0]) / params.a; 
		var y = (coords[1] - params.center[1]) / params.b; 
		var z = (coords[2] - params.center[2]) / params.c; 
        return 1 - (x * x) - (y * y) - (z * z); 
	});
};

CSG.cylinderY = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = coords[0] - params.center[0]; 
		var z = coords[2] - params.center[2]; 
        return (params.R * params.R) - (x * x) - (z * z); 
	});
};

CSG.cylinderX = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2]; 
        return (params.R * params.R) - (y * y) - (z * z); 
	});
};

CSG.cylinderZ = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
        return (params.R * params.R) - (x * x) - (y * y); 
	});
};



CSG.heart = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = (coords[0] - params.center[0]); 
		var y = (coords[1] - params.center[1]); 
		var z = (coords[2] - params.center[2]); 
		var pow = Math.pow; 
		return pow(pow(x,2)+(9/4)*pow(y,2)+pow(z,2)-1,3) - pow(x,2)*pow(z,3)-(9/80)*pow(y,2)*pow(z,3); 
	});
};
