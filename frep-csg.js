var ALPHA = 0.0;
var EPS = 0.0;

var MAX_BOUNDING_BOX = {min:{x:-200.0,y:-200.0,z:-200.0},max:{x:200.0,y:200.0,z:200.0}};
var MIN_BOUNDING_BOX = {min:{x:-5.0,y:-5.0,z:-5.0},max:{x:5.0,y:5.0,z:5.0}};
var DEFAULT_GRID_SIZE = 50;

CSG = function(params, attrs, func) {
	this.params = params;
	this.attrs = attrs||{};
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
	polygonise: function(grid, boundingBox, isosurface, numWorkers, callback) {

		var grid = grid || {x:DEFAULT_GRID_SIZE,y:DEFAULT_GRID_SIZE,z:DEFAULT_GRID_SIZE}
		var boundingBox = boundingBox || {min:{x:-5.0,y:-5.0,z:-5.0},max:{x:5.0,y:5.0,z:5.0}};
		var isosurface = isosurface || 0.0;

		var numWorkers = numWorkers?numWorkers:4;
		var polygoniserWorkers = new Array(numWorkers);
		var resultsCount = 0;
		
		this.verticesArray = new Array(numWorkers);
		this.normalsArray = new Array(numWorkers);
		this.indicesArray = new Array(numWorkers);
		this.colorsArray = new Array(numWorkers);
		this.vertices = new Array();
		this.normals = new Array();
		this.indices = new Array();
		this.colors = new Array();
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
					that.colorsArray[workerId] = e.data.results.colors;

					resultsCount++;
					
					if (resultsCount == numWorkers){

						var offset = 0

						for(var h=0; h < numWorkers; h++){
							that.vertices = that.vertices.concat(that.verticesArray[h]);
							that.normals = that.normals.concat(that.normalsArray[h]);
							that.colors = that.colors.concat(that.colorsArray[h]);

							var indexArray = that.indicesArray[h];
							for (var k in indexArray){
								var newIndex = parseInt(indexArray[k]) + offset
								that.indices.push(newIndex)
							}
							offset += that.verticesArray[h].length
						}

						notify("Vertices: "+that.vertices.length+ "; Normals: "+ that.normals.length+ "; Indices: "+ that.indices.length+ "; Colors: "+ that.colors.length+";");

						mesh = new GL.Mesh({ normals: true, colors: true });

						for (var i = 2; i < that.indices.length; i+=3) {
							mesh.triangles.push([that.indices[i-2], that.indices[i - 1], that.indices[i]]);
						}

						mesh.vertices = that.vertices;
						mesh.normals = that.normals;
						mesh.colors = that.colors;
						mesh.computeWireframe();

						callback(mesh);
					}
				}
			}
			var bb = $.extend(true, {}, boundingBox);
			bb.min.z = boundingBox.min.z + (division*i);
			bb.max.z = bb.min.z + division;

			var subGrid = {x:grid.x, y:grid.y, z:(grid.z/numWorkers)};
			polygoniserWorkers[i].postMessage({'worker':i, 'boundingBox':bb, 'grid':subGrid, 'isosurface':isosurface, 'funcDef': this.funcDef, 'params':this.params, 'attrs':this.attrs})
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
		var funcDef = "	var f1 = "+this.func.toString()+"(coords,params[0],attrs[0]);\
						var result1 = f1[0]; \
						var f2 = "+otherCsg.func.toString()+"(coords,params[1],attrs[1]); \
						var result2 = f2[0]; \
						var result = (1/(1+ALPHA))*(result2 + result1 + Math.sqrt(Math.abs(result2*result2 + result1*result1-2*ALPHA*result1*result2))); \
								\
						var attrs1 = f1[1];\
						var attrs2 = f2[1];\
						var attrs = {}; \
								\
						var col1 = attrs1.color || [1,1,1]; \
						var col2 = attrs2.color || [1,1,1]; \
						var unionCol = [1,1,1]; \
						if (result2 >= 0.0 && result1 >= 0.0) { \
							unionCol[0] = (col1[0] + col2[0]) / 2; \
							unionCol[1] = (col1[1] + col2[1]) / 2; \
							unionCol[2] = (col1[2] + col2[2]) / 2; \
						} else if (result1 > result2){ \
							unionCol = col1; \
						} else { \
							unionCol = col2; \
						} \
						attrs.color = unionCol; \
						return [result, attrs];";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	intersect: function(otherCsg){
		var params = [this.params, otherCsg.params];
		var attrs =  [this.attrs, otherCsg.attrs];
		var funcDef = "	var f1 = "+this.func.toString()+"(coords,params[0],attrs[0]);\
						var result1 = f1[0]; \
						var f2 = "+otherCsg.func.toString()+"(coords,params[1],attrs[1]); \
						var result2 = f2[0]; \
						var result = (1/(1+ALPHA))*(result2 + result1 - Math.sqrt(Math.abs(result2*result2 + result1*result1 - 2*ALPHA*result1*result2))); \
						var attrs = f1[1]; \
						_.extend(true, attrs, f2[1]); \
						return [result, attrs]";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	subtract: function(otherCsg){
		var params = [this.params, otherCsg.params];
		var attrs =  [this.attrs, otherCsg.attrs];
		var funcDef = "	var f1 = "+this.func.toString()+"(coords,params[0],attrs[0]);\
						var result1 = f1[0]; \
						var f2 = "+otherCsg.func.toString()+"(coords,params[1],attrs[1]); \
						var result2 = f2[0]; \
						result1 = -result1; \
						var result = (1/(1+ALPHA))*(result2 + result1 - Math.sqrt(Math.abs(result2*result2 + result1*result1 - 2*ALPHA*result1*result2))); \
						var attrs = f1[1]; \
						_.extend(true, attrs, f2[1]); \
						return [result, attrs]; \
						";
		return new CSG(params, attrs, new Function('coords', 'params', 'attrs', funcDef));
	},
	//	RotateX
	//	Definition: inverse mapping
    //		y'=y*cos(theta)+z*sin(theta)
    //		z'=-y*sin(theta)+z*cos(theta)
    //	Parameters:
    //		theta - rotation angle in radians
	rotateX: function(p){
		return CSG.methodFactory(this, p, 
			"	var theta = params[1].theta; \
				var ct = Math.cos(theta); \
				var st = Math.sin(theta); \
				var yr = coords[1] * ct + coords[2] * st;\
				var zr = -coords[1] * st + coords[2] * ct;\
				coords[1] = yr;\
				coords[2] = zr;\
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
				");
	},
	//	RotateY
	//	Definition: inverse mapping
    //		z'=z*cos(theta)+x*sin(theta)
    //		x'=-z*sin(theta)+z*cos(theta)
    //	Parameters:
    //		theta - rotation angle in radians
	rotateY: function(p){
		return CSG.methodFactory(this, p, "	var theta = params[1].theta; \
						var ct = Math.cos(theta); \
						var st = Math.sin(theta); \
						var zr = coords[2] * ct + coords[0] * st; \
						var xr = -coords[2] * st + coords[0] * ct; \
						coords[0] = xr; \
						coords[2] = zr; \
						return "+this.func.toString()+"(coords,params[0],attrs[0]);\
						");
	},
	//	RotateZ
	//	Definition: inverse mapping
    //		x'=x*cos(theta)+y*sin(theta)
    //		y'=-x*sin(theta)+y*cos(theta)
    //	Parameters:
    //		theta - rotation angle in radians
	rotateZ: function(p){
		return CSG.methodFactory(this, p, "	var theta = params[1].theta; \
						var ct = Math.cos(theta); \
						var st = Math.sin(theta); \
						var xr = coords[0] * ct + coords[1] * st; \
						var yr = -coords[0] * st + coords[1] * ct; \
						coords[0] = xr; \
						coords[1] = yr; \
						return "+this.func.toString()+"(coords,params[0],attrs[0]);\
						");
	},
	//	TwistX
	//	Definition: inverse mapping
    //		t = (x-x1)/(x2-x1)
    //		theta = (1-t)*theta1 + t*theta2
    //		y'=y*cos(theta)+z*sin(theta)
    //		z'=-y*sin(theta)+z*cos(theta)
	//	Parameters:
	//		x1, x2 - end points of x-interval
 	//		theta1, theta2 - rotation angles in radians for end points
	twistX: function(p){
		return CSG.methodFactory(this, p, "	var theta1 = params[1].theta1; \
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
						");
	},
	//	TwistY
	//	Definition: inverse mapping
	//		t = (y-y1)/(y2-y1)
	//		theta = (1-t)*theta1 + t*theta2
	//		z'=z*cos(theta)+x*sin(theta)
	//		x'=-z*sin(theta)+x*cos(theta)
	//	Parameters:
	//		y1, y2 - end points of y-interval
	//		theta1, theta2 - rotation angles in radians for end points
	twistY: function(p){
		return CSG.methodFactory(this, p, 
			"	var theta1 = params[1].theta1; \
				var theta2 = params[1].theta2; \
				var y1 = params[1].y1; \
				var y2 = params[1].y2; \
				var t = (coords[1]-y1)/(y2-y1); \
				var theta = (1-t)*theta1 + t*theta2; \
				var ct = Math.cos(theta); \
				var st = Math.sin(theta); \
				var zr = coords[2] * ct + coords[0] * st; \
				var xr = -coords[2] * st + coords[0] * ct; \
				coords[0] = xr; \
				coords[2] = zr; \
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
				");
	},
	//	TwistZ
	//	Definition: inverse mapping
	//		t = (z-z1)/(z2-z1)
	//		theta = (1-t)*theta1 + t*theta2
	//		x'=x*cos(theta)+y*sin(theta)
	//		y'=-x*sin(theta)+y*cos(theta)
	//	Parameters:
	//		z1, z2 - end points of z-interval
	//		theta1, theta2 - rotation angles in radians for end points
	twistZ: function(p){
		return CSG.methodFactory(this, p, 
			"	var theta1 = params[1].theta1; \
				var theta2 = params[1].theta2; \
				var z1 = params[1].z1; \
				var z2 = params[1].z2; \
				var t = (coords[2]-z1)/(z2-z1); \
				var theta = (1-t)*theta1 + t*theta2; \
				var ct = Math.cos(theta); \
				var st = Math.sin(theta); \
				var xr = coords[0] * ct + coords[1] * st; \
				var yr = -coords[0] * st + coords[1] * ct; \
				coords[0] = xr; \
				coords[1] = yr; \
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
				");
	},
	//	Shift
	//	Definition: x'=x+dx
	//	Parameters:
 	//		dx,dy,dz - scaling factors along axes
	shift: function(p){
		return CSG.methodFactory(this, p, "var dx = params[1].dx; \
						var dy = params[1].dy; \
						var dz = params[1].dz; \
						coords[0] = coords[0] - dx; \
						coords[1] = coords[1] - dy; \
						coords[2] = coords[2] - dz; \
						return "+this.func.toString()+"(coords,params[0],attrs[0]);\
						");
	},
	//	Stretch
	//	Definition: x'=x0+(x-x0)/scale  (inverse mapping)
	//	Parameters:
	//		x0 - reference point for stretching
 	//		sx,sy,sz - scaling factors along axes
 	stretch: function(p){
		return CSG.methodFactory(this, p, 
			"	var x0 = params[1].x0; \
				var sx = params[1].sx; \
				var sy = params[1].sy; \
				var sz = params[1].sz; \
				coords[0] = x0[0] + (coords[0] - x0[0]) / sx; \
				coords[1] = x0[1] + (coords[1] - x0[1]) / sy; \
				coords[2] = x0[2] + (coords[2] - x0[2]) / sz; \
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
			");
	},
	//	Scale
	//	Definition: x'=sx*x
	//	Parameters:
	//		sx,sy,sz - scaling factors along axes
	scale: function(p){
		return CSG.methodFactory(this, p, 
			"	coords[0] = coords[0] / params[1].sx; \
				coords[1] = coords[1] / params[1].sy; \
				coords[2] = coords[2] / params[1].sz; \
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
			");
	},
	//	TaperX
	//	Definition: inverse mapping
	//		x1<= x <= x2
	//			t = (x-x1)/(x2-x1)
	//			scale = (1-t)*s1 + t*s2
	//			y'=y/scale
	//			z'=z/scale
	//		x < x1   scale = s1
	//		x > x2   scale = s2
	//	Parameters:
	//		x1, x2 - end points of x-interval, x2 > x1
 	//		s1, s2 - scaling factors for end points
	taperX: function(p){
		return CSG.methodFactory(this, p, 
			"	var scale, t;            \
				var s2 = params[1].s2;   \
				var s1 = params[1].s1;   \
				var x2 = params[1].x2;   \
				var x1 = params[1].x1;   \
									     \
				if (coords[0] < x1) {    \
					scale = s1;          \
				} else {                 \
					if(coords[0] > x2) { \
						scale = s2;      \
					} else {             \
						t = (coords[0] - x1) / (x2 - x1); \
						scale = (1-t)*s1 + t*s2; \
					}                    \
				}                        \
				if(Math.abs(scale) < EPS) scale = 1.0; \
				coords[1] = coords[1] / scale; \
				coords[2] = coords[2] / scale; \
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
			");
	},
	//	TaperY
	//	Definition: inverse mapping
	//		y1<= y <= y2
	//			t = (y-y1)/(y2-y1)
	//			scale = (1-t)*s1 + t*s2
	//			z'=z/scale
	//			x'=x/scale
	//		y < y1   scale = s1
	//		y > y2   scale = s2
	//	Parameters:
	//		y1, y2 - end points of y-interval, y2 > y1
 	//		s1, s2 - scaling factors for end points
	taperY: function(p){
		return CSG.methodFactory(this, p, 
			"	var scale, t;            \
				var s2 = params[1].s2;   \
				var s1 = params[1].s1;   \
				var y2 = params[1].y2;   \
				var y1 = params[1].y1;   \
									     \
				if (coords[1] < y1) {    \
					scale = s1;          \
				} else {                 \
					if(coords[1] > y2) { \
						scale = s2;      \
					} else {             \
						t = (coords[1] - y1) / (y2 - y1); \
						scale = (1-t)*s1 + t*s2; \
					}                    \
				}                        \
				if(Math.abs(scale) < EPS) scale = 1.0; \
				coords[2] = coords[2] / scale; \
				coords[0] = coords[0] / scale; \
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
			");
	},
	//	TaperZ
	//	Definition: inverse mapping
	//		z1<= z <= z2
	//			t = (z-z1)/(z2-z1)
	//			scale = (1-t)*s1 + t*s2
	//			x'=x/scale
	//			y'=y/scale
	//		z < z1   scale = s1
	//		z > z2   scale = s2
	//	Parameters:
	//		z1, z2 - end points of z-interval, z2 > z1
 	//		s1, s2 - scaling factors for end points
	taperZ: function(p){
		return CSG.methodFactory(this, p, 
			"	var scale, t;            \
				var s2 = params[1].s2;   \
				var s1 = params[1].s1;   \
				var z2 = params[1].z2;   \
				var z1 = params[1].z1;   \
									     \
				if (coords[2] < z1) {    \
					scale = s1;          \
				} else {                 \
					if(coords[2] > z2) { \
						scale = s2;      \
					} else {             \
						t = (coords[2] - z1) / (z2 - z1); \
						scale = (1-t)*s1 + t*s2; \
					}                    \
				}                        \
				if(Math.abs(scale) < EPS) scale = 1.0; \
				coords[0] = coords[0] / scale; \
				coords[1] = coords[1] / scale; \
				return "+this.func.toString()+"(coords,params[0],attrs[0]);\
			");
	},
};

CSG.methodFactory = function(csg, params, funcDef) {
	return new CSG([csg.params, params], [csg.attrs||{}], new Function('coords', 'params', 'attrs', funcDef));
};

/**Shapes **/

//	Block
//	Definition: x:[vertex[1], vertex[1]+dx], ...
//	Parameters:
//		vertex - block vertex coordinates array
//		dx,dy,dz - edge lengths along x,y,z
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
		result = i0 + z0 - Math.sqrt(i0 * i0 + z0 * z0); 
		return [result, attrs];
	});
};


//	Sphere
//	Definition: R^2-(x-x0)^2-(y-y0)^2-(z-z0)^2
//	Parameters:
//		center - sphere center array
//		R - sphere radius
CSG.sphere = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.radius;
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
		var result = (R * R) - (x * x) - (y * y) - (z * z);
		return [result, attrs];
		});
};

//	Torus with X-axis
//	Parameters:
//		center - center array
//		R - radius of revolution
//		r0 - disk radius	
CSG.torusX = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.R; 
		var r0 = params.r0; 
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
		var result = (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((y * y) + (z * z)); 
		return [result, attrs];
	});
};

//	Torus with Y-axis
//	Parameters:
//		center - center array
//		R - radius of revolution
//		r0 - disk radius
CSG.torusY = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.R; 
		var r0 = params.r0; 
		var x = coords[0] - params.center[0];
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2]; 
		var result =  (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((x * x) + (z * z));
		return [result, attrs];
	});
};

//	Torus with Z-axis
//	Parameters:
//		center - center array
//		R - radius of revolution
//		r0 - disk radius
CSG.torusZ = function(params, attrs) {
	return new CSG(params, attrs, function(coords, params, attrs){
		var R = params.R; 
		var r0 = params.r0; 
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
		var z = coords[2] - params.center[2]; 
		var result =  (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((x * x) + (y * y)); 
		return [result, attrs];
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
        var result =  (Math.cos(x / v) * Math.sin(y / v) + Math.cos(y / v) * Math.sin(z / v)  + Math.cos(z / v) * Math.sin(x / v) + 1.0) - 0.1 * (1 - 0.016 * (r - 10 / r)); 
		return [result, attrs];
	});
};

//	EllipticCylinderX
//	Definition: 1-((y-y0)/a)^2-((z-z0)/b)^2
//	Parameters:
//		center - sphere center array
//		a,b - elliptic half-axes along y,z		
CSG.ellipticCylinderX = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var yt = (coords[1] - params.center[1]) / params.a;
		var zt = (coords[2] - params.center[2]) / params.b;
		var result =  1.0 - yt * yt - zt * zt;
		return [result, attrs];
	});
};

//	EllipticCylinderY
//	Definition: 1-((x-x0)/a)^2-((z-z0)/b)^2
//	Parameters:
//		center - sphere center array
//		a,b - elliptic half-axes along x,z
CSG.ellipticCylinderY = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var xt = (coords[0] - params.center[0]) / params.a;
		var zt = (coords[2] - params.center[2]) / params.b;
		var result =  1.0 - xt * xt - zt * zt;
		return [result, attrs];
	});
};

//	EllipticCylinderZ
//	Definition: 1-((x-x0)/a)^2-((y-y0)/b)^2
//	Parameters:
//		center - sphere center array
//		a,b - elliptic half-axes along x,y
CSG.ellipticCylinderZ = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var xt = (coords[0] - params.center[0]) / params.a;
		var yt = (coords[1] - params.center[1]) / params.b;
		var result =  1.0 - xt * xt - yt * yt;
		return [result, attrs];
	});
};

//	Ellipsoid
//	Definition: 1-((x-x0)/a)^2-((y-y0)/b)^2-((z-z0)/c)^2
//	Parameters:
//		center - sphere center array
//		a,b,c - ellipsoid half-axes along x,y,z
CSG.ellipsoid = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = (coords[0] - params.center[0]) / params.a; 
		var y = (coords[1] - params.center[1]) / params.b; 
		var z = (coords[2] - params.center[2]) / params.c; 
        var result =  1 - (x * x) - (y * y) - (z * z); 
		return [result, attrs];
	});
};

// CylinderX
// Definition: R^2-(y-y0)^2-(z-z0)^2
// Parameters:
//		center - sphere center array
//		R - cylinder radius
CSG.cylinderX = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var y = coords[1] - params.center[1];
		var z = coords[2] - params.center[2]; 
        var result =  (params.R * params.R) - (y * y) - (z * z); 
		return [result, attrs];
	});
};

// CylinderY
// Definition: R^2-(x-x0)^2-(z-z0)^2
// Parameters:
//		center - sphere center array
//		R - cylinder radius
CSG.cylinderY = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = coords[0] - params.center[0]; 
		var z = coords[2] - params.center[2]; 
        var result =  (params.R * params.R) - (x * x) - (z * z); 
		return [result, attrs];
	});
};

// CylinderZ
// Definition: R^2-(x-x0)^2-(y-y0)^2
// Parameters:
//		center - sphere center array
//		R - cylinder radius
CSG.cylinderZ = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = coords[0] - params.center[0]; 
		var y = coords[1] - params.center[1]; 
        var result =  (params.R * params.R) - (x * x) - (y * y); 
		return [result, attrs];
	});
};

CSG.heart = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var x = (coords[0] - params.center[0]); 
		var y = (coords[1] - params.center[1]); 
		var z = (coords[2] - params.center[2]); 
		var pow = Math.pow; 
		var result =  pow(pow(x,2)+(9/4)*pow(y,2)+pow(z,2)-1,3) - pow(x,2)*pow(z,3)-(9/80)*pow(y,2)*pow(z,3); 
		return [result, attrs];
	});
};

// Primitive: Cone with x-axis 
// Definition: (x-x0)^2-((y-y0)/R)^2-((z-z0)/R)^2
// Parameters:
//		center - sphere center array
//		R - radius at height 1 
CSG.coneX = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var xt = (coords[0] - params.center[0]);
		var yt = (coords[1] - params.center[1]) / params.R;
		var zt = (coords[2] - params.center[2]) / params.R;
		var result =  rn (xt*xt) - (yt*yt) - (zt*zt);
		return [result, attrs];
	});
};

// Primitive: Cone with y-axis 
// Definition: (y-y0)^2-((x-x0)/R)^2-((z-z0)/R)^2
// Parameters:
//		center - sphere center array
//		R - radius at height 1 
CSG.coneY = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var xt = (coords[0] - params.center[0]) / params.R;
		var yt = (coords[1] - params.center[1]);
		var zt = (coords[2] - params.center[2]) / params.R;
		var result =  (yt*yt) - (xt*xt) - (zt*zt);
		return [result, attrs];
	});
};

// Primitive: Cone with z-axis 
// Definition: (z-z0)^2-((x-x0)/R)^2-((y-y0)/R)^2
// Parameters:
//		center - sphere center array
//		R - radius at height 1 
CSG.coneZ = function(params, attrs){
	return new CSG(params, attrs, function(coords, params, attrs){
		var xt = (coords[0] - params.center[0]) / params.R;
		var yt = (coords[1] - params.center[1]) / params.R;
		var zt = (coords[2] - params.center[2]);
		var result =  (zt*zt) - (xt*xt) - (yt*yt);
		return [result, attrs];
	});
};
