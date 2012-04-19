
var ALPHA = 0.0;

CSG = function(p, a, f) {
	this.func = f;
	this.params = p;
	this.attrs = a;
	this.vertices = new Array();
	this.normals = new Array();
	this.indices = new Array();
}


CSG.prototype = {

	toMesh: function(g,bb,iso) {

		var grid = g?g:50;
		var boundingBox = bb?bb:{min:[-1,-1,-1],max:[5,5,5]};
		var isosurface = iso?iso:0.0;

		var mesh = new GL.Mesh({ normals: true, colors: true });
		var indexer = new GL.Indexer();

		this.vertices = new Array();
		this.normals = new Array();
		this.indices = new Array();

		var polygonizer = new CSG.Polygonizer(
				boundingBox.min, 
				boundingBox.max,
				[grid,grid,grid],
				isosurface,
				this.func,
				this.params,
				this.attrs);

		polygonizer.polygonize(this.vertices, this.normals, this.indices);

		notify("Vertices: "+this.vertices.length+ "; Normals: "+ this.normals.length+ "; Indices: "+ this.indices.length+";");

		for (var i = 2; i < this.indices.length; i+=3) {
			mesh.triangles.push([this.indices[i-2], this.indices[i - 1], this.indices[i]]);
		}

		mesh.vertices = this.vertices;
		mesh.normals = this.normals;
		mesh.colors = this.vertices.map(function(v){return [0.85,0.85,0.85,0.85];});
		mesh.computeWireframe();

		return mesh;
	},
	toStl: function(){
		if (this.vertices == undefined || this.vertices.length == 0){
			notify("No vertices found, Polygonise!")
			return;
		}

		var X = 0;
		var Y = 1;
		var Z = 2;
		
		var CB = new Array(3);
		var CA = new Array(3);
		var vec_length = 1;

		// Normal Vector for a facet
		var nx, ny, nz;

		Point3D = function(x,y,z){
			this.x = x;
			this.y = y;
			this.z = z;
		}

		notify("solid\n");

	    for (var i = 0; i < this.indices.length/3; i++) {

console.log("stl: "+i)

			// Calculate a normal vector from the cross product
			CB[X] = this.vertices[this.indices[3*i + 1]][0] - this.vertices[this.indices[3*i + 2]][0];
			CB[Y] = this.vertices[this.indices[3*i + 1]][1] - this.vertices[this.indices[3*i + 2]][1];
			CB[Z] = this.vertices[this.indices[3*i + 1]][2] - this.vertices[this.indices[3*i + 2]][2];
			CA[X] = this.vertices[this.indices[3*i]][0] - this.vertices[this.indices[3*i + 2]][0];
			CA[Y] = this.vertices[this.indices[3*i]][1] - this.vertices[this.indices[3*i + 2]][1];
			CA[Z] = this.vertices[this.indices[3*i]][2] - this.vertices[this.indices[3*i + 2]][2];
			nx = CB[Y]*CA[Z] - CB[Z]*CA[Y];
			ny = CB[Z]*CA[X] - CB[X]*CA[Z];
			nz = CB[X]*CA[Y] - CB[Y]*CA[X];
			
			// Normalize the calculated normal vector
			vec_length = Math.sqrt(nx*nx + ny*ny + nz*nz);
			nx = nx / vec_length;
			ny = ny / vec_length;
			nz = nz / vec_length;
			
			notify(" facet normal " + nx + " " + ny + " " + nz + "\n");
			notify("  outer loop" + "\n"); 
			notify("   vertex ");
			notify(this.vertices[this.indices[3*i]][0] + " ");
			notify(this.vertices[this.indices[3*i]][1] + " ");
			notify(this.vertices[this.indices[3*i]][2] + "\n");
			notify("   vertex ");
			notify(this.vertices[this.indices[3*i + 1]][0] + " ");
			notify(this.vertices[this.indices[3*i + 1]][1] + " ");
			notify(this.vertices[this.indices[3*i + 1]][2] + "\n");
			notify("   vertex ");
			notify(this.vertices[this.indices[3*i + 2]][0] + " ");
			notify(this.vertices[this.indices[3*i + 2]][1] + " ");
			notify(this.vertices[this.indices[3*i + 2]][2] + "\n");
			notify("  endloop" + "\n"); 
			notify(" endfacet" + "\n"); 
	    }
	    
	    notify("endsolid\n");




	},
	union: function(otherCsg){
		var csg = new CSG();
		var that = this;
		csg.func = function(coords) {
			var result1 = that.func(coords, that.params);
			var result2 = otherCsg.func(coords, otherCsg.params);
			//return (1/(1+ALPHA))*(result2 + result1 + Math.sqrt(Math.abs(result2*result2 + result1*result1-2*ALPHA*result1*result2)));
			return Math.max(result2, result1);
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


/** Taken from Hyperfun Applet Polygonizer by Yuichiro Goto **/
CSG.Polygonizer = function(min,max,div,isovalue,func,params,attrs) {
	// Lower left front corner of the bounding box
	this.xMin = min[0];
	this.yMin = min[1];
	this.zMin = min[2];

	// Upper right back corner of the bounding box
	this.xMax = max[0];
	this.yMax = max[1];
	this.zMax = max[2];

	// Number of divisions along each axis of the bounding box
	this.xDiv = div[0];
	this.yDiv = div[1];
	this.zDiv = div[2];
	

	// Isovalue of the isosurface
	this.isovalue = isovalue;
	// Function defining the isosurface
	this.func = func;
	this.params = params;
	this.attrs = attrs;

	var dx, dy, dz;
	var ndx, ndy, ndz;
}


CSG.Polygonizer.prototype = {
	lerp: function(t, v0, v1) {
		return [ 	v0[0] + t * (v1[0] - v0[0]),
					v0[1] + t * (v1[1] - v0[1]), 
					v0[2] + t * (v1[2] - v0[2]) 
				];
	},
	calcNormal: function(v) {
		var x = v[0];
		var y = v[1];
		var z = v[2];

		var f = this.func([x, y, z]);
		var nx = -(this.func([x + this.ndx, y, z]) - f) / this.ndx;
		var ny = -(this.func([x, y + this.ndy, z]) - f) / this.ndy;
		var nz = -(this.func([x, y, z + this.ndz]) - f) / this.ndz;

		var len = Math.sqrt(nx * nx + ny * ny + nz * nz);
		if (len > 0.0) {
			nx /= len;
			ny /= len;
			nz /= len;
		}
		return [ nx, ny, nz ];
	},
	sample: function(plane, z) {
		for (var j = 0; j <= this.yDiv; j++) {
			var y = this.yMin + j * this.dy;
			for (var i = 0; i <= this.xDiv; i++) {
				var x = this.xMin + i * this.dx;
				plane[j][i] = this.func([x,y,z]);
			}
		}
	},

	polygonize: function(vertices, normals, indices) {

		var eps = (this.isovalue == 0.0) ? 1.0E-5 : this.isovalue * 1.0E-5;

		var values = new Array(8);
		for (var i = 0; i < values.length; i++) {
			values[i] = 0.0;
		}
		
		var positionsD = new Array(8)
		var positionsI = new Array(8)
		for (var i = 0; i < positionsD.length; i++) {
			positionsD[i] = new Array(3);
			positionsI[i] = new Array(3);
		}

		var upperPlane = new Array(this.yDiv + 1);
		var lowerPlane = new Array(this.yDiv + 1);
		for (var i = 0; i < upperPlane.length; i++) {
			upperPlane[i] = new Array(this.xDiv + 1);
			lowerPlane[i] = new Array(this.xDiv + 1);
		}

		var connectionSwitches = new Array(6);
		for (var i = 0; i < connectionSwitches.length; i++) {
			connectionSwitches[i] = 0;
		}
		var edgeToIndex = new Array(12);
		for (var i = 0; i < edgeToIndex.length; i++) {
			edgeToIndex[i] = 0;
		}
		this.indexTable = {};

		this.dx = (this.xMax - this.xMin) / this.xDiv;
		this.dy = (this.yMax - this.yMin) / this.yDiv;
		this.dz = (this.zMax - this.zMin) / this.zDiv;

		this.ndx = 0.001 * this.dx;
		this.ndy = 0.001 * this.dy;
		this.ndz = 0.001 * this.dz;

		this.sample(lowerPlane, this.zMin);

		for (var k = 0; k < this.zDiv; k++) {
			var zLower = this.zMin + k * this.dz;
			var zUpper = this.zMin + (k+1) * this.dz;
			
			this.sample(upperPlane, zUpper);

			console.log("z: "+k);

			for (var j = 0; j < this.yDiv; j++) {
				var yLower = this.yMin + j * this.dy;
				var yUpper = this.yMin + (j+1) * this.dy;

				for (var i = 0; i < this.xDiv; i++) {
					var xLower = this.xMin + i * this.dx;
					var xUpper = this.xMin + (i+1) * this.dx;

					// Set sampled function values on each corner of the cube
					values[0] = lowerPlane[j][i];
					values[1] = lowerPlane[j + 1][i];
					values[2] = lowerPlane[j + 1][i + 1];
					values[3] = lowerPlane[j][i + 1];
					values[4] = upperPlane[j][i];
					values[5] = upperPlane[j + 1][i];
					values[6] = upperPlane[j + 1][i + 1];
					values[7] = upperPlane[j][i + 1];

					// Adjust the function values which are almost same as the isovalue
					for (var v in values){
						if (Math.abs(values[v] - this.isovalue) < eps) {
							values[v] += 10.0 * eps;
						}
					}

					// Calculate index into the lookup table
					var cubeIndex = 0;
					if (values[0] > this.isovalue) {
						cubeIndex += 1;
					}
					if (values[1] > this.isovalue) {
						cubeIndex += 2;
					}
					if (values[2] > this.isovalue) {
						cubeIndex += 4;
					}
					if (values[3] > this.isovalue) {
						cubeIndex += 8;
					}
					if (values[4] > this.isovalue) {
						cubeIndex += 16;
					}
					if (values[5] > this.isovalue) {
						cubeIndex += 32;
					}
					if (values[6] > this.isovalue) {
						cubeIndex += 64;
					}
					if (values[7] > this.isovalue) {
						cubeIndex += 128;
					}

					// Skip the empty cube
					if (cubeIndex == 0 || cubeIndex == 255) {
						//console.log("Skip the empty cube");
						continue;
					}
					
					var cube = CSG.LookupTable.getCube(cubeIndex);
					// Set up corner positions of the cube
					positionsD[0][0] = xLower;
					positionsD[0][1] = yLower;
					positionsD[0][2] = zLower;
					positionsD[1][0] = xLower;
					positionsD[1][1] = yUpper;
					positionsD[1][2] = zLower;
					positionsD[2][0] = xUpper;
					positionsD[2][1] = yUpper;
					positionsD[2][2] = zLower;
					positionsD[3][0] = xUpper;
					positionsD[3][1] = yLower;
					positionsD[3][2] = zLower;
					positionsD[4][0] = xLower;
					positionsD[4][1] = yLower;
					positionsD[4][2] = zUpper;
					positionsD[5][0] = xLower;
					positionsD[5][1] = yUpper;
					positionsD[5][2] = zUpper;
					positionsD[6][0] = xUpper;
					positionsD[6][1] = yUpper;
					positionsD[6][2] = zUpper;
					positionsD[7][0] = xUpper;
					positionsD[7][1] = yLower;
					positionsD[7][2] = zUpper;

					positionsI[0][0] = i;
					positionsI[0][1] = j;
					positionsI[0][2] = k;
					positionsI[1][0] = i;
					positionsI[1][1] = j + 1;
					positionsI[1][2] = k;
					positionsI[2][0] = i + 1;
					positionsI[2][1] = j + 1;
					positionsI[2][2] = k;
					positionsI[3][0] = i + 1;
					positionsI[3][1] = j;
					positionsI[3][2] = k;
					positionsI[4][0] = i;
					positionsI[4][1] = j;
					positionsI[4][2] = k + 1;
					positionsI[5][0] = i;
					positionsI[5][1] = j + 1;
					positionsI[5][2] = k + 1;
					positionsI[6][0] = i + 1;
					positionsI[6][1] = j + 1;
					positionsI[6][2] = k + 1;
					positionsI[7][0] = i + 1;
					positionsI[7][1] = j;
					positionsI[7][2] = k + 1;

					// Find the cube edges which have intersection points with the isosurface
					for (var edgeIndex = 0; edgeIndex < 12; edgeIndex++) {
						var edge = cube.getEdge(edgeIndex);

						if (edge.getConnectedEdge(0) !== undefined) {
							var key = new CSG.EdgeKey(positionsI[edge.getStartVertexIndex()],positionsI[edge.getEndVertexIndex()]);
							if (this.indexTable.hasOwnProperty(key)) {
								edgeToIndex[edgeIndex] = this.indexTable[key];
							} else {
								var t = (this.isovalue - values[edge.getStartVertexIndex()]) / (values[edge.getEndVertexIndex()] - values[edge.getStartVertexIndex()]);
								var v = this.lerp(t, positionsD[edge.getStartVertexIndex()], positionsD[edge.getEndVertexIndex()]);
								vertices.push(v);
								if (normals !== undefined) {
									normals.push(this.calcNormal(v));
								}
								this.indexTable[key] = (edgeToIndex[edgeIndex] = vertices.length - 1);
							}
						}
					}

					// Resolve topological ambiguity on cube faces
					for (var faceIndex = 0; faceIndex < 6; faceIndex++) {
						var face = cube.getFace(faceIndex);
						if (face.isAmbiguous()) {
							var d0 = values[face.getEdge(0).getEndVertexIndex()] - values[face.getEdge(0).getStartVertexIndex()];
							var d1 = values[face.getEdge(2).getEndVertexIndex()] - values[face.getEdge(2).getStartVertexIndex()];
							var t = (this.isovalue - values[face.getEdge(1).getStartVertexIndex()]) / (values[face.getEdge(1).getEndVertexIndex()] - values[face.getEdge(1).getStartVertexIndex()]);
							connectionSwitches[faceIndex] = (t > -d0 / (d1 - d0)) ? 1 : 0;
						} else {
							connectionSwitches[faceIndex] = 0;
						}
					}

					// Get the connectivity graph of the cube edges and trace it to generate triangles
					var connectivity = cube.getEdgeConnectivity(connectionSwitches);

					for (var edgeIndex = 0; edgeIndex < 12;) {
						if (connectivity[edgeIndex] != -1) {
							var index0 = edgeIndex;
							var index1 = connectivity[index0];
							var index2 = connectivity[index1];

							indices.push(edgeToIndex[index0]);
							indices.push(edgeToIndex[index1]);
							indices.push(edgeToIndex[index2]);

							connectivity[index0] = -1;
							connectivity[index1] = -1;
							if (connectivity[index2] != index0) {
								connectivity[index0] = index2;
								continue;
							}
							connectivity[index2] = -1;
						}
						edgeIndex++;
					}

				};
			};
			var tmp = lowerPlane;
			lowerPlane = upperPlane;
			upperPlane = tmp;
		};
	}
};

CSG.Edge = function(index){
	var EDGE_VERTICES = [[0, 1], [1, 2], [3, 2], [0, 3],
						 [4, 5], [5, 6], [7, 6], [4, 7],
						 [0, 4], [1, 5], [2, 6], [3, 7]];
	var index = index;
    var startVertexIndex = EDGE_VERTICES[index][0];
    var endVertexIndex = EDGE_VERTICES[index][1];
    var connectedEdge0 = undefined;
    var connectedEdge1 = undefined;

    return {
		getIndex: function(){
			return index;		
		},
		getStartVertexIndex: function() {
			return startVertexIndex;
		},
		getEndVertexIndex: function() {
			return endVertexIndex;
		},
		setConnectedEdge: function(index, edge){
			if (index != 0 && index != 1) {
			    console.error("Edge.setConnectedEdge: IndexOutOfBoundsException!");
			}
			if (index == 0) {
			    connectedEdge0 = edge;
			} else {
			    connectedEdge1 = edge;
			}
		},
		getConnectedEdge: function(index) {
			if (index != 0 && index != 1) {
			    console.error("Edge.getConnectedEdge: IndexOutOfBoundsException!");
			}
			return (index == 0) ? connectedEdge0 : connectedEdge1;
		},
		toString: function(){
			return "Edge" + index + "[" + startVertexIndex + "," + endVertexIndex + "]";
		}
    }
};


CSG.EdgeKey = function(p0, p1){
	var BIT_SHIFT = 10;
	var BIT_MASK = (1 << BIT_SHIFT) - 1;

	var i0 = p0[0]; 
	var j0 = p0[1];
	var k0 = p0[2];
	var i1 = p1[0]; 
	var j1 = p1[1];
	var k1 = p1[2];

	if (i0 < i1 || (i0 == i1 && (j0 < j1 || (j0 == j1 && k0 < k1)))) {
		// do nothing....
	} else {
		i0 = i1;
		j0 = j1;
		k0 = k1;
		i1 = i0;
		j1 = j0;
		k1 = k0;
	}

	return {
		equals: function(obj) {
			if (this == obj) {
				return true;
			}
			if (obj instanceof EdgeKey) {
				var key = obj;
				if (i0 == key.i0 && j0 == key.j0 && k0 == key.k0
						&& i1 == key.i1 && j1 == key.j1 && k1 == key.k1) {
					return true;
				}
			}
			return false;
		},
		hashCode: function() {
			return (((((i0 & BIT_MASK) << BIT_SHIFT) | (j0 & BIT_MASK)) << BIT_SHIFT) | (k0 & BIT_MASK))
					+ (((((i1 & BIT_MASK) << BIT_SHIFT) | (j1 & BIT_MASK)) << BIT_SHIFT) | (k1 & BIT_MASK));
		},
		toString: function(){
			return "EdgeKey ["+i0+","+j0+","+k0+", "+i1+","+j1+","+k1+"]"
		}
	}
};

CSG.Face = function(index, edges, ambiguous){

	var index = index;
    var edges = edges;
    var ambiguous = ambiguous;

	return {
		getIndex: function() {
			return index;
		},
		getEdge: function(index) {
			return edges[index];
		},
		getEdgeCount: function() {
			return edges.length;
		},
		isAmbiguous: function() {
			return ambiguous;
		},
		contains: function(edge) {
			return (edge == edges[0] || edge == edges[1] ||
			edge == edges[2] || edge == edges[3]) ? true : false;
		},
		toString: function() {
			return "Face" + index + "[" + edges[0] + "," + edges[1] + "," +
		      edges[2] + "," + edges[3] + "]" + (ambiguous ? "*" : "");
		}
	}
};

CSG.FaceFactory = function(){
	var FACE_VERTICES = [
			[0, 1, 2, 3],
			[0, 1, 5, 4],
			[0, 3, 7, 4],
			[4, 5, 6, 7],
			[3, 2, 6, 7],
			[1, 2, 6, 5]
		];

	var FACE_EDGES = [  
			[0,  1,  2,  3],
			[0,  9,  4,  8],
			[3, 11,  7,  8],
			[4,  5,  6,  7],
			[2, 10,  6, 11],
			[1, 10,  5,  9]
		];

	var EDGE_CONNECTIVITY_ON_FACE = [
			[[-1,-1,-1,-1], undefined],
			[[-1,-1,-1, 0], undefined],
			[[ 1,-1,-1,-1], undefined],
			[[-1,-1,-1, 1], undefined],
			[[-1, 2,-1,-1], undefined],
			[[-1, 0,-1, 2], [-1, 2,-1, 0]],
			[[ 2,-1,-1,-1], undefined],
			[[-1,-1,-1, 2], undefined],
			[[-1,-1, 3,-1], undefined],
			[[-1,-1, 0,-1], undefined],
			[[ 1,-1, 3,-1], [ 3,-1, 1,-1]],
			[[-1,-1, 1,-1], undefined],
			[[-1, 3,-1,-1], undefined],
			[[-1, 0,-1,-1], undefined],
			[[ 3,-1,-1,-1], undefined],
			[[-1,-1,-1,-1], undefined] 
		];

    var CW  = 1;
    var CCW = 0;
	var FACE_ORIENTATION = [CW, CCW, CW, CCW, CW, CCW];

	function isAmbiguousBitPattern(bitPatternOnFace) {
		return (bitPatternOnFace == 5 || bitPatternOnFace == 10) ? true : false;
	};

	function isBitOn(bitPatternOnCube, vertexIndex) {
		return ((bitPatternOnCube & (1 << vertexIndex)) != 0) ? true : false;
    };

    function buildBitPatternOnFace(bitPatternOnCube, faceIndex) {
		var bitPatternOnFace = 0;
		for (var vertexIndex = 0; vertexIndex < 4; vertexIndex++) {
		    if (isBitOn(bitPatternOnCube, FACE_VERTICES[faceIndex][vertexIndex])) {
				bitPatternOnFace |= 1 << vertexIndex;
		    }
		}
		return bitPatternOnFace;
    };

	return {
		createFace: function(faceIndex, bitPatternOnCube, edges) {
			if (faceIndex < 0 || faceIndex > 5) {
				console.error("IllegalArgumentException - faceIndex must be in the range between 0 and 5");
				return;
			}
			if (bitPatternOnCube < 0 || bitPatternOnCube > 255) {
				console.error("IllegalArgumentException - bitPatternOnCube must be in the range between 0 and 255");
				return;
			}
			if (edges.length != 12) {
				console.error("IllegalArgumentException - length of edges must be 12");
				return;
			}
			var bitPatternOnFace = buildBitPatternOnFace(bitPatternOnCube, faceIndex);


			

			var face = new CSG.Face(faceIndex, [edges[FACE_EDGES[faceIndex][0]],
												 edges[FACE_EDGES[faceIndex][1]],
												 edges[FACE_EDGES[faceIndex][2]],
												 edges[FACE_EDGES[faceIndex][3]]], 
												 isAmbiguousBitPattern(bitPatternOnFace));

			var connectivity = EDGE_CONNECTIVITY_ON_FACE[bitPatternOnFace];
			for (var i = 0; i < 2; i++) {
				if (connectivity[i] !== undefined) {
					for (var vertexIndex = 0; vertexIndex < 4; vertexIndex++) {
						if (connectivity[i][vertexIndex] != -1) {
							if (FACE_ORIENTATION[faceIndex] == CW) {
								var edge = face.getEdge(vertexIndex);
								edge.setConnectedEdge(i, face.getEdge(connectivity[i][vertexIndex]));
							} else {
								var edge = face.getEdge(connectivity[i][vertexIndex]);
								edge.setConnectedEdge(i, face.getEdge(vertexIndex));
							}
						}
					}
				}
			}
			return face;
		}
	}
}();

CSG.Cube = function(index){

	var index = index;
	var edges = new Array(12);
	for (var i = 0; i < 12; i++) {
		edges[i] = new CSG.Edge(i);
	};
	var faces = new Array(6);
	for (var i = 0; i < 6; i++) {
		faces[i] = CSG.FaceFactory.createFace(i, index, edges);
	};

	return {
		getIndex: function() {
			return index;
		},
		getEdge: function(index) {
			return edges[index];
		},
		getEdgeCount: function() {
			return edges.length;
		},
		getFace: function(index) {
			return faces[index];
    	},
		getFaceCount: function() {
			return faces.length;
    	},
    	indexToString: function(index) {
			return  (((index & (1<<7)) != 0) ? "1" : "0") +
					(((index & (1<<6)) != 0) ? "1" : "0") +
					(((index & (1<<5)) != 0) ? "1" : "0") +
					(((index & (1<<4)) != 0) ? "1" : "0") +
					(((index & (1<<3)) != 0) ? "1" : "0") +
					(((index & (1<<2)) != 0) ? "1" : "0") +
					(((index & (1<<1)) != 0) ? "1" : "0") +
					(((index & (1<<0)) != 0) ? "1" : "0");
		},
		getEdgeConnectivity: function(connectionSwitches) {
			var connectivity = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
			for (var faceIndex = 0; faceIndex < 6; faceIndex++) {
				var face = faces[faceIndex];
				if (face.isAmbiguous() == false && connectionSwitches[faceIndex] != 0) {
					console.error("Face.getEdgeConnectivity FaceNotAmbiguousException");
				}
				for (var edgeIndex = 0; edgeIndex < 4; edgeIndex++) {
					var edge = face.getEdge(edgeIndex);
					if (edge.getConnectedEdge(0) !== undefined && face.contains(edge.getConnectedEdge(0))) {
						connectivity[edge.getIndex()] = edge.getConnectedEdge(connectionSwitches[faceIndex]).getIndex();
					}
				}
			}
			return connectivity;
	    },
		toString: function(){
			return "Cube" + index + "[" + faces[0] + "," + faces[1] + "," + faces[2] + "," +
				      faces[3] + "," + faces[4] + "," + faces[5] + "]";
		}
	}
};

CSG.LookupTable = function(){
	var cubes = new Array(256);
	for (var i = 0; i < 256; i++) {
		cubes[i] = new CSG.Cube(i);
	};

	return {
		getCube: function(cubeIndex){
			return cubes[cubeIndex];
		},
		getCubeCount: function(){
			return cubes.length;
		}
	}
}();


/**Shapes **/

CSG.block = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var vertex = this.params.vertex;
		var dx = this.params.dx;
		var dy = this.params.dy;
		var dz = this.params.dz;

		var x0 = -(coords[0] - vertex[0]) * (coords[0] - (vertex[0] + dx));
		var y0 = -(coords[1] - vertex[1]) * (coords[1] - (vertex[1] + dy));
		var z0 = -(coords[2] - vertex[2]) * (coords[2] - (vertex[2] + dz));
		var i0 = x0 + y0 - Math.sqrt(x0 * x0 + y0 * y0);
		return i0 + z0 - Math.sqrt(i0 * i0 + z0 * z0);
	});
};


CSG.sphere = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var R = this.params.radius;
		var x = coords[0] - this.params.center[0];
		var y = coords[1] - this.params.center[1];
		var z = coords[2] - this.params.center[2];
		return (R * R) - (x * x) - (y * y) - (z * z);
	});
};

CSG.torusX = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var R = this.params.R;
		var r0 = this.params.r0;

		var x = coords[0] - this.params.center[0];
		var y = coords[1] - this.params.center[1];
		var z = coords[2] - this.params.center[2];

		return (r0 * r0) - (x * x) - (y * y) - (z * z) - (R * R) + 2 * R * Math.sqrt((y * y) + (z * z));
	});
};

CSG.gyroid = function(params, attrs) {
	return new CSG(params, attrs, function(coords) {
		var x = coords[0] - this.params.center[0];
		var y = coords[1] - this.params.center[1];
		var z = coords[2] - this.params.center[2];
        var r = x * x + y * y + z * z;
        var	ti = Math.abs(Math.sin(this.params.t / 100000)) / 10 + 0.06;
        var v = ti * r;
        return (Math.cos(x / v) * Math.sin(y / v) + Math.cos(y / v) * Math.sin(z / v) 
        	+ Math.cos(z / v) * Math.sin(x / v) + 1.0) - 0.1 * (1 - 0.016 * (r - 10 / r));
	});
};

CSG.ellipsoid = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var x = (coords[0] - this.params.center[0]) / this.params.a;
		var y = (coords[1] - this.params.center[1]) / this.params.b;
		var z = (coords[2] - this.params.center[2]) / this.params.c;
        
        return 1 - (x * x) - (y * y) - (z * z);
	});
};

CSG.cylinderY = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var x = coords[0] - this.params.center[0];
		var z = coords[2] - this.params.center[2];
        return (this.params.R * this.params.R) - (x * x) - (z * z);
	});
};

CSG.cylinderX = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var y = coords[1] - this.params.center[1];
		var z = coords[2] - this.params.center[2];
        return (this.params.R * this.params.R) - (y * y) - (z * z);
	});
};

CSG.cylinderZ = function(params, attrs){
	return new CSG(params, attrs, function(coords) {
		var x = coords[0] - this.params.center[0];
		var y = coords[1] - this.params.center[1];
        return (this.params.R * this.params.R) - (x * x) - (y * y);
	});
};