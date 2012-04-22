
importScripts('frep-csg.js');

var MAX_BOUNDING_BOX = {min:{x:-200.0,y:-200.0,z:-200.0},max:{x:200.0,y:200.0,z:200.0}};
var BOUNDING_BOX_ACCURACY_COARSE = 1;
var BOUNDING_BOX_ACCURACY_FINE = 0.1;
var BOUNDING_BOX_FIXED_PRECISION = 2;

function inc(val, by){
	return parseFloat((val+by).toFixed(BOUNDING_BOX_FIXED_PRECISION));
}
function dec(val, by){
	return parseFloat((val-by).toFixed(BOUNDING_BOX_FIXED_PRECISION));
}

function refineMaxLimit(max, upperLimit, lower1, upper1, lower2, upper2, funcMap){
	var searching = true;
	while (searching) {
		// use full MAX_BOUNDING_BOX as there may be thin spikes outside of the 1/2 region
		for (var i = lower1; i <= upper1; i+=BOUNDING_BOX_ACCURACY_FINE){
			for (var j = lower2; j <= upper2; j+=BOUNDING_BOX_ACCURACY_FINE){
				if (funcMap(i,j,max) >=0){
					max = inc(max, BOUNDING_BOX_ACCURACY_FINE);
				} else {
					searching = false;
				}
			}
		}
		if (max >= upperLimit){
			searching = false;	
		}
	}
	return max;
}

function refineMinLimit(min, lowerLimit, lower1, upper1, lower2, upper2, funcMap){
	var searching = true;
	while (searching) {
		// use full MAX_BOUNDING_BOX as there may be thin spikes outside of the 1/2 region
		for (var i = lower1; i <= upper1; i+=BOUNDING_BOX_ACCURACY_FINE){
			for (var j = lower2; j <= upper2; j+=BOUNDING_BOX_ACCURACY_FINE){
				if (funcMap(i,j,min) >=0){
					min = dec(min, BOUNDING_BOX_ACCURACY_FINE);
				} else {
					searching = false;
				}
			}
		}
		if (min <= lowerLimit){
			searching = false;	
		}
	}
	return min;
}

self.addEventListener('message', function(e) {
	if (!e.data){
		return
	}

	var csgFunction = new Function('', e.data.func)
	var csg = csgFunction();
	var f = csg.func;
	var axis = e.data.axis
	var dir = e.data.dir
	var val = 0;
	var center = [0,0,0];

	if (e.data.params && e.data.params.center) center = e.data.params.center;


	switch (dir){
		case 'min':
			switch (axis){
				case 'x': 
					var minX = center[0]
					while (minX >= MAX_BOUNDING_BOX.min.x && f([minX,center[1],center[2]]) >=0) minX = dec(minX, BOUNDING_BOX_ACCURACY_COARSE);

					val = refineMinLimit(minX, MAX_BOUNDING_BOX.min.x, MAX_BOUNDING_BOX.min.y, MAX_BOUNDING_BOX.max.y, MAX_BOUNDING_BOX.min.z, MAX_BOUNDING_BOX.max.z, function(y,z,x){return f([x,y,z])});
					break;
				case 'y': 
					var minY = center[1]
					while (minY >= MAX_BOUNDING_BOX.min.y && f([center[0],minY,center[2]]) >=0) minY = dec(minY, BOUNDING_BOX_ACCURACY_COARSE);
		
					val = refineMinLimit(minY, MAX_BOUNDING_BOX.min.y, MAX_BOUNDING_BOX.min.x, MAX_BOUNDING_BOX.max.x, MAX_BOUNDING_BOX.min.z, MAX_BOUNDING_BOX.max.z, function(x,z,y){return f([x,y,z])});
					break;
				case 'z':
					var minZ = center[2]
					while (minZ >= MAX_BOUNDING_BOX.min.z && f([center[0],center[1],minZ]) >=0) minZ = dec(minZ, BOUNDING_BOX_ACCURACY_COARSE);
		
					val = refineMinLimit(minZ, MAX_BOUNDING_BOX.min.z, MAX_BOUNDING_BOX.min.x, MAX_BOUNDING_BOX.max.x, MAX_BOUNDING_BOX.min.y, MAX_BOUNDING_BOX.max.y, function(x,y,z){return f([x,y,z])});
					break;
			}
			break;
		case 'max':
			switch (axis){
				case 'x':
					var maxX = center[0]
					while (maxX <= MAX_BOUNDING_BOX.max.x && f([maxX,center[1],center[2]]) >=0) maxX = inc(maxX, BOUNDING_BOX_ACCURACY_COARSE);

					val = refineMaxLimit(maxX, MAX_BOUNDING_BOX.max.x, MAX_BOUNDING_BOX.min.y, MAX_BOUNDING_BOX.max.y, MAX_BOUNDING_BOX.min.z, MAX_BOUNDING_BOX.max.z, function(y,z,x){return f([x,y,z])});
					break;
				case 'y':
					var maxY = center[1]
					while (maxY <= MAX_BOUNDING_BOX.max.y && f([center[0],maxY,center[2]]) >=0) maxY = inc(maxY, BOUNDING_BOX_ACCURACY_COARSE);
 
					val = refineMaxLimit(maxY, MAX_BOUNDING_BOX.max.y, MAX_BOUNDING_BOX.min.x, MAX_BOUNDING_BOX.max.x, MAX_BOUNDING_BOX.min.z, MAX_BOUNDING_BOX.max.z, function(x,z,y){return f([x,y,z])});
					break;
				case 'z':
					var maxZ = center[2]
					while (maxZ <= MAX_BOUNDING_BOX.max.z && f([center[0],center[1],maxZ]) >=0) maxZ = inc(maxZ, BOUNDING_BOX_ACCURACY_COARSE);

					val = refineMaxLimit(maxZ, MAX_BOUNDING_BOX.max.z, MAX_BOUNDING_BOX.min.x, MAX_BOUNDING_BOX.max.x, MAX_BOUNDING_BOX.min.y, MAX_BOUNDING_BOX.max.y, function(x,y,z){return f([x,y,z])});
					break;
			}
			break;
	}

	self.postMessage({'axis':axis,'dir':e.data.dir,'val':val});


	
}, false);