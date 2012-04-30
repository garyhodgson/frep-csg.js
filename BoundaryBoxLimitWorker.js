
importScripts('frep-csg.js', 'js/underscore-min.js');

var BOUNDING_BOX_ACCURACY = 0.5;
var BOUNDING_BOX_FIXED_PRECISION = 2;

/*
Algorithm: 
	check min
	if !found return min
	if found check max
	if found return max

	if !found check middle

	if found
		if distance since last check <= accuracy defaultValue?
			yes -> return last value
			no -> repeat in direction away from center
	if !found 
		if distance since last check <= accuracy defaultValue?
			yes -> return value
			no -> repeat in direction towards center
*/

onmessage = function(e) {
	if (!e.data){
		return
	}
	var f = eval("("+e.data.funcDef+")");
	localCsg = new CSG(e.data.params, e.data.attrs, f);
	
	var axis = e.data.axis
	var dir = e.data.dir
	var val = 0;
	var center = [0,0,0];

	if (e.data.params && e.data.params.center) center = e.data.params.center;

	switch (dir){
		case 'min':
			switch (axis){
				case 'x':
					val = planeSearch('y', 'z', 'x', dir);
					break;
				case 'y': 
					val = planeSearch('x', 'z', 'y', dir);
					break;
				case 'z':
					val = planeSearch('x', 'y', 'z', dir);
					break;
			}
			break;
		case 'max':
			switch (axis){
				case 'x':
					val = planeSearch('y', 'z', 'x', dir);
					break;
				case 'y':
					val = planeSearch('x', 'z', 'y', dir);
					break;
				case 'z':
					val = planeSearch('x', 'y', 'z', dir);
					break;
			}
			break;
	}
	self.postMessage({'axis':axis,'dir':dir,'val':val});

}

function samplePlane(k, func, plane, step){

	//shortcut - check middle point
	var val = func(plane.i.max/2,plane.j.max/2,k)
	if (val >= 0){
		//console.log(""+k+"-> hit (*) ");
		return true;
	}

	for (var i = plane.i.min; i <= plane.i.max; i+=step){
		for (var j = plane.j.min; j <= plane.j.max; j+=step){
			var val2 = func(i,j,k)
			if (val2 >= 0){
				//console.log(""+k+"-> hit ");
				return true;
			}
		}
	}
	//console.log(""+k+"-> miss ");
	return false
}

function recursivePlaneSearch(targetPoint, f, plane, delta, lastMiss, lastHit){
	//console.log("tp:",targetPoint, "d:", delta, "lM:", lastMiss, "lH:", lastHit);
	
	found = samplePlane(targetPoint, f, plane, BOUNDING_BOX_ACCURACY);

	if (found){
		if (Math.abs(delta) <= BOUNDING_BOX_ACCURACY){
			return lastMiss.toFixed(BOUNDING_BOX_FIXED_PRECISION);
		} else {
			var newTargetPointAwayFromCenter = (lastMiss + targetPoint)/2;
			var newDelta = targetPoint - newTargetPointAwayFromCenter;
			return recursivePlaneSearch(newTargetPointAwayFromCenter, f, plane, newDelta, lastMiss, targetPoint);
		}

	} else {
		if (Math.abs(delta) <= BOUNDING_BOX_ACCURACY){
			return targetPoint.toFixed(BOUNDING_BOX_FIXED_PRECISION);
		} else {
			var newTargetPointTowardsCenter = (lastHit + targetPoint)/2;
			var newDelta = targetPoint + newTargetPointTowardsCenter;
			return recursivePlaneSearch(newTargetPointTowardsCenter, f, plane, newDelta, targetPoint, lastHit);
		}
	}

}

function planeSearch(i, j, k, dir){
	
	var plane = {i:{min:MIN_BOUNDING_BOX.min[i],max:MIN_BOUNDING_BOX.max[i]}, j:{min:MIN_BOUNDING_BOX.min[j],max:MIN_BOUNDING_BOX.max[j]}};

	var f = new Function(i, j, k, "return localCsg.call([x,y,z])[0]")

	var targetPoint = MIN_BOUNDING_BOX[dir][k];

	var found = samplePlane(targetPoint, f, plane, BOUNDING_BOX_ACCURACY);

	if (!found) return targetPoint; // model within min boundary

	targetPoint = MAX_BOUNDING_BOX[dir][k];

	found = samplePlane(targetPoint, f, plane, BOUNDING_BOX_ACCURACY);

	if (found) return targetPoint; // model touches max boundary

	var lastTargetPoint = targetPoint;

	targetPoint = lastTargetPoint/2;

	var delta = lastTargetPoint - targetPoint;

	return recursivePlaneSearch(targetPoint, f, plane, delta, lastTargetPoint, null);
}

