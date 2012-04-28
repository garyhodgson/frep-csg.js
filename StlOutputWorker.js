importScripts('frep-csg.js');

onmessage = function(e){
	if (e.data.vertices == undefined || e.data.indices == undefined){
		return;
	}

	var vertices = e.data.vertices;
	var indices = e.data.indices;

	var X = 0;
	var Y = 1;
	var Z = 2;
	
	var CB = new Array(3);
	var CA = new Array(3);
	var vec_length = 1;

	// Normal Vector for a facet
	var nx, ny, nz;

	var stlOutput = new Array()
	stlOutput.push("solid\n");

    for (var i = 0; i < indices.length/3; i++) {
    	postMessage({'progress':i});
		// Calculate a normal vector from the cross product
		CB[X] = vertices[indices[3*i + 1]][0] - vertices[indices[3*i + 2]][0];
		CB[Y] = vertices[indices[3*i + 1]][1] - vertices[indices[3*i + 2]][1];
		CB[Z] = vertices[indices[3*i + 1]][2] - vertices[indices[3*i + 2]][2];
		CA[X] = vertices[indices[3*i]][0] - vertices[indices[3*i + 2]][0];
		CA[Y] = vertices[indices[3*i]][1] - vertices[indices[3*i + 2]][1];
		CA[Z] = vertices[indices[3*i]][2] - vertices[indices[3*i + 2]][2];
		nx = CB[Y]*CA[Z] - CB[Z]*CA[Y];
		ny = CB[Z]*CA[X] - CB[X]*CA[Z];
		nz = CB[X]*CA[Y] - CB[Y]*CA[X];
		
		// Normalize the calculated normal vector
		vec_length = Math.sqrt(nx*nx + ny*ny + nz*nz);
		nx = nx / vec_length;
		ny = ny / vec_length;
		nz = nz / vec_length;
		
		stlOutput.push(" facet normal " + nx + " " + ny + " " + nz + "\n");
		stlOutput.push("  outer loop" + "\n"); 
		stlOutput.push("   vertex ");
		stlOutput.push(vertices[indices[3*i]][0] + " ");
		stlOutput.push(vertices[indices[3*i]][1] + " ");
		stlOutput.push(vertices[indices[3*i]][2] + "\n");
		stlOutput.push("   vertex ");
		stlOutput.push(vertices[indices[3*i + 1]][0] + " ");
		stlOutput.push(vertices[indices[3*i + 1]][1] + " ");
		stlOutput.push(vertices[indices[3*i + 1]][2] + "\n");
		stlOutput.push("   vertex ");
		stlOutput.push(vertices[indices[3*i + 2]][0] + " ");
		stlOutput.push(vertices[indices[3*i + 2]][1] + " ");
		stlOutput.push(vertices[indices[3*i + 2]][2] + "\n");
		stlOutput.push("  endloop" + "\n"); 
		stlOutput.push(" endfacet" + "\n"); 
    }
    
    stlOutput.push("endsolid\n");

    function errorHandler(e) {
	    var msg = '';
	    switch (e.code) {
	      case FileError.QUOTA_EXCEEDED_ERR:
	        msg = 'QUOTA_EXCEEDED_ERR';
	        break;
	      case FileError.NOT_FOUND_ERR:
	        msg = 'NOT_FOUND_ERR';
	        break;
	      case FileError.SECURITY_ERR:
	        msg = 'SECURITY_ERR';
	        break;
	      case FileError.INVALID_MODIFICATION_ERR:
	        msg = 'INVALID_MODIFICATION_ERR';
	        break;
	      case FileError.INVALID_STATE_ERR:
	        msg = 'INVALID_STATE_ERR';
	        break;
	      default:
	        msg = 'Unknown Error';
	        break;
	    };
	    postMessage({'msg':'STL write error. ' + msg});
  }

	webkitRequestFileSystem(TEMPORARY, 0, function(fs) {
	    fs.root.getFile('out.stl', {create: true}, function(fileEntry) {
	        fileEntry.createWriter(function(fileWriter) {

	        	
				fileWriter.onwriteend = function(e) {
					fileWriter.onwriteend = null; // Avoid an infinite loop.
					fileEntry.createWriter(function(fw) {

						var builder = new WebKitBlobBuilder();

		                for (var i = 0; i < stlOutput.length; i++){
		                	builder.append(stlOutput[i]);
		                }

						fw.onwriteend = function(e) {
							postMessage({'msg':'STL write completed.', 'url':fileEntry.toURL()});
						};

	  					fw.write(builder.getBlob('text/plain'));						
					});
	            };

				fileWriter.onerror = function(e) {
					postMessage({'msg':'STL write failed: ' + e.toString()});
				};
				fileWriter.truncate(0);					
	            
	        }, function(error) {errorHandler});
	    }, function(error) {errorHandler});
	}, function(error) {errorHandler});

}

