
var angleX = 291;
var angleY = 0;
var angleZ = 337;
var viewer;
var depth = 20;
var xPos = -3;
var yPos = -2;

// A viewer is a WebGL canvas that lets the user view a mesh. The user can
// tumble it around by dragging the mouse.
function Viewer(width, height, depth, id) {
  viewer = this;
  depth = depth;
  gl = GL.create();
  this.gl = gl;

  // Set up the viewport
  gl.canvas.width = width;
  gl.canvas.height = height;
  gl.viewport(0, 0, width, height);
  gl.matrixMode(gl.PROJECTION);
  gl.loadIdentity();
  gl.perspective(45, width / height, 0.5, 1000);
  gl.matrixMode(gl.MODELVIEW);

  // Set up WebGL state
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.clearColor(0.93, 0.93, 0.93, 1);
  gl.enable(gl.DEPTH_TEST);
  gl.enable(gl.CULL_FACE);
  gl.polygonOffset(1, 1);

  var that = this;

  // Black shader for wireframe
  this.blackShader = new GL.Shader('\
    void main() {\
      gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\
    }\
  ', '\
    void main() {\
      gl_FragColor = vec4(0.0, 0.0, 0.0, 0.1);\
    }\
  ');

  // Shader with diffuse and specular lighting
  this.lightingShader = new GL.Shader('\
    varying vec3 color;\
    varying vec3 normal;\
    varying vec3 light;\
    void main() {\
      const vec3 lightDir = vec3(1.0, 2.0, 3.0) / 3.741657386773941;\
      light = (gl_ModelViewMatrix * vec4(lightDir, 0.0)).xyz;\
      color = gl_Color.rgb;\
      normal = gl_NormalMatrix * gl_Normal;\
      gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\
    }\
  ', '\
    varying vec3 color;\
    varying vec3 normal;\
    varying vec3 light;\
    void main() {\
      vec3 n = normalize(normal);\
      float diffuse = max(0.0, dot(light, n));\
      float specular = pow(max(0.0, -reflect(light, n).z), 32.0) * sqrt(diffuse);\
      gl_FragColor = vec4(mix(color * (0.55 + 0.7 * diffuse), vec3(1.0), specular), 1.0);\
    }\
  ');

  $('#viewer1').bind('mousewheel', function(e, delta, deltaX, deltaY) {
    if (that.hasMesh()){
      if (e.shiftKey){
        if (e.altKey){
          depth -= delta*0.5;
        } else {
          depth -= delta;
        }
      } else {
        depth -= delta*5;  
      }      
      viewer.gl.ondraw();
      e.preventDefault();
    }
  });

  gl.onmousemove = function(e) {
    if (that.hasMesh()){
      if (e.dragging) {
        if (e.button == 2){
          xPos += e.deltaX/10;
          yPos -= e.deltaY/10;
          e.preventDefault();
        } else if (e.shiftKey){
          xPos += e.deltaX/10;
          yPos -= e.deltaY/10;
        } else if (e.ctrlKey){
          angleY += e.deltaX * 2;
          if (angleY>=360)angleY=0;
          if (angleY<0)angleY=360;          
        } else {
          angleZ += e.deltaX * 2;
          if (angleZ>=360)angleZ=0;
          if (angleZ<0)angleZ=360;
          angleX += e.deltaY * 2;
          if (angleX>=360)angleX=0;
          if (angleX<0)angleX=360;
        }

        viewer.gl.ondraw();
      }
    }
  };
  
  gl.ondraw = function() {
    gl.makeCurrent();

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.loadIdentity();
    gl.translate(xPos, yPos, -depth);
    gl.rotate(angleX, 1, 0, 0);
    gl.rotate(angleY, 0, 1, 0);
    gl.rotate(angleZ, 0, 0, 1);
    $('#currAngleX').html(angleX);
    $('#currAngleY').html(angleY);
    $('#currAngleZ').html(angleZ);
    $('#currPosX').html(xPos.toFixed(2));
    $('#currPosY').html(yPos.toFixed(2));
    $('#currDepth').html(depth);

    if (!Viewer.lineOverlay) gl.enable(gl.POLYGON_OFFSET_FILL);
    that.lightingShader.draw(that.mesh, gl.TRIANGLES);
    if (!Viewer.lineOverlay) gl.disable(gl.POLYGON_OFFSET_FILL);

    if (Viewer.lineOverlay) gl.disable(gl.DEPTH_TEST);
    gl.enable(gl.BLEND);
    if (showOutlines){
      that.blackShader.draw(that.mesh, gl.LINES);      
    }
    gl.disable(gl.BLEND);
    if (Viewer.lineOverlay) gl.enable(gl.DEPTH_TEST);

    if (showNormals) {
      gl.lineWidth(2);
      gl.begin(gl.LINES);
      gl.color(0, 0, 0); 
      for (var i = 0; i < mesh.vertices.length; i++) {
        var v = mesh.vertices[i]
        var n = mesh.normals[i]
        gl.vertex(v[0],v[1],v[2]);
        gl.vertex(v[0]+n[0],v[1]+n[1],v[2]+n[2]);
      };
      gl.end();
    }

    var minX = parseFloat($('#boundingBoxMinX').val())
    var maxX = parseFloat($('#boundingBoxMaxX').val())
    var minY = parseFloat($('#boundingBoxMinY').val())
    var maxY = parseFloat($('#boundingBoxMaxY').val())
    var minZ = parseFloat($('#boundingBoxMinZ').val())
    var maxZ = parseFloat($('#boundingBoxMaxZ').val())
    var gridX = parseInt($('#gridX').val());
    var gridY = parseInt($('#gridY').val());
    var gridZ = parseInt($('#gridZ').val());
    
    if (showGrid) {
      gl.lineWidth(2);
      gl.begin(gl.LINES);
      gl.color(0.8,0.8,0.8); 
      for (var i = gridMin; i <= gridMax; i++) {
        gl.vertex(i,gridMin,0);
        gl.vertex(i,gridMax,0);
        gl.vertex(gridMin,i,0);
        gl.vertex(gridMax,i,0);
      };
      gl.end();
    }

    if (showBoundingBox) {
      gl.lineWidth(2);
      gl.begin(gl.LINES);
      gl.color(0.8,0.8,0.8);
      

      gl.vertex(minX,minY,minZ);
      gl.vertex(minX,minY,maxZ);

      gl.vertex(minX,minY,minZ);
      gl.vertex(minX,maxY,minZ);

      gl.vertex(minX,minY,minZ);
      gl.vertex(maxX,minY,minZ);

      gl.vertex(maxX,maxY,maxZ);
      gl.vertex(maxX,maxY,minZ);

      gl.vertex(maxX,maxY,maxZ);
      gl.vertex(maxX,minY,maxZ);

      gl.vertex(maxX,maxY,maxZ);
      gl.vertex(minX,maxY,maxZ);

      gl.vertex(maxX,minY,minZ);
      gl.vertex(maxX,minY,maxZ);

      gl.vertex(maxX,minY,minZ);
      gl.vertex(maxX,maxY,minZ);

      gl.vertex(minX,minY,maxZ);
      gl.vertex(maxX,minY,maxZ);

      gl.vertex(minX,maxY,minZ);
      gl.vertex(maxX,maxY,minZ);

      gl.vertex(minX,maxY,minZ);
      gl.vertex(minX,maxY,maxZ);

      gl.vertex(minX,minY,maxZ);
      gl.vertex(minX,maxY,maxZ);

      gl.end();
    }

    if (showAxis) {
      gl.begin(gl.LINES);
      gl.lineWidth(5.0);
      gl.color(1,0,0); 
      gl.vertex(gridMin,0,0);
      gl.vertex(gridMax,0,0);
      gl.vertex(gridMax-0.25,-0.25,0);
      gl.vertex(gridMax-0.25,0.25,0);
      for (var i = gridMin; i < gridMax; i++){
        gl.vertex(i,-0.1,0);
        gl.vertex(i,0.1,0);
      }
      gl.color(0,1.0,0); 
      gl.vertex(0,gridMin,0);
      gl.vertex(0,gridMax,0);
      gl.vertex(-0.25,gridMax-0.25,0);
      gl.vertex(0.25,gridMax-0.25,0);
      for (var i = gridMin; i < gridMax; i++){
        gl.vertex(-0.1,i,0);
        gl.vertex(0.1,i,0);
      }
      gl.color(0,0,1.0); 
      gl.vertex(0,0,gridMin);
      gl.vertex(0,0,gridMax);
      gl.vertex(-0.25,0,gridMax-0.25);
      gl.vertex(0.25,0,gridMax-0.25);
      for (var i = gridMin; i < gridMax; i++){
        gl.vertex(-0.1,0,i);
        gl.vertex(0.1,0,i);
      }
      gl.end();
    }
    
  };

  this.updateMesh = function(m){
    this.mesh = m;
    gl.ondraw();
  }

  this.hasMesh = function(){
    return this.mesh != undefined; 
  }

  this.cameraAngle = function(X,Y,Z){
    angleX = X;
    angleY = Y;
    angleZ = Z;
    gl.ondraw();
  }

  this.cameraPosition = function(X,Y){
    xPos = X;
    yPos = Y;
    gl.ondraw();
  }

  this.cameraZoom = function(d){
    depth = d || 20;
    gl.ondraw();
  }

  this.canvasResize = function(width,height){
    gl.canvas.width = width;
    gl.canvas.height = height;
    gl.viewport(0, 0, width, height);
    gl.matrixMode(gl.PROJECTION);
    gl.loadIdentity();
    gl.perspective(45, width / height, 0.5, 1000);
    gl.matrixMode(gl.MODELVIEW);
    gl.ondraw();
  }

  $('#'+id).append(this.gl.canvas);
}
