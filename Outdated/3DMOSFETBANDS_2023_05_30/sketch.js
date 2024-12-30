function preload() {
  homeIcon = loadImage('Assets/HomeIcon.png');
  homeIconHover = loadImage('Assets/HomeIconHover.png');
  transverseIcon = loadImage('Assets/TransverseIcon.png');
  transverseIconHover = loadImage('Assets/TransverseIconHover.png');
  longitudinalIcon = loadImage('Assets/LongitudinalIcon.png');
  longitudinalIconHover = loadImage('Assets/LongitudinalIconHover.png');

  tox = 3e-9;
  na =1e17;
  na = na*1e6; 
  Nv = 9.84e18;
  Nv=Nv*1e6; 
  Nc = 2.78e19;
  Nc = Nc*1e6; 
  T = 300;
  Eg = 1.12;
  phi_m = 4.08;
  chi_s = 4.05;
  eps = 8.854187817E-12*12;
  eps_ox = 8.854187817E-12*4;
  W = 10e-6;
  L = 1e-6;
  u = 450/(10**4);
  C = (eps_ox)/tox;
  kp = ((W*C*u)/L)
  Vds = 1.9;
  Vgs = 3;
  Va = Vgs;
  
  Ev = [];  //valence band energy
  // d arrays used for plotting
  d1 = [];  //Ev
  d2 = [];  //EF
  d3 = [];  //Ec
  d4 = [];  //n
  d5 = [];  //p
  d6 = [];  //rho
  d7 = [];  //Electric field
  d8 = [];  //Ev_ch
  d9 = [];  //Ec_ch
  transverse_valence_bands = [];
  transverse_conduction_bands = [];
  echarge = 1.60217733E-19;
  kb = 1.380658E-23/echarge;  //eV/K
  
  ni=Math.sqrt(Nc*Nv*Math.pow(T/300,3))*Math.exp(1.6022E-19*(-Eg)/(2*1.380658E-23*T));
  VT = 2*tox*Math.sqrt(eps*na*1.380658E-23*T*Math.log(na/ni))/eps_ox+2*kb*T*Math.log(na/ni);  //threshold voltage in depletion approximation
  Ev0 = kb*T*Math.log(na/Nv); //valence band far from the oxide

  vertices = 10;
  num_points = 10*vertices;
}

function setup() {
  createCanvas(600, 600, WEBGL);
  moscap_first();
  calculate_Vch();
  
  for (iterator=0; iterator<vertices+1; iterator++){
    Ev_goal = d8[iterator][1];
    moscap();
    transverse_valence_bands[iterator] = d1;
    transverse_conduction_bands[iterator] = d3;
  }
  // print(transverse_valence_bands);
  
  cam = createEasyCam();
  document.oncontextmenu = ()=>false;
  angle = 35.264/180*3.14;
  homeSize = 70;
  bandButtonX = 90;
  bandButtonY = 50;
  homeIcon.resize(homeSize, homeSize);
  homeIconHover.resize(homeSize, homeSize);
  transverseIcon.resize(bandButtonX, bandButtonY);
  transverseIconHover.resize(bandButtonX, bandButtonY);
  longitudinalIcon.resize(bandButtonX, bandButtonY);
  longitudinalIconHover.resize(bandButtonX, bandButtonY);
  homeOver = false;
  transverseOver = false;
  longitudinalOver = false;
  homeX = width-homeSize-25;
  homeY = 20;
  transverseX = width-homeSize-bandButtonX-50;
  transverseY = homeY + 0.5*homeSize - bandButtonY/2;
  longitudinalX = homeX + homeSize/2 - bandButtonX/2;
  longitudinalY = 50 + homeSize;
  frameRate(30); 
  cam.setDistanceMin(10);
  cam.setDistanceMax(1000);
  cam.panX(100);
  cam.panY(40);
  cam.mouseDragPan();
  cam.rotateX(-HALF_PI)
  longitudinal = cam.getState();
  cam.rotateY(-HALF_PI);
  transverse = cam.getState();
  cam.rotateY(HALF_PI);
  cam.rotateX(HALF_PI);
  cam.rotateX(-angle);  // rotate around the x-axis passing through the subject
  cam.rotateY(-angle);  // rotate around the y-axis passing through the subject
  cam.rotateZ(QUARTER_PI);  // rotate around the z-axis passing through the subject
  home = cam.getState();
}

function draw() {
  update(mouseX, mouseY);
  background(100);
  lights();
  translate(0,0,100);
  push();
  for(i = 0; i < vertices; i++){ 
    for(j = num_points - vertices; j < num_points; j++){
      beginShape();
      fill(0,0,255,150);
      vertex(200*d8[i][0], 10000*transverse_valence_bands[i][j][0], 50*transverse_valence_bands[i][j][1]);
      vertex(200*d8[i][0], 10000*transverse_valence_bands[i][j+1][0], 50*transverse_valence_bands[i][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_valence_bands[i+1][j+1][0], 50*transverse_valence_bands[i+1][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_valence_bands[i+1][j][0], 50*transverse_valence_bands[i+1][j][1]);
      vertex(200*d8[i][0], 10000*transverse_valence_bands[i][j][0], 50*transverse_valence_bands[i][j][1]);
      endShape();
    }
  }
  
  for(i = 0; i < vertices; i++){ 
    for(j = num_points - vertices; j < num_points; j++){
      beginShape();
      fill(255,0,0,150);
      vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][j][0], 50*transverse_conduction_bands[i][j][1]);
      vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][j+1][0], 50*transverse_conduction_bands[i][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_conduction_bands[i+1][j+1][0], 50*transverse_conduction_bands[i+1][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_conduction_bands[i+1][j][0], 50*transverse_conduction_bands[i+1][j][1]);
      vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][j][0], 50*transverse_conduction_bands[i][j][1]);
      endShape();
    }
  }
  
  fill(255,255,255);
   for(i = 0; i < vertices; i++){ 
    for(j = num_points; j < transverse_valence_bands[i].length - 1; j++){
      beginShape();
      vertex(200*d8[i][0], 10000*transverse_valence_bands[i][j][0], 50*transverse_valence_bands[i][j][1]);
      vertex(200*d8[i][0], 10000*transverse_valence_bands[i][j+1][0], 50*transverse_valence_bands[i][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_valence_bands[i+1][j+1][0], 50*transverse_valence_bands[i+1][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_valence_bands[i+1][j][0], 50*transverse_valence_bands[i+1][j][1]);
      vertex(200*d8[i][0], 10000*transverse_valence_bands[i][j][0], 50*transverse_valence_bands[i][j][1]);
      endShape();
    }
  }
  for(i = 0; i < vertices; i++){ 
    for(j = num_points; j < transverse_conduction_bands[i].length -1; j++){
      beginShape();
      vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][j][0], 50*transverse_conduction_bands[i][j][1]);
      vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][j+1][0], 50*transverse_conduction_bands[i][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_conduction_bands[i+1][j+1][0], 50*transverse_conduction_bands[i+1][j+1][1]);
      vertex(200*d8[i+1][0], 10000*transverse_conduction_bands[i+1][j][0], 50*transverse_conduction_bands[i+1][j][1]);
      vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][j][0], 50*transverse_conduction_bands[i][j][1]);
      endShape();
    }
  }
  
  for(i = 0; i < vertices; i++){
      beginShape();
    vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][num_points][0], 50*transverse_conduction_bands[i][num_points][1]);
        vertex(200*d8[i+1][0], 10000*transverse_conduction_bands[i+1][num_points][0], 50*transverse_conduction_bands[i+1][num_points][1]);
       vertex(200*d8[i+1][0], 10000*transverse_valence_bands[i+1][num_points][0], 50*transverse_valence_bands[i+1][num_points][1]);
        vertex(200*d8[i][0], 10000*transverse_valence_bands[i][num_points][0], 50*transverse_valence_bands[i][num_points][1]);
        vertex(200*d8[i][0], 10000*transverse_conduction_bands[i][num_points][0], 50*transverse_conduction_bands[i][num_points][1]);
    endShape();
  }
  
  for (i = 0; i<2; i++){
    beginShape();
    vertex(200*d8[vertices*i][0], 10000*transverse_conduction_bands[vertices*i][num_points+1][0], 50*transverse_conduction_bands[vertices*i][num_points+1][1]);
        vertex(200*d8[vertices*i][0], 10000*transverse_conduction_bands[vertices*i][num_points+2][0], 50*transverse_conduction_bands[vertices*i][num_points+2][1]);
       vertex(200*d8[vertices*i][0], 10000*transverse_valence_bands[vertices*i][num_points+2][0], 50*transverse_valence_bands[vertices*i][num_points+2][1]);
        vertex(200*d8[vertices*i][0], 10000*transverse_valence_bands[vertices*i][num_points+1][0], 50*transverse_valence_bands[vertices*i][num_points+1][1]);
        vertex(200*d8[vertices*i][0], 10000*transverse_conduction_bands[vertices*i][num_points+1][0], 50*transverse_conduction_bands[vertices*i][num_points+1][1]);
    endShape();
  }
  
  
  
  

  cam.beginHUD();  
   if (homeOver) {
    image(homeIconHover, homeX, homeY);
  } else {
    image(homeIcon, homeX, homeY);
  }
  if (longitudinalOver) {
    image(longitudinalIconHover, longitudinalX, longitudinalY);
  } else {
    image(longitudinalIcon, longitudinalX, longitudinalY);
  }
  if (transverseOver) {
    image(transverseIconHover, transverseX, transverseY);
  } else {
    image(transverseIcon, transverseX, transverseY);
  }
  cam.endHUD();
}

function update(x, y) {
  if (over(homeX, homeY, homeSize, homeSize) ) {
    homeOver = true;
  } else {
    homeOver = false;
  }
  if (over(longitudinalX, longitudinalY, bandButtonX, bandButtonY) ) {
    longitudinalOver = true;
  } else {
    longitudinalOver = false;
  }
  if (over(transverseX, transverseY, bandButtonX, bandButtonY) ){
    transverseOver = true;
  } else {
    transverseOver = false;
  }
}

function mousePressed() {
  if (homeOver) {
    cam.setState(home);
  }
  if (longitudinalOver) {
    cam.setState(longitudinal);
  }
  if (transverseOver) {
    cam.setState(transverse);
  }
}

function over(x, y, width, height)  {
  if (mouseX >= x && mouseX <= x+width && 
      mouseY >= y && mouseY <= y+height) {
    return true;
  } else {
    return false;
  }
}

function moscap_first() { 
  phi_s = chi_s + Eg + Ev0; // work function of semiconductor 
  Vfb = phi_m  - phi_s; // flat band voltage
  xp = 2*Math.sqrt(eps*kb*T*Math.log(na/ni)/(echarge*na)); // maximum depletion width in depletion approximation
  cox= eps_ox/tox;
  V = Va - phi_m + phi_s; //zero bias should result in a built-in voltage corresponding to the work function difference

  /* This section starts with somewhat arbitray initial conditions for Ev(x) and its derivative
   dEv(x)/dx inside the bulk of the semiconductor and then integrates the Possion equation
   to the left to calculate the voltage at the semiconductor/oxide interface and then the 
   voltage at the gate. The calculated voltage at the gate is compared to the desired voltage at the gate 
   the initial conditions are adjusted to get the calculated value closer to the desired value.
  */
  xr = 3.6*xp;   //right position for binary search
  xm = 1.8*xp; //middle position for binary search
  xl = 0;      //left position for binary search
  for (j=0; j<16; j++) {
    Ev = Ev0;
    dEvdx = 0;
    if (V>0) {dEvdx = 1e3; } 
    if (V<0) {dEvdx = -1e3; } 
    dx = xm/num_points;
    for (i=0; i<num_points+1; i++) { //numerical integration using the midpoint method
      x = xm*(1-i/num_points);
      Ev_mid = Ev - dEvdx*dx/2;  //Ev at the midpoint
      // if (j == 0){
      //   print(i);
      //   print(Ev_mid);
      // }
      rho_mid = -na - Nc*(1/(1+Math.exp((Eg+Ev_mid)/(kb*T)))) + Nv*(1/(1+Math.exp(-Ev_mid/(kb*T))));
      //print(Math.exp(-(Eg+Ev_mid)/(kb*T)));
      // print(Eg+Ev_mid)
      dEvdx_mid = dEvdx - echarge*rho_mid*0.5*dx/eps;
      Ev = Ev - dEvdx_mid*dx;
      dEvdx = dEvdx - echarge*rho_mid*dx/eps;
      // if (j==2){
      //   print(xm);
      //   print(dx);
      //   print(dEvdx_mid);
      // }
      // print(dEvdx)
    }
    Vs = Ev0 - Ev;
    Es = dEvdx;
    Eox = Es*eps/eps_ox;
    Vtmp = Vs + Eox*tox;

  //binary search
    if (V<0) {  
      if (Vtmp > V) {xl = xm; xm = (xr + xl)/2; }
      if (Vtmp < V) {xr = xm; xm = (xr + xl)/2; } 
    }
    else {
      if (Vtmp < V) {xl = xm; xm = (xr + xl)/2; }
      if (Vtmp > V) {xr = xm; xm = (xr + xl)/2; } 
    }
    // print(xm);
    // print(Vs);
    // print(Vtmp)
  } // end j

  /* This section solves the Possion equation one more time using the optimized initial conditions
   and saves the data in arrays so it can be plotted.
  */
    Ev = Ev0;
    dEvdx = 0;
    if (V>0) {dEvdx = 1e3; } 
    if (V<0) {dEvdx = -1e3; } 
    dx = xm/num_points;
    for (i=0; i<num_points + 1; i++) { //numerical integration using the midpoint method
      x = xm*(1-i/num_points);
      Ev_mid = Ev - dEvdx*dx/2;  //Ev at the midpoint
      rho_mid = -na - Nc*(1/(1+Math.exp((Eg+Ev_mid)/(kb*T)))) + Nv*(1/(1+Math.exp(-Ev_mid/(kb*T))));
      dEvdx_mid = dEvdx - echarge*rho_mid*0.5*dx/eps;
      Ev = Ev - dEvdx_mid*dx;
      dEvdx = dEvdx - echarge*rho_mid*dx/eps;

      d1[i] = [1E6*x,Ev];
      d3[i] = [1E6*x,Ev+Eg];
      n = Nc*Math.exp(-(Eg+Ev)/(kb*T))/1E21;
      p = Nv*Math.exp(Ev/(kb*T))/1E21
      d4[i] = [1E6*x,n];
      d5[i] = [1E6*x,p];
      d6[i] = [1E6*x,p-n-na/1E21];
      d7[i] = [1E6*x,dEvdx*1E-6];
    }
    d4[num_points + 1] = [0,0]; d4[num_points + 2] = [-1E6*(0.25*xp+tox),0];
    d5[num_points + 1] = [0,0]; d5[num_points + 2] = [-1E6*(0.25*xp+tox),0];
    d6[num_points + 1] = [0,0]; d6[num_points + 2] = [-1E6*tox,0]; d6[num_points + 3] = [-1E6*tox,-d6[num_points][1]];  d6[num_points + 4] = [-1E6*(tox-xm/100),-0.92*d6[num_points][1]]; 
    d6[num_points + 5] = [-1E6*(tox+xm/100),-0.92*d6[num_points][1]]; d6[num_points + 6] = [-1E6*tox,-d6[num_points][1]]; d6[num_points + 7] = [-1E6*tox,0]; d6[num_points + 8] = [-1E6*(0.25*xp+tox),0];
    d2[0] = [1E6*xm,0];
    d2[1] = [0,0];
    VIB = Ev0 - Ev;
    Es = dEvdx;
    Eox = Es*eps/eps_ox;
    Q = -Es*eps; //the charge on the semiconductor
    d7[num_points + 1] = [0,Eox*1E-6]; d7[num_points + 2]=[-tox*1E6,Eox*1E-6];  d7[num_points + 3] = [-tox*1E6,0]; d7[num_points + 4] = [-1E6*(0.25*xp+tox),0];
    Vtmp = VIB + Eox*tox;
    d2[3] = null;
    d2[4] = [-1E6*tox,-Va];
    d2[5] = [-1E6*(0.25*xp+tox),-Va];
    // d3[num_points + 1] = [0,d3[num_points][1]+1];
    // d3[num_points + 2] = [-1E6*tox,d3[num_points][1]+1- Eox*tox];
    // d3[num_points + 3] = [-1E6*tox,-Vtmp];
    // d1[num_points + 1] = [0,d1[num_points][1]-1];
    // d1[num_points + 2] = [-1E6*tox,d1[num_points][1]-1- Eox*tox];
    // d1[num_points + 3] = [-1E6*tox,-Vtmp];
    if (V==0) { //
      d1[0] = [1E6*xp,Ev0]; d2[0] = [1E6*xp,0]; d3[0] = [1E6*xp,Ev0+Eg];
      d4[0] = [1E6*xp,0]; d5[0] = [1E6*xp,na/1E21]; d6[0] = [1E6*xp,0]; 
      d7.length = 0;
      d7[0] = [1E6*xp,0]; d7[1] = [-1E6*(0.25*xp+tox),0];
    }
}

function moscap() { 
  var valence = [];
  var conduction = [];
  phi_s = chi_s + Eg + Ev0; // work function of semiconductor 
  Vfb = phi_m  - phi_s; // flat band voltage
  xp = 2*Math.sqrt(eps*kb*T*Math.log(na/ni)/(echarge*na)); // maximum depletion width in depletion approximation
  cox= eps_ox/tox;
  V = Va - phi_m + phi_s; //zero bias should result in a built-in voltage corresponding to the work function difference

  /* This section starts with somewhat arbitray initial conditions for Ev(x) and its derivative
   dEv(x)/dx inside the bulk of the semiconductor and then integrates the Possion equation
   to the left to calculate the voltage at the semiconductor/oxide interface and then the 
   voltage at the gate. The calculated voltage at the gate is compared to the desired voltage at the gate 
   the initial conditions are adjusted to get the calculated value closer to the desired value.
  */
  xr = 3.6*xp;   //right position for binary search
  xm = 1.8*xp; //middle position for binary search
  xl = 0;      //left position for binary search
  for (j=0; j<16; j++) {
    Ev = Ev0;
    dEvdx = 0;
    if (V>0) {dEvdx = 1e3; } 
    if (V<0) {dEvdx = -1e3; } 
    dx = xm/num_points;
    for (i=0; i<num_points+1; i++) { //numerical integration using the midpoint method
      x = xm*(1-i/num_points);
      Ev_mid = Ev - dEvdx*dx/2;  //Ev at the midpoint
      // if (j == 0){
      //   print(i);
      //   print(Ev_mid);
      // }
      rho_mid = -na - Nc*(1/(1+Math.exp((Eg+Ev_mid)/(kb*T)))) + Nv*(1/(1+Math.exp(-Ev_mid/(kb*T))));
      //print(Math.exp(-(Eg+Ev_mid)/(kb*T)));
      // print(Eg+Ev_mid)
      dEvdx_mid = dEvdx - echarge*rho_mid*0.5*dx/eps;
      Ev = Ev - dEvdx_mid*dx;
      dEvdx = dEvdx - echarge*rho_mid*dx/eps;
      // if (j==2){
      //   print(xm);
      //   print(dx);
      //   print(dEvdx_mid);
      // }
      // print(dEvdx)
    }
    Vs = Ev0 - Ev;
    Es = dEvdx;
    Eox = Es*eps/eps_ox;
    Vtmp = Vs + Eox*tox;

  //binary search
    if (V<0) {  
      if (Ev < Ev_goal) {xl = xm; xm = (xr + xl)/2; }
      if (Ev > Ev_goal) {xr = xm; xm = (xr + xl)/2; } 
    }
    else {
      if (Ev > Ev_goal) {xl = xm; xm = (xr + xl)/2; }
      if (Ev < Ev_goal) {xr = xm; xm = (xr + xl)/2; } 
    }
    // print(xm);
    // print(Vs);
    // print(Vtmp)
  } // end j

  /* This section solves the Possion equation one more time using the optimized initial conditions
   and saves the data in arrays so it can be plotted.
  */
    Ev = Ev0;
    dEvdx = 0;
    if (V>0) {dEvdx = 1e3; } 
    if (V<0) {dEvdx = -1e3; } 
    dx = xm/num_points;
    for (i=0; i<num_points + 1; i++) { //numerical integration using the midpoint method
      x = xm*(1-i/num_points);
      Ev_mid = Ev - dEvdx*dx/2;  //Ev at the midpoint
      rho_mid = -na - Nc*(1/(1+Math.exp((Eg+Ev_mid)/(kb*T)))) + Nv*(1/(1+Math.exp(-Ev_mid/(kb*T))));
      dEvdx_mid = dEvdx - echarge*rho_mid*0.5*dx/eps;
      Ev = Ev - dEvdx_mid*dx;
      dEvdx = dEvdx - echarge*rho_mid*dx/eps;

      valence[i] = [1E6*x,Ev];
      conduction[i] = [1E6*x,Ev+Eg];
      n = Nc*Math.exp(-(Eg+Ev)/(kb*T))/1E21;
      p = Nv*Math.exp(Ev/(kb*T))/1E21
      d4[i] = [1E6*x,n];
      d5[i] = [1E6*x,p];
      d6[i] = [1E6*x,p-n-na/1E21];
      d7[i] = [1E6*x,dEvdx*1E-6];
    }
    d4[num_points + 1] = [0,0]; d4[num_points + 2] = [-1E6*(0.25*xp+tox),0];
    d5[num_points + 1] = [0,0]; d5[num_points + 2] = [-1E6*(0.25*xp+tox),0];
    d6[num_points + 1] = [0,0]; d6[num_points + 2] = [-1E6*tox,0]; d6[num_points + 3] = [-1E6*tox,-d6[num_points][1]];  d6[num_points + 4] = [-1E6*(tox-xm/100),-0.92*d6[num_points][1]]; 
    d6[num_points + 5] = [-1E6*(tox+xm/100),-0.92*d6[num_points][1]]; d6[num_points + 6] = [-1E6*tox,-d6[num_points][1]]; d6[num_points + 7] = [-1E6*tox,0]; d6[num_points + 8] = [-1E6*(0.25*xp+tox),0];
    d2[0] = [1E6*xm,0];
    d2[1] = [0,0];
    VIB = Ev0 - Ev;
    Es = dEvdx;
    Eox = Es*eps/eps_ox;
    Q = -Es*eps; //the charge on the semiconductor
    d7[num_points + 1] = [0,Eox*1E-6]; d7[num_points + 2]=[-tox*1E6,Eox*1E-6];  d7[num_points + 3] = [-tox*1E6,0]; d7[num_points + 4] = [-1E6*(0.25*xp+tox),0];
    Vtmp = VIB + Eox*tox;
    d2[3] = null;
    d2[4] = [-1E6*tox,-Va];
    d2[5] = [-1E6*(0.25*xp+tox),-Va];
    conduction[num_points + 1] = [0,d3[num_points][1]+2];
    conduction[num_points + 2] = [-1E6*tox,d3[num_points][1]+2- Eox*tox];
    conduction[num_points + 3] = [-1E6*tox,-Vtmp];
    valence[num_points + 1] = [0,d1[num_points][1]-2];
    valence[num_points + 2] = [-1E6*tox,d1[num_points][1]-2- Eox*tox];
    valence[num_points + 3] = [-1E6*tox,-Vtmp];
    if (V==0) { //
      valence[0] = [1E6*xp,Ev0]; d2[0] = [1E6*xp,0]; conduction[0] = [1E6*xp,Ev0+Eg];
      d4[0] = [1E6*xp,0]; d5[0] = [1E6*xp,na/1E21]; d6[0] = [1E6*xp,0]; 
      d7.length = 0;
      d7[0] = [1E6*xp,0]; d7[1] = [-1E6*(0.25*xp+tox),0];
    }
  
  d1 = valence;
  d3 = conduction;
}

function calculate_Vch(){
  Ev_ch0 = d1[num_points][1];
  if (Vds<Vgs-VT){
    Id = kp*(Vgs - VT - (Vds/2))*Vds;
  }
  else{
    Id = kp*(Vgs - VT)**2;
  }  v_term = Vgs-VT;
 v_term = Vgs-VT;
  if (v_term>0){
    for (i=0; i<vertices+1; i++) { 
      y = (L*(i/vertices));
      Vch = v_term-sqrt(((v_term)**2)-(2*Id*y)/(W*C*u));
      if (isNaN(Vch)){
        Vch = Vds;
      }
      Ev_ch = Ev_ch0 - Vch;
      d8[i] = [1E6*y,Ev_ch];
      d9[i] = [1E6*y,Ev_ch+Eg];
    }
  }
  else {
    for (i=0; i<vertices+1; i++) { 
      y = (L*(i/vertices));
      d8[i] = [1E6*y,Ev_ch0];
      d9[i] = [1E6*y,Ev_ch0+Eg];
    }
  }
}










