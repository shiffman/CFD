// Attempting to port
// Real-Time Fluid Dynamics for Games
// Jos Stam
// http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf

final int N = 64;
final int W = N+2;
final int H = N+2;
final int size = W*H;
float[] u = new float[size];
float[] v = new float[size];
float[] u_prev = new float[size];
float[] v_prev = new float[size];
float[] dens = new float[size];
float[] dens_prev = new float[size];

float dt = 0.1f;
float diff = 0.0f;
float visc = 0.0f;
float force = 1.0f;
float source = 100.0f;

void settings() {
  size(W*10, H*10);
  for (int i = 0; i < size; i++) {
    //PVector vec = PVector.random2D();
    //u[i] = vec.x;
    //v[i] = vec.y;
    //u_prev[i] = vec.x;
    //v_prev[i] = vec.y;
    //dens[i] = random(1);
    //dens_prev[i] = random(1);
  }
}

void setup() {
}


void get_from_UI ( float[] d, float[] u, float[] v ) {

  //for (int i=0; i<size; i++ ) {
  //  u[i] = v[i] = d[i] = 0.0f;
  //}

  if (!mousePressed)  return;

  int i = mouseX/10;
  int j = mouseY/10;

  if ( i<1 || i>N || j<1 || j>N ) return;

  if (mousePressed) {
    u[IX(i, j)] = force * (mouseX-pmouseX);
    v[IX(i, j)] = force * (mouseY-pmouseY);
    //println(u[IX(i, j)], v[IX(i, j)]);
  }
}

void mouseReleased() {
  int i = mouseX/10;
  int j = mouseY/10;
  dens_prev[IX(i, j)] = source;
  println(dens_prev[IX(i, j)]);
}

void draw() {
  background(0);
  get_from_UI ( dens_prev, u_prev, v_prev );
  vel_step ( N, u, v, u_prev, v_prev, visc, dt );
  dens_step ( N, dens, dens_prev, u, v, diff, dt );
  velocity();
  //density();
  //println(frameRate);
}

void density() {
  float h = 10;
  beginShape(QUADS);
  for (int i=0; i<=N; i++ ) {
    float x = i*h + h/2;
    for (int j=0; j<=N; j++ ) {
      float y =j*h+j/2;
      noStroke();
      fill(255, dens[IX(i, j)] * 255);
      rect(x, y, h, h);
    }
  }
  endShape();
}

void velocity() {
  float h = 10;
  stroke(255);
  noFill();
  for (int i=1; i<=N; i++ ) {
    float x = i*h+h/2;
    for (int  j=1; j<=N; j++ ) {
      float y = j*h+h/2;
      float x2 = x+u[IX(i, j)]*10;
      float y2 = y+v[IX(i, j)]*10;
      line(x,y,x2,y2);
    }
  }
}
