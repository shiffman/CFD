

int IX(int x, int y) {
  x = constrain(x, 0, NX-1);
  y = constrain(y, 0, NY-1);
  return x + y * NX;
}
class Fluid {
  int size;
  float dt;
  float diff;
  float visc;
  float[] s;
  float[] density;
  float[] vx;
  float[] vy;
  float[] vx0;
  float[] vy0;

  Fluid(float diffusion, float viscosity, float dt) {
    this.dt = dt;
    this.diff = diffusion;
    this.visc = viscosity;
    this.s = new float[NX*NY];
    this.density = new float[NX*NY];
    this.vx = new float[NX*NY];
    this.vy = new float[NX*NY];
    this.vx0 = new float[NX*NY];
    this.vy0 = new float[NX*NY];
  }

  void step() {
    float visc     = this.visc;
    float diff     = this.diff;
    float dt       = this.dt;
    float[] vx      = this.vx;
    float[] vy      = this.vy;
    float[] vx0     = this.vx0;
    float[] vy0     = this.vy0;
    float[] s       = this.s;
    float[] density = this.density;

    diffuse(1, vx0, vx, visc, dt, iter);
    diffuse(2, vy0, vy, visc, dt, iter);

    project(vx0, vy0, vx, vy, iter);

    advect(1, vx, vx0, vx0, vy0, dt);
    advect(2, vy, vy0, vx0, vy0, dt);

    project(vx, vy, vx0, vy0, iter);

    diffuse(0, s, density, diff, dt, iter);
    advect(0, density, s, vx, vy, dt);
  }

  void addDensity(int x, int y, float amount) {
    this.density[IX(x, y)] += amount;
  }

  void addVelocity(int x, int y, float amountX, float amountY) {
    int index = IX(x, y);
    this.vx[index] += amountX;
    this.vy[index] += amountY;
  }

  void renderV() {
    for (int j = 0; j < NY; j++) {
      for (int i = 0; i < NX; i++) {
        float x = i * scl;
        float y = j * scl;
        stroke(255);
        strokeWeight(1);
        float vx = this.vx[IX(i, j)];
        float vy = this.vy[IX(i, j)];
        if (vx + vy > 0.05) {
          line(x, y, x + scl*vx, y + scl*vy);
        }
      }
    }
  }

  void renderD() {
    for (int j = 0; j < NY; j++) {
      for (int i = 0; i < NX; i++) {
        float x = i * scl;
        float y = j * scl;
        noStroke();
        fill(255, this.density[IX(i, j)]);
        rect(x, y, scl, scl);
      }
    }
  }
}


void set_bnd ( int b, float[] x ) {
  for (int i=1; i<NY-1; i++ ) {
    x[IX(0, i)] = b==1 ? -x[IX(1, i)] : x[IX(1, i)];
    x[IX(NX-1, i)] = b==1 ? -x[IX(NX-1, i)] : x[IX(NX-1, i)];
  }
  for (int i=1; i<NX-1; i++ ) {
    x[IX(i, 0  )] = b==2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, NY-1)] = b==2 ? -x[IX(i, NY-1)] : x[IX(i, NY-1)];
  }
  x[IX(0, 0  )] = 0.5f*(x[IX(1, 0  )]+x[IX(0, 1)]);
  x[IX(0, NY-1)] = 0.5f*(x[IX(1, NY-1)]+x[IX(0, NY-1)]);
  x[IX(NX-1, 0  )] = 0.5f*(x[IX(NX-2, 0  )]+x[IX(NX-1, 1)]);
  x[IX(NX-1, NY-1)] = 0.5f*(x[IX(NX-2, NY-1)]+x[IX(NX-1, NY-2)]);
}

void lin_solve(int b, float[]x, float[] x0, float a, float c, int iter) {
  float cRecip = 1.0 / c;
  for (int k = 0; k < iter; k++) {
    for (int j = 1; j < NY - 1; j++) {
      for (int i = 1; i < NX - 1; i++) {
        x[IX(i, j)] =
          (x0[IX(i, j)]
          + a*
          (x[IX(i+1, j)]
          +x[IX(i-1, j)]
          +x[IX(i, j+1)]
          +x[IX(i, j-1)]
          )) * cRecip;
      }
    }
  }
  set_bnd(b, x);
}


void diffuse (int b, float[] x, float[] x0, float diff, float dt, int iter) {
  float a = dt * diff * (NX - 2) * (NY - 2);
  lin_solve(b, x, x0, a, 1 + 6 * a, iter);
}

void project(float[] velocX, float[] velocY, float[] p, float[] div, int iter) {
  for (int j = 1; j < NY - 1; j++) {
    for (int i = 1; i < NX - 1; i++) {
      div[IX(i, j)] = -0.5f*(
        velocX[IX(i+1, j)]
        -velocX[IX(i-1, j)]
        +velocY[IX(i, j+1)]
        -velocY[IX(i, j-1)]
        )/((NX+NY)*0.5);
      p[IX(i, j)] = 0;
    }
  }
  set_bnd(0, div); 
  set_bnd(0, p);
  lin_solve(0, p, div, 1, 6, iter);

  for (int j = 1; j < NY - 1; j++) {
    for (int i = 1; i < NX - 1; i++) {
      velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
        -p[IX(i-1, j)]) * NX;
      velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
        -p[IX(i, j-1)]) * NY;
    }
  }

  set_bnd(1, velocX);
  set_bnd(2, velocY);
}

void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt) {
  float i0, i1, j0, j1;

  float dtx = dt * (NX - 2);
  float dty = dt * (NY - 2);

  float s0, s1, t0, t1;
  float tmp1, tmp2, x, y;

  float NXfloat = NX;
  float NYfloat = NY;
  float ifloat, jfloat;
  int i, j;

  for (j = 1, jfloat = 1; j < NY - 1; j++, jfloat++) { 
    for (i = 1, ifloat = 1; i < NX - 1; i++, ifloat++) {
      tmp1 = dtx * velocX[IX(i, j)];
      tmp2 = dty * velocY[IX(i, j)];
      x    = ifloat - tmp1; 
      y    = jfloat - tmp2;

      if (x < 0.5f) x = 0.5f; 
      if (x > NXfloat + 0.5f) x = NXfloat + 0.5f; 
      i0 = floor(x); 
      i1 = i0 + 1.0f;
      if (y < 0.5f) y = 0.5f; 
      if (y > NYfloat + 0.5f) y = NYfloat + 0.5f; 
      j0 = floor(y);
      j1 = j0 + 1.0f; 

      s1 = x - i0; 
      s0 = 1.0f - s1; 
      t1 = y - j0; 
      t0 = 1.0f - t1;

      int i0i = int(i0);
      int i1i = int(i1);
      int j0i = int(j0);
      int j1i = int(j1);

      d[IX(i, j)] = 
        s0 * ( t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) + 
        s1 * ( t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
    }
  }

  set_bnd(b, d);
}
