// Fluid Simulation
// Daniel Shiffman
// https://thecodingtrain.com/CodingChallenges/132-fluid-simulation.html
// https://youtu.be/alhpH6ECFvQ

// This would not be possible without:
// Real-Time Fluid Dynamics for Games by Jos Stam
// http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
// Fluid Simulation for Dummies by Mike Ash
// https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

int IX(int x, int y) {
  x = constrain(x, 0, N-1);
  y = constrain(y, 0, N-1);
  return x + (y * N);
}

class Velocity {
  PVector vel;
  PVector prev_vel;
  Velocity() {
    this.vel = new PVector();
    this.prev_vel = new PVector();
  }

  void add(float x, float y) {
    vel.add(x, y);
  }
}

class Dye {
  float dye, prev_dye;

  Dye() {
    this.dye = 0;
    this.prev_dye = 0;
  }

  void add(float amt) {
    this.dye += amt;
  }
} 

class Fluid {
  float dt;
  float diff;
  float visc;

  Velocity[] velocity;
  Dye[] density;

  Fluid() {
    this.velocity = new Velocity[N*N];
    this.density = new Dye[N*N];
    for (int i = 0; i < velocity.length; i++) {
      velocity[i] = new Velocity();
      density[i] = new Dye();
    }

    for (int i = 0; i < velocity.length; i++) {
      velocity[i] = new Velocity();
    }
  }

  void step() {
    diffuse(velocity);
    project(velocity);
    advect(velocity);
    project(velocity);
    diffuse(density);
    advect(density);
  }

  void addDensity(int x, int y, float amount) {
    int index = IX(x, y);
    this.density[index].add(amount);
  }

  void addVelocity(int x, int y, float amountX, float amountY) {
    int index = IX(x, y);
    this.velocity[index].add(amountX, amountY);
  }

  void renderD() {
    colorMode(HSB, 255);

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        float x = i * SCALE;
        float y = j * SCALE;
        float d = this.density[IX(i, j)].dye;
        fill((d + 50) % 255, 200, d);
        noStroke();
        square(x, y, SCALE);
      }
    }
  }

  void renderV() {

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        float x = i * SCALE;
        float y = j * SCALE;
        float vx = this.velocity[IX(i, j)].vel.x;
        float vy = this.velocity[IX(i, j)].vel.y;
        stroke(255);
        if (!(abs(vx) < 0.1 && abs(vy) <= 0.1)) {
          line(x, y, x+vx*SCALE, y+vy*SCALE );
        }
      }
    }
  }

  void fadeD() {
    for (int i = 0; i < this.density.length; i++) {
      float d = density[i].dye;
      density[i].dye = constrain(d - 0.02, 0, 255);
    }
  }
}
