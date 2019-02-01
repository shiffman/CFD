final int scl = 4;
final int NX =128;
final int NY =128;
final int iter = 4;

Fluid fluid;
void settings() {
  size(NX*scl, NY*scl);
}

void setup() {
  fluid = new Fluid(0, 0.000001, 0.01);
}


//void mouseDragged() {
//  float dx = mouseX - pmouseX;
//  float dy = mouseY - pmouseY;
//  fluid.addVelocity(mouseX/scl, mouseY/scl, dx*5, dy*5);
//}


void draw() {
  background(0);
  if (mousePressed) {
    for (int i = 0; i < 5; i++) {
      //float a = random(-PI/4, PI/4);
      //PVector v = PVector.fromAngle(a);
      PVector v = new PVector(mouseX - pmouseX, mouseY-pmouseY);
      v.mult(2);
      int x = mouseX/scl + int(random(-2, 3));
      int y = mouseY/scl + int(random(-2, 3));
      fluid.addVelocity(x, y, v.x, v.y);
    }

    for (int x = mouseX-2; x < mouseX+2; x++) {
      for (int y = mouseY-2; y < mouseY+2; y++) {
        fluid.addDensity(x/scl, y/scl, random(10, 25));
      }
    }
  }

  fluid.step();
  //fluid.renderV();
  fluid.renderD();
  
  fluid.fadeD();
}
