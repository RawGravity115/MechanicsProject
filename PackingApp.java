
package org.opensourcephysics.sip.ch08;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.DisplayFrame;

import java.io.*;
import java.util.*;


public class PackingApp extends AbstractSimulation {

  
  double dt,                 // time step Δt
         k,                  // spring (repulsion) constant
         gamma,              // linear air-resistance coefficient
         Lx, Ly;             // box size
  int    N;                  // number of particles
  int    iterations;         // T – Verlet steps carried out *every* doStep call
  boolean readFromFile;      // read initial config instead of random?
  String  inFile, outFile;   // file names
  int maxGUIsteps = 2000;   // ≈ 2000 × T × dt physics steps
  int guiSteps = 0;


  final List<Disk> disks = new ArrayList<>();
  final Random rng = new Random();
  double t = 0;              // current time
  PrintWriter writer;        // csv log


  DisplayFrame frame = new DisplayFrame("x", "y", "Packing Problem");


  class Disk extends Circle {
    double r;                // radius
    double vx, vy;           // velocity
    double ax, ay;           // acceleration

    Disk(double x, double y, double r) {
      this.r = r;
      setXY(x, y);
      color = java.awt.Color.BLUE;
    }
  }


  public void initialize() {

    dt         = control.getDouble("dt");
    iterations = control.getInt  ("T");
    N          = control.getInt  ("N");
    k          = control.getDouble("spring k");
    gamma      = control.getDouble("gamma");
    Lx         = control.getDouble("Lx");
    Ly         = control.getDouble("Ly");
    readFromFile = control.getBoolean("read config from file");
    inFile     = control.getString("input file");
    outFile    = control.getString("output file");


    disks.clear();
    frame.clearDrawables();
    if (readFromFile) {
      loadConfiguration(inFile);
    } else {
      randomConfiguration();
    }
    for (Disk d : disks) frame.addDrawable(d);
    frame.setPreferredMinMax(0, Lx, 0, Ly);
    frame.setSquareAspect(true);

    try {
      writer = new PrintWriter(outFile);
      writer.println("# t,i,x,y,r");
    } catch (IOException e) {
      control.println("Cannot write "+outFile+": "+e);
    }
    t = 0;
  }

  public void reset() {
    control.setValue("dt",             0.01);
    control.setValue("T",              1);
    control.setValue("N",              30);
    control.setValue("spring k",       1.0e3);
    control.setValue("gamma",          1.0);
    control.setAdjustableValue("Lx",   10.0);
    control.setAdjustableValue("Ly",   10.0);
    control.setValue("read config from file", false);
    control.setValue("input file",  "packing_in.txt");
    control.setValue("output file", "packing_out.csv");
    enableStepsPerDisplay(true);
    initialize();
  }

  public void stop() {
    if (writer != null) writer.close();
  }

  protected void doStep() {

	  for (int step = 0; step < iterations; step++) {
	    computeAccelerations();

	    for (Disk d : disks) {
	      d.vx += 0.5 * d.ax * dt;
	      d.vy += 0.5 * d.ay * dt;
	      d.setXY(d.getX() + d.vx * dt, d.getY() + d.vy * dt);
	    }

	    computeAccelerations();

	    for (Disk d : disks) {
	      d.vx += 0.5 * d.ax * dt;
	      d.vy += 0.5 * d.ay * dt;
	    }

	    t += dt;
	    dumpState();
	  }

	  frame.setMessage(String.format("t = %.3f", t));

	  if (++guiSteps >= maxGUIsteps) {
	    stopSimulation();
	  }
	}


  void computeAccelerations() {
	  for (Disk d : disks) {
	    d.ax = 0;
	    d.ay = 0;
	  }

	  int n = disks.size();
	  for (int i = 0; i < n; i++) {
	    Disk A = disks.get(i);
	    for (int j = i + 1; j < n; j++) {
	      Disk B = disks.get(j);
	      double dx = B.getX() - A.getX();
	      double dy = B.getY() - A.getY();
	      double dist = Math.hypot(dx, dy);
	      double overlap = (A.r + B.r) - dist;
	      if (overlap > 0) {
	        double f  = k * overlap;         // Hooke’s law
	        double fx = f * dx / dist;
	        double fy = f * dy / dist;
	        A.ax -= fx;  A.ay -= fy;
	        B.ax += fx;  B.ay += fy;
	      }
	    }


	    double left   = A.r - A.getX();
	    double right  = A.getX() - (Lx - A.r);
	    double bottom = A.r - A.getY();
	    double top    = A.getY() - (Ly - A.r);
	    if (left   > 0) A.ax += k * left;
	    if (right  > 0) A.ax -= k * right;
	    if (bottom > 0) A.ay += k * bottom;
	    if (top    > 0) A.ay -= k * top;
	  }

	  for (Disk d : disks) {
	    d.ax -= gamma * d.vx;
	    d.ay -= gamma * d.vy;
	  }
	}



  void dumpState() {
    if (writer == null) return;
    for (int i = 0; i < disks.size(); i++) { // csv
      Disk d = disks.get(i);
      writer.printf(Locale.US, "%.5f,%d,%.5f,%.5f,%.4f%n",
                    t, i, d.getX(), d.getY(), d.r);
    }
  }

  void randomConfiguration() {
    while (disks.size() < N) {
      double r = 0.2 + 0.3 * rng.nextDouble();  
      boolean good;
      double x=0, y=0;
      do {
        good = true;
        x = r + (Lx - 2*r) * rng.nextDouble();
        y = r + (Ly - 2*r) * rng.nextDouble();
        for (Disk d : disks) {
          if (Math.hypot(x - d.getX(), y - d.getY()) < r + d.r) { good = false; break; }
        }
      } while (!good);
      
      Disk d = new Disk(x, y, r);
      d.vx = (rng.nextDouble() - 0.5);     
      d.vy = (rng.nextDouble() - 0.5);    
      disks.add(d);
    }
  }

  void loadConfiguration(String file) {
    try (BufferedReader br = new BufferedReader(new FileReader(file))) {
      String s;
      while ((s = br.readLine()) != null) {
        s = s.trim();
        if (s.isEmpty() || s.startsWith("#")) continue;
        String[] tok = s.split("[,\\s]+");
        if (tok.length < 3) continue;
        double x = Double.parseDouble(tok[0]);
        double y = Double.parseDouble(tok[1]);
        double r = Double.parseDouble(tok[2]);
        disks.add(new Disk(x, y, r));
      }
    } catch (IOException e) {
      control.println("Could not read "+file+": "+e);
    }
  }

  public static void main(String[] args) {
    SimulationControl.createApp(new PackingApp());
  }
}
