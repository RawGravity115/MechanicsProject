package org.opensourcephysics.sip.ch08;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.DisplayFrame;

import java.io.*;
import java.util.*;
import java.awt.Graphics;

public class PackingApp extends AbstractSimulation {

  double dt;               // time step
  int    iterations;       // Verlet sub‐steps per doStep()
  int    N;                // number of disks
  double k;                // repulsion spring constant
  double gamma;            // linear drag
  double Lx, Ly;           // box dimensions
  double stopTime;         // maximum sim time (seconds)

  boolean readFromFile;    // read initial state?
  String  inFile, outFile; // file names for config and log


  final List<Disk> disks = new ArrayList<>();
  final Random rng      = new Random();
  double t = 0;          
  PrintWriter writer;    

  DisplayFrame frame = new DisplayFrame("x", "y", "Packing Problem");

  class Disk extends Circle {
    double r, vx, vy, ax, ay;
    Disk(double x, double y, double r) {
      this.r = r;
      setXY(x, y);
      color = java.awt.Color.BLUE;
    }
    @Override
    public void draw(DrawingPanel panel, Graphics g) {
      int px = panel.xToPix(getX());
      int py = panel.yToPix(getY());
      int rr = panel.xToPix(r) - panel.xToPix(0);
      g.setColor(color);
      g.drawOval(px - rr, py - rr, 2 * rr, 2 * rr);
    }
  }


  @Override
  public void reset() {
    control.setValue("dt",                0.01);
    control.setValue("T",                 1);
    control.setValue("N",                 30);
    control.setValue("spring k",          1.0e3);
    control.setValue("gamma",             1.0);
    control.setAdjustableValue("Lx",      10.0);
    control.setAdjustableValue("Ly",      10.0);
    control.setValue("stop time (s)",     10.0);
    control.setValue("read config from file", false);
    control.setValue("input file",        "packing_in.txt");
    control.setValue("output file",       "packing_out.csv");
    enableStepsPerDisplay(true);
    initialize();
  }

  @Override
  public void initialize() {
    
    dt           = control.getDouble("dt");
    iterations   = control.getInt   ("T");
    N            = control.getInt   ("N");
    k            = control.getDouble("spring k");
    gamma        = control.getDouble("gamma");
    Lx           = control.getDouble("Lx");
    Ly           = control.getDouble("Ly");
    stopTime     = control.getDouble("stop time (s)");
    readFromFile = control.getBoolean("read config from file");
    inFile       = control.getString ("input file");
    outFile      = control.getString ("output file");

    
    disks.clear();
    frame.clearDrawables();
    if (readFromFile) {
      loadConfiguration(inFile);
    } else {
      randomConfiguration();
    }
    for (Disk d : disks) {
      frame.addDrawable(d);
    }
    frame.setPreferredMinMax(0, Lx, 0, Ly);
    frame.setSquareAspect(true);

  
    try {
      writer = new PrintWriter(outFile);
      writer.printf("k=%.3e, gamma=%.3f, dt=%.4f, T=%d, stopTime=%.4f%n",
                    k, gamma, dt, iterations, stopTime);
      writer.println("t,i,x,y,r,U,xMin,xMax,yMin,yMax,width,height,area");
    } catch (IOException e) {
      control.println("Cannot write " + outFile + ": " + e);
      writer = null;
    }

    t = 0;  
  }

  @Override
  public void stop() {
    if (writer != null) {
      writer.close();
      writer = null;
    }
  }

  
  @Override
  protected void doStep() {
    for (int step = 0; step < iterations; step++) {
      if (t >= stopTime) {          
        stopSimulation();
        return;
      }
    
      double dtStep = dt;
      if (t + dtStep > stopTime) {
        dtStep = stopTime - t;
      }
     
      computeAccelerations();
      for (Disk d : disks) {
        d.vx += 0.5 * d.ax * dtStep;
        d.vy += 0.5 * d.ay * dtStep;
        d.setXY(d.getX() + d.vx * dtStep,
                d.getY() + d.vy * dtStep);
      }
      computeAccelerations();
      for (Disk d : disks) {
        d.vx += 0.5 * d.ax * dtStep;
        d.vy += 0.5 * d.ay * dtStep;
      }
      t += dtStep;                  
      dumpState();
      if (dtStep != dt) {          
        frame.setMessage(String.format("t = %.4f (stopped)", t));
        stopSimulation();
        return;
      }
    }
    frame.setMessage(String.format("t = %.4f", t));
  }


  void computeAccelerations() {
    for (Disk d : disks) {
      d.ax = d.ay = 0;
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
          double f  = k * overlap;
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
    // linear drag
    for (Disk d : disks) {
      d.ax -= gamma * d.vx;
      d.ay -= gamma * d.vy;
    }
  }


  double computeTotalPotentialEnergy() {
    double U = 0;
    for (int i = 0; i < disks.size(); i++) {
      Disk A = disks.get(i);
      for (int j = i + 1; j < disks.size(); j++) {
        Disk B = disks.get(j);
        // manual distance
        double dx = A.getX() - B.getX();
        double dy = A.getY() - B.getY();
        double d  = Math.hypot(dx, dy);
        double ov = (A.r + B.r) - d;
        if (ov > 0) {
          U += 0.5 * k * ov * ov;
        }
      }

      double left   = A.r - A.getX();
      double right  = A.getX() - (Lx - A.r);
      double bottom = A.r - A.getY();
      double top    = A.getY() - (Ly - A.r);
      if (left   > 0) U += 0.5 * k * left * left;
      if (right  > 0) U += 0.5 * k * right * right;
      if (bottom > 0) U += 0.5 * k * bottom * bottom;
      if (top    > 0) U += 0.5 * k * top * top;
    }
    return U;
  }

 
  double[] computeBoundingBox() {
    double xMin = Double.POSITIVE_INFINITY, xMax = Double.NEGATIVE_INFINITY;
    double yMin = Double.POSITIVE_INFINITY, yMax = Double.NEGATIVE_INFINITY;
    for (Disk d : disks) {
      xMin = Math.min(xMin, d.getX() - d.r);
      xMax = Math.max(xMax, d.getX() + d.r);
      yMin = Math.min(yMin, d.getY() - d.r);
      yMax = Math.max(yMax, d.getY() + d.r);
    }
    double width  = xMax - xMin;
    double height = yMax - yMin;
    double area   = width * height;
    return new double[]{xMin, xMax, yMin, yMax, width, height, area};
  }

  
  void dumpState() {
    if (writer == null) return;
    for (int i = 0; i < disks.size(); i++) {
      Disk d = disks.get(i);
      writer.printf(Locale.US,
        "%.5f,%d,%.5f,%.5f,%.4f,",
        t, i, d.getX(), d.getY(), d.r);
      if (i == 0) {
        double U  = computeTotalPotentialEnergy();
        double[] bb = computeBoundingBox();
        writer.printf(Locale.US,
          "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f%n",
          U, bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);
      } else {
        writer.println();
      }
    }
  }

  void randomConfiguration() {
    while (disks.size() < N) {
      double r = 0.2 + 0.3 * rng.nextDouble();
      boolean good;
      double x, y;
      do {
        good = true;
        x = r + (Lx - 2 * r) * rng.nextDouble();
        y = r + (Ly - 2 * r) * rng.nextDouble();
        for (Disk d : disks) {
          if (Math.hypot(x - d.getX(), y - d.getY()) < r + d.r) {
            good = false;
            break;
          }
        }
      } while (!good);
      Disk d = new Disk(x, y, r);
      d.vx = rng.nextDouble() - 0.5;
      d.vy = rng.nextDouble() - 0.5;
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
      control.println("Could not read " + file + ": " + e);
    }
  }

  public static void main(String[] args) {
    SimulationControl.createApp(new PackingApp());
  }
}
