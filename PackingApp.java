package org.opensourcephysics.sip.ch08;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.DisplayFrame;

import java.io.*;
import java.util.*;
import java.util.Locale;
import java.awt.Graphics;
import java.awt.Color;

public class PackingApp extends AbstractSimulation {


  double dt;               // time step
  int    iterations;       // Verlet sub‐steps per doStep()
  int    N;                // number of disks
  double k;                // repulsion spring constant
  double gamma;            // linear drag
  double Lx, Ly;           // box dimensions
  double stopTime;         // max sim time per cycle (s)


  double rMin, rMax;       // min and max disk radii


  double prevArea     = Double.POSITIVE_INFINITY;
  double epsilon      = 1e-3;    // relative area tolerance
  int    cycleCount   = 0;
  int    maxCycles    = 1000;    // hard cap
  int    stableCount  = 0;
  int    stableThreshold = 15;    // require 5 consecutive small changes


  boolean     readFromFile;
  String      inFile, outFile;
  final List<Disk> disks = new ArrayList<>();
  final Random      rng   = new Random();
  double     t       = 0;        // simulation clock
  PrintWriter writer;            // CSV logger


  DisplayFrame frame = new DisplayFrame("x", "y", "Packing Problem");

  class Disk extends Circle {
    double r, vx, vy, ax, ay;
    Disk(double x, double y, double r) {
      this.r = r;
      setXY(x, y);
      color = Color.BLUE;
    }
    @Override
    public void draw(DrawingPanel panel, Graphics g) {
      int px = panel.xToPix(getX()), py = panel.yToPix(getY());
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
    control.setValue("min radius",        0.2);
    control.setValue("max radius",        0.5);

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
    rMin         = control.getDouble("min radius");
    rMax         = control.getDouble("max radius");
    readFromFile = control.getBoolean("read config from file");
    inFile       = control.getString ("input file");
    outFile      = control.getString ("output file");

    disks.clear();
    if (readFromFile) {
      loadConfiguration(inFile);
    } else {
      randomConfiguration();
    }


    try {
      writer = new PrintWriter(outFile);

      writer.printf("k=%.3e, gamma=%.3f, dt=%.4f, T=%d, stopTime=%.4f%n",
                    k, gamma, dt, iterations, stopTime);
      writer.printf("min radius=%.3f, max radius=%.3f%n", rMin, rMax);
      writer.printf("N=%d, Lx=%.4f, Ly=%.4f%n", N, Lx, Ly);
      writer.println();


      writer.println("# initial positions: i, x, y, r");
      writer.println("i,x,y,r");
      for (int i = 0; i < disks.size(); i++) {
        Disk d = disks.get(i);
        writer.printf(Locale.US, "%d,%.5f,%.5f,%.5f%n",
                      i, d.getX(), d.getY(), d.r);
      }
      writer.println();

  
      writer.println("# cycle,area");
      writer.println("cycle,area");
      double[] bb0 = computeBoundingBox();
      prevArea = bb0[6];
      writer.printf("0,%.5f%n", prevArea);

    } catch (IOException e) {
      control.println("Cannot write " + outFile + ": " + e);
      writer = null;
    }


    frame.clearDrawables();
    frame.setPreferredMinMax(0, Lx, 0, Ly);
    frame.setSquareAspect(true);


    cycleCount  = 0;
    stableCount = 0;


    startCycle();
  }


  private void startCycle() {
    t = 0;
    frame.clearDrawables();
    Trail bounds = new Trail();
    bounds.color = Color.BLACK;
    bounds.addPoint(0,  0);
    bounds.addPoint(Lx, 0);
    bounds.addPoint(Lx, Ly);
    bounds.addPoint(0,  Ly);
    bounds.addPoint(0,  0);
    frame.addDrawable(bounds);
    for (Disk d : disks) {
      d.vx = rng.nextDouble() - 0.5;
      d.vy = rng.nextDouble() - 0.5;
      frame.addDrawable(d);
    }
    frame.setMessage("Starting cycle " + (cycleCount+1));
  }


  @Override
  public void doStep() {
    for (int step = 0; step < iterations; step++) {
      double dtStep = dt;
      if (t + dtStep > stopTime) dtStep = stopTime - t;


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

      if (t >= stopTime) {
        handleCycleEnd();
        return;
      }
    }
    frame.setMessage(String.format("t = %.4f (cycle %d)", t, cycleCount));
  }


  private void handleCycleEnd() {
    double[] bb   = computeBoundingBox();  
    double   area = bb[6];


    cycleCount++;
    if (writer != null) {
      writer.printf("%d,%.5f%n", cycleCount, area);
    }

    frame.clearDrawables();
    for (Disk d : disks) {
      double x = d.getX() - bb[0], y = d.getY() - bb[2];
      d.setXY(x, y);
      frame.addDrawable(d);
    }
    Trail box = new Trail();
    box.color = Color.RED;
    box.addPoint(0, 0);
    box.addPoint(bb[4], 0);
    box.addPoint(bb[4], bb[5]);
    box.addPoint(0, bb[5]);
    box.addPoint(0, 0);
    frame.addDrawable(box);


    double relChange = Math.abs(prevArea - area)/prevArea;
    prevArea = area;

    if (relChange < epsilon) {
      stableCount++;
    } else {
      stableCount = 0;
    }


    if (cycleCount >= maxCycles || stableCount >= stableThreshold) {
      if (writer != null) {
        writer.println();
        writer.println("# final positions: i, x, y, r");
        writer.println("i,x,y,r");
        for (int i = 0; i < disks.size(); i++) {
          Disk d = disks.get(i);
          writer.printf(Locale.US, "%d,%.5f,%.5f,%.5f%n",
                        i, d.getX(), d.getY(), d.r);
        }
        writer.printf("%nfinal,%.5f%n", area);
        writer.close();
      }
      frame.setMessage(
        String.format("Converged after %d cycles: area=%.5f", cycleCount, area)
      );
      stopSimulation();
    } else {

      Lx = bb[4];
      Ly = bb[5];
      frame.setMessage(
        String.format("Cycle %d end: area=%.5f (stableCount=%d)",
                      cycleCount, area, stableCount)
      );
      startCycle();
    }
  }



  private void computeAccelerations() {
    for (Disk d : disks) { d.ax = d.ay = 0; }
    int n = disks.size();
    for (int i = 0; i < n; i++) {
      Disk A = disks.get(i);
      // pairwise repulsion
      for (int j = i+1; j < n; j++) {
        Disk B = disks.get(j);
        double dx = B.getX() - A.getX(), dy = B.getY() - A.getY();
        double dist = Math.hypot(dx, dy), overlap = (A.r + B.r) - dist;
        if (overlap > 0) {
          double f = k * overlap, fx = f * dx / dist, fy = f * dy / dist;
          A.ax -= fx; A.ay -= fy;
          B.ax += fx; B.ay += fy;
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

  private double[] computeBoundingBox() {
    double xMin = Double.POSITIVE_INFINITY, xMax = Double.NEGATIVE_INFINITY;
    double yMin = Double.POSITIVE_INFINITY, yMax = Double.NEGATIVE_INFINITY;
    for (Disk d : disks) {
      xMin = Math.min(xMin, d.getX() - d.r);
      xMax = Math.max(xMax, d.getX() + d.r);
      yMin = Math.min(yMin, d.getY() - d.r);
      yMax = Math.max(yMax, d.getY() + d.r);
    }
    double w = xMax - xMin, h = yMax - yMin;
    double area = w * h;
    return new double[]{xMin, xMax, yMin, yMax, w, h, area};
  }

  private void randomConfiguration() {
    while (disks.size() < N) {
      double r = rMin + (rMax - rMin) * rng.nextDouble();
      double x, y; boolean good;
      do {
        good = true;
        x = r + (Lx - 2*r) * rng.nextDouble();
        y = r + (Ly - 2*r) * rng.nextDouble();
        for (Disk d : disks) {
          if (Math.hypot(x - d.getX(), y - d.getY()) < r + d.r) {
            good = false;
            break;
          }
        }
      } while (!good);
      disks.add(new Disk(x, y, r));
    }
  }

  private void loadConfiguration(String file) {
    try (BufferedReader br = new BufferedReader(new FileReader(file))) {
      String line;
      while ((line = br.readLine()) != null) {
        line = line.trim();
        if (line.isEmpty() || line.startsWith("#")) continue;
        String[] tok = line.split("[,\\s]+");
        if (tok.length < 3) continue;
        disks.add(new Disk(
          Double.parseDouble(tok[0]),
          Double.parseDouble(tok[1]),
          Double.parseDouble(tok[2])
        ));
      }
    } catch (IOException e) {
      control.println("Could not read " + file + ": " + e);
    }
  }

  @Override
  public void stop() {
    if (writer != null) {
      writer.close();
      writer = null;
    }
  }

  public static void main(String[] args) {
    SimulationControl.createApp(new PackingApp());
  }
}
