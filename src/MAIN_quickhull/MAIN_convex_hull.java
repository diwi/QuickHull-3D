package MAIN_quickhull;



import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Locale;
import peasy.*;
import processing.core.PApplet;
import processing.core.PGraphics3D;


public class MAIN_convex_hull extends PApplet {
  
//----------------------------------------------------------------------------
//
// author: (c) thomas diewald
// date: 17.01.2013
// shadow mapping
//----------------------------------------------------------------------------

  DwConvexHull3D conv_hull;
  
  int num_points = 1024*256;
  
  
  float[][] parsePointListFile(String filename) throws IOException{
      BufferedReader br = new BufferedReader(new FileReader(filename));
  
      float[][] points = null;
      
      try {
          String line;
          int idx = 0;
          while ((line = br.readLine()) != null) {
              line = line.trim();
              String[] tokens = line.split(" ");
              if(tokens.length == 1){
                num_points = Integer.parseInt(tokens[0]);
                points = new float[num_points][3];
                System.out.println("num_points = "+num_points);
              }
              if(tokens.length == 3){
                 points[idx][0] = Float.parseFloat(tokens[0]);
                 points[idx][1] = Float.parseFloat(tokens[1]);
                 points[idx][2] = Float.parseFloat(tokens[2]);
                 idx++;
                 
               
              }
          }
          if(idx != num_points){
            return null;
          }
          System.out.println("idx = "+idx);
      } finally {
          br.close();
      }
      
      
      return points;
  }
  
  
  

  public void setup(){
    size(800, 800, P3D);
    
    PeasyCam cam = new PeasyCam(this, 0, 0, 0, 1500);
    cam.setRotations( -1.3968815, 0.71367186, -0.12572986);
    
    
    PointCloud pc = createPointCloud(2, num_points);
    float[][] points = pc.points;
    
    
//    try {
////      points = parsePointListFile("data/uniform_points.8192.in");
//      points = parsePointListFile("data/uniform_points.4194304.in");
//    } catch (IOException e) {
//      // TODO Auto-generated catch block
//      e.printStackTrace();
//    }
    
    
    
    
//    points = new float[][]{
////        { 1, -1, 2},
////        {-1,  1, 2},
////        {-1, -1, 2},
////        {-1,  0, 2},
////        { 0,  0, 2},
////        { 0, -1, 2},
////        { 0,  1, 2},
////        { 1,  1, 2},
////        { 1,  0, 2},
//        { 0, 0, 0},
//        { 1, 0, 0},
//        { 1, 0, 0},
//        { 1, 0, 0},
//
//    };
    
    



    try {
      long timer = System.currentTimeMillis();
      
      conv_hull = new DwConvexHull3D();
      
      conv_hull.DEBUG_MODE    = false; 
      conv_hull.DEBUG_VERBOSE = false;
      
      conv_hull.LIFO(!true );
      
      conv_hull.cvhInit(points);
      conv_hull.cvhBuild(-1);

      timer = System.currentTimeMillis()-timer;
      System.out.println("");
      System.out.println("  total time  "+timer+" ms");      
      System.out.println("  iterations  "+conv_hull.numIterations());
      System.out.println("  num faces   "+conv_hull.numFaces() +" (exact="+conv_hull.numFacesExact()+")");
      System.out.println("  LIFO        "+conv_hull.LIFO());  
     
    } catch (Exception e) {
      e.printStackTrace();
      return;
    }

//    
//    conv_hull.DEBUG_MODE    = false; 
//    conv_hull.DEBUG_VERBOSE = false;
    
//   Verdana
    textFont( createFont("Calibri", 12));
//    textFont( createFont("CenturyGothic", 12));
//    textFont( createFont("LucidaSans", 12));
//    textFont( createFont("Simplex", 12));
//    textFont( createFont("Verdana", 11));
    textMode( SCREEN );
  }
  
  
  PointCloud createPointCloud(int type, int num_points){
//    long timer = System.currentTimeMillis();

    PointCloud point_cloud = new PointCloud(200); 
    
    switch(type){
      case 1: point_cloud.randomCube         ( num_points, 1000, 1000, 1000); break;
      case 2: point_cloud.centerSampleSphere ( num_points,  500,  500,  500); break;
      case 3: point_cloud.uniformSampleSphere( num_points,  500,  500,  500); break; // check: F = 2*V-4 (number of faces = 2*number of vertices - 4)
    }
    
//    timer = System.currentTimeMillis()-timer;
//    System.out.println("created pointcloud "+timer+" ms, points = "+point_cloud.num_points);
    return point_cloud;
  }
  

  @Override
  public void draw() {
    background(255);
    lights();

    conv_hull.draw( (PGraphics3D) this.g );

    updateInfoText();
    updateWindowTitle();
  }
  
  void updateInfoText(){

    int dy = 14;
    int x = 20;
    int y = 20;

    fill(0);
    
    String txt_what       = "Convex Hull 3D: Quickhull";
    
    String txt_num_points = "num points: "+conv_hull.numPoints();
    String txt_num_faces  = "num faces: "+conv_hull.numFaces();
    String txt_num_iter   = "num iterations: "+conv_hull.numIterations();
  
    String txt_iterations = "iterations/step [ UP, DOWN ]: "+iterations;
    String txt_stack      = "stack [ x,y ]: "+ (conv_hull.LIFO() ?  "depth first (LIFO)" : "breadth first (FIFO)");

    String txt_finished   = "finish [ f ]: "+conv_hull.finished();
    String txt_restart    = "restart [1,2,3]";
    

    text(txt_what      , x, y); y+= dy;
    y+= dy;
    text(txt_num_points, x, y); y+= dy;
    text(txt_num_faces , x, y); y+= dy;
    text(txt_num_iter  , x, y); y+= dy;
    y+= dy;
    text(txt_iterations, x, y); y+= dy;
    text(txt_stack     , x, y); y+= dy;
    y+= dy;
    text(txt_finished  , x, y); y+= dy;
    text(txt_restart   , x, y); y+= dy;
  }

  void updateWindowTitle(){
    String txt_title     = "QuickHull";
    String txt_framerate = String.format(Locale.ENGLISH, "fps %6.2f", frameRate);
    txt_title += " | " + txt_framerate;
    frame.setTitle(txt_title);
  }
  
  
  @Override
  public void keyPressed(){
    if( key == ESC) key = 0;
    
    if(  key == 'x' || key == 'y' ){  
      boolean LIFO = (key == 'x');
      conv_hull.LIFO( LIFO );
      conv_hull.cvhBuild(iterations);
    }
    
    if( key == CODED){
      if( keyCode == UP   ) iterations++;
      if( keyCode == DOWN ) iterations--;
      iterations = Math.max(iterations, 1);
    }
  }
  
  int iterations = 1;
  

  int screen_shot_nr = 0;
  @Override
  public void keyReleased(){
    if (key == 's' ){
      String stack = ""+ (conv_hull.LIFO() ?  "LIFO" : "FIFO");
      String iter  = ""+ conv_hull.numIterations();
      String nr    = String.format("%02d", screen_shot_nr++);
      String filename = "data/screenshot/diewald_quickhull_"+iter+"_"+stack+"_"+nr+".png";
      save(filename);
      System.out.println("saved image: "+filename);
    }
    

    if( key == '1' || key == '2' || key == '3' ){
      try {
        conv_hull.cvhInit( createPointCloud( key-'0', num_points ).points );
      } catch (Exception e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    
    if( key == 'f'){
      conv_hull.cvhBuild(-1);
    }
  }
  
  

  public static void main(final String args[]) {
    PApplet.main(new String[] { MAIN_convex_hull.class.getName() });
  }
}

