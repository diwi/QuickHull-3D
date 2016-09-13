package MAIN_quickhull;

import java.util.Random;

import processing.core.PConstants;
import processing.core.PGraphics3D;


public class PointCloud {
  public static final double PI_TWO = Math.PI*2.0;
  
  float[][] points;
  int num_points;
  Random rand;
  
  public PointCloud(int seed){
    rand = new Random(seed);
  }
  
  
  public static float[] float3(double x, double y, double z){
    return new float[]{(float)x, (float)y, (float)z};
  }
  public static float[] float3(float x, float y, float z){
    return new float[]{x, y, z};
  }
  

  public void randomCube(int num_points, float sx, float sy, float sz){
    points = new float[num_points][3];
    for (int i = 0; i < num_points; i++ ) {
      points[i][0] = (rand.nextFloat()-0.5f)*sx; //(-size,size); 
      points[i][1] = (rand.nextFloat()-0.5f)*sy;
      points[i][2] = (rand.nextFloat()-0.5f)*sz;
    }
    this.num_points = points.length;
  }
  
  public void ifsFile(IFS ifs){
    points =  ifs.v;
    this.num_points = points.length;
  }
  
  public void uniformSampleSphere(int num_points, float sx, float sy, float sz){
    points = new float[num_points][3];
    for (int i = 0; i < num_points; i++ ) {
      points[i] = uniformSampleSphere(sx, sy, sz);
    }
    this.num_points = points.length;
  }
  public void centerSampleSphere(int num_points, float sx, float sy, float sz){
    points = new float[num_points][3];
    for (int i = 0; i < num_points; i++ ) {
      points[i] = centerSampleSphere(sx, sy, sz);
//      points[i][1] *= 1.5f;
    }
    this.num_points = points.length;
  }
  

  private float[] uniformSampleSphere(double sx, double sy, double sz) {
    double phi = rand.nextFloat() * PI_TWO;
    double rnd = rand.nextFloat() * 2.0 - 1.0;
    double rad = Math.sqrt(1.0 - rnd*rnd);
    double X   = sx * Math.cos(phi) * rad;
    double Y   = sy * Math.sin(phi) * rad;
    double Z   = sz * rnd;
    return float3(X,Y,Z);
  }
  
  private float[] centerSampleSphere(double sx, double sy, double sz) {
    sx *= Math.sqrt(rand.nextFloat());
    sy *= Math.sqrt(rand.nextFloat());
    sz *= Math.sqrt(rand.nextFloat());
    double phi = rand.nextFloat() * PI_TWO;
    double rnd = rand.nextFloat() * 2.0 - 1.0;
    double rad = Math.sqrt(1.0 - rnd*rnd);
    double X   = sx * Math.cos(phi) * rad;
    double Y   = sy * Math.sin(phi) * rad;
    double Z   = sz * rnd;
    return float3(X,Y,Z);
  }
  
  
  public void draw(PGraphics3D g){
    g.beginShape(PConstants.POINTS);
      for( int i = 0; i < num_points; i++){
        g.vertex(points[i][0], points[i][1], points[i][2] );
      }
    g.endShape();
  }
  
  


}
