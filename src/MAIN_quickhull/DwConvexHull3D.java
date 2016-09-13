/**
 * 
 *   author: (c) Thomas Diewald, http://thomasdiewald.com
 *   date: 27.02.2013
 *   
 *
 * This source is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This code is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * A copy of the GNU General Public License is available on the World
 * Wide Web at <http://www.gnu.org/copyleft/gpl.html>. You can also
 * obtain it by writing to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



package MAIN_quickhull;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import processing.core.PConstants;
import processing.core.PGraphics3D;

/**
 * 
 * 
 * Creates a 3D Convex Hull of a given list of points.<br>
 * <br>
 * For fast generation, the QuickHull Algorithm is used (32 bit floating 
 * point precision).<br>
 * The Hull internally is represented by the HalfEdge-DataStructure, but can 
 * also be returned as an Indexed-Face-Set (IFS).<br>
 * <br>
 * note: increase memory (-Xms512M -Xmx1024M) when used for a huge number
 * of points (10.000.000 points and more).
 * 
 * @author (c) Thomas Diewald, 2013
 * @web http://thomasdiewald.com
 *
 */
public class DwConvexHull3D {

  /**
   * change this, depending on the extends of the given point-cloud.
   */
  public double EPSILON = 0.00001;

  // good to keep references? memory?
  private ArrayList<float[][]> POINT_CLOUDS = new ArrayList<float[][]>();
  
  private HalfEdge.Mesh mesh = null;  // convex hull mesh
  
  private int NUM_FACES_CVH   = 0; // general counter for the number of mesh-faces
  private int NUM_FACES_LIGHT = 0; // points to last+1 of faces_light
  private int NUM_FACES_NEW   = 0; // points to last+1 of faces_new
  private HalfEdge.Face[] faces_light; // face-buffer, reallocated when needed
  private HalfEdge.Face[] faces_new;   // face-buffer, reallocated when needed
  
  /**
   * dynamic 2d-matrix -> [faces][points as float[3]]<br>
   * <br>
   * each face has its own set of points, that need to be considered during <br>
   * building the hull.
   */
  private final ArrayList< PointSet<float[]> > point_sets = new ArrayList<PointSet<float[]>>();

  /**
   * stack for iterations. contains faces with points.
   * two ways, for creating the hull: breadth-first / depth-first.
   */
  public DwStack STACK = new DwStack();

  private int iterations = 0; // iteration counter, used for debugging
  
  public boolean DEBUG_MODE    = true; // TODO
  public boolean DEBUG_VERBOSE = true; // TODO
  

  /**
   * empty constructor.
   */
  public DwConvexHull3D(){
  }

  

  /**
   * initializes the convex hull using the gven pointcloud.<br>
   * creates the initial pyramid (tetrahedron).<br>
   * assigns the points to its faces, and pushes the faces on the stack.<br>
   * 
   * 
   * @param point_cloud
   * @return null, on failure
   * @throws Exception 
   */
  public DwConvexHull3D cvhInit( float[][] point_cloud ) throws Exception{
    
    if( !validInput(point_cloud) ) return null;

    clearPointClouds();
    POINT_CLOUDS.add(point_cloud);
    
    faces_light = new HalfEdge.Face[100];
    faces_new   = new HalfEdge.Face[100];
    
    long timer;
    
    // init hull
    timer = System.currentTimeMillis();
    float[][] exclude = initMesh(point_cloud);
    timer = System.currentTimeMillis()-timer;
    if(DEBUG_VERBOSE | DEBUG_MODE) System.out.printf("  %-20s (%5d ms)\n", "pointcloud extrema", timer);  
    
    // check if pointcloud is degenerate
    if(mesh == null){
      throw new Exception("QuickHull Error: Pointcloud Degenerate");
    }
    
//    reorderFaces(point_cloud);

    // assign points to new face's point-sets
    timer = System.currentTimeMillis();
    assignPointSet(faces_new, NUM_FACES_CVH, point_cloud, exclude);
    timer = System.currentTimeMillis()-timer;
    if(DEBUG_VERBOSE | DEBUG_MODE) System.out.printf("  %-20s (%5d ms)\n", "assign points", timer);

    // push faces with points on the stack for further subdivisions.
    STACK.reset(500);
    pushFacesOnStack( faces_new, NUM_FACES_CVH);
    
    iterations = 0;
    return this;
  }
  
  
  
  
  void reorderFaces(float[][] point_cloud){
    
    HalfEdge.Face[] faces =  faces_new;
    
    for(int j = 0; j < NUM_FACES_CVH; j++) faces[j].data.tmp_val = 0;
    
    for(int i = 0; i < point_cloud.length; i++){
      for(int j = 0; j < NUM_FACES_CVH; j++){
        HalfEdge.FaceData face_data = faces[j].data;
        float distance = face_data.distance(point_cloud[i]);
        if( distance > EPSILON ){ 
          face_data.tmp_val++;
//          faces[j].data.tmp_val += (int)(distance);
//          if(distance > face_data.tmp_val){
//            face_data.tmp_val = (int)(distance);
//          }
        }
      }
    }
    
//    for(int j = 0; j < NUM_FACES_CVH; j++) System.out.println("face "+j+": "+faces[j].data.tmp_val);
    
    
    Arrays.sort(faces, 0, 4, new Comparator<HalfEdge.Face>() {
      @Override
      public int compare(HalfEdge.Face first, HalfEdge.Face second) {
        return second.data.tmp_val - first.data.tmp_val;
//        return first.data.tmp_val - second.data.tmp_val;
      }
    });
    
    for(int j = 0; j < NUM_FACES_CVH; j++) faces[j].data.tmp_val = 0;
    

    System.out.println("");
    for(int i = 0; i < point_cloud.length; i++){
      for(int j = 0; j < NUM_FACES_CVH; j++){
        HalfEdge.FaceData face_data = faces[j].data;
        float distance = face_data.distance(point_cloud[i]);
        if( distance > EPSILON ){ 
          face_data.tmp_val++;
//          faces[j].data.tmp_val += (int)(distance);
//          if(distance > face_data.tmp_val){
//            face_data.tmp_val = (int)(distance);
//          }
        }
      }
    }
    
//    for(int j = 0; j < NUM_FACES_CVH; j++) System.out.println("face "+j+": "+faces[j].data.tmp_val);
    for(int j = 0; j < NUM_FACES_CVH; j++) faces[j].data.tmp_val = 0;
  }
  
  
  /**
   * extend this hull (finished or not, doesnt't matter) by a given point-cloud.<br>
   * depending on how many faces the hull has at that moment, this may take<br>
   * longer than building a hull from scratch, because each point is tested for<br>
   * all faces, until it is seen by one.<br>
   * 
   * if extending is used a lot, its better to collect points first, and then<br>
   * call cvhExtend().<br>
   * 
   * @param point_cloud
   * @return null, on failure
   */
  public DwConvexHull3D cvhExtend( float[][] point_cloud ){
    if( !validInput(point_cloud) ) 
      return null;
    
    POINT_CLOUDS.add(point_cloud);

    HalfEdge.Face[] faces = HalfEdge.Collect.faces(null, mesh); 
    NUM_FACES_CVH = faces.length;
    
    assignPointSet(faces_new, NUM_FACES_CVH, point_cloud, new float[][]{} );

    STACK.reset(500);
    pushFacesOnStack(faces_new, NUM_FACES_CVH);
    return this;
  }
  

  

  /**
   * 
   * build the convex hull. the number of iterations to is given by "num_iter".
   * if num_iter == -1, the complete hull is generated. otherwise, the given
   * number of iterations are processed.
   * 
   * during calls for cvhBuild, the stack-traversal can be changed by changing LIFO.
   * 
   * @param num_iter number of iterations to do.
   */ 
  public void cvhBuild(int num_iter){
    if(STACK.isEmpty()) return ;
    try {

      long timer = System.currentTimeMillis();

      if( num_iter == -1 ) { 
        while( !STACK.isEmpty() ){
          quickHullIteration();
        }
      } else{
        while( !STACK.isEmpty() && (--num_iter >= 0) ){
          quickHullIteration();
        }
      }
      
      timer = System.currentTimeMillis()-timer;
      if(DEBUG_VERBOSE | DEBUG_MODE) System.out.printf("  %-20s (%5d ms)\n", "building", timer);

    } catch (Exception e) {
      e.printStackTrace();
    }
    if(STACK.isEmpty()) cleanUp(false);
  }
  
  
  
  
  /**
   * creates the initial volume (tetrahedron) and  returns a list of 
   * current hull-vertices (4 points, float[4][3])
   * 
   * @param point_cloud
   * @return current hull-vertices (4 points, float[4][3])
   */
  private float[][] initMesh(float[][] point_cloud){
    
    float[][] v = extremePoints(point_cloud);
    if(v == null){
      mesh = null;
      faces_new = null;
      NUM_FACES_CVH = 0;
      return null;
    }
    
    ArrayList<HalfEdge.Face> faces_tmp = new ArrayList<HalfEdge.Face>(4);
//    mesh = HalfEdge.Creator.newPyramid(faces_tmp, v[0], new float[][]{ v[1], v[2], v[3] });
    mesh = HalfEdge.Creator.newPyramid(faces_tmp, v[3], new float[][]{ v[0], v[1], v[2] });
    
    faces_new     = faces_tmp.toArray(faces_new);
    NUM_FACES_CVH = faces_tmp.size();

    return v;
  }


  /**
   * compute 4 points, for the initial hull (tetrahedron).
   * 
   * @param point_cloud
   * @return [4][3] - 4 points
   */
  private float[][] extremePoints(float[][] point_cloud){
    
    final int NUM_POINTS = point_cloud.length;
    float[] v0, v1, v2, v3; // vertices of initial pyramid (tetrahedron)
    float dmax, dcur;       // tmp max/current distance
    int i, j;               // loop-counter
   
    // 1) extreme points in x, y and z (min and max)
    float[][] ep = {  point_cloud[0], point_cloud[0],    // { [0]x_min, [1]x_max 
                      point_cloud[0], point_cloud[0],    //   [2]y_min, [3]y_max 
                      point_cloud[0], point_cloud[0]  }; //   [4]z_min, [5]z_max }

    for(i = 1; i < NUM_POINTS; i++){
      float[] p = point_cloud[i];                    
      if( p[0] <= ep[0][0] ) { ep[0] = p; /* continue; */ } 
      if( p[0] >= ep[1][0] ) { ep[1] = p; /* continue; */ } 
      if( p[1] <= ep[2][1] ) { ep[2] = p; /* continue; */ } 
      if( p[1] >= ep[3][1] ) { ep[3] = p; /* continue; */ } 
      if( p[2] <= ep[4][2] ) { ep[4] = p; /* continue; */ } 
      if( p[2] >= ep[5][2] ) { ep[5] = p; /* continue; */ } 
    }

    // 2) line segment: find the 2 most distant points.
    for(v0 = v1 = null, dmax = 0, i = 0; i < ep.length; i++){
      for(j = i+1; j < ep.length; j++){
        if( dmax < (dcur = DwVec3.distSq(ep[i], ep[j])) ){
          dmax = dcur;
          v0 = ep[ i ]; 
          v1 = ep[ j ]; 
        }
      }
    }

    if(v0 == null && v1 == null){
      System.out.println("QuickHull Error: pointcloud degenerate");
      return null;
    }
    
    // 3) triangle: find the most distant point to the line segment.
    for(v2 = null, dmax = 0, i = 0; i < ep.length; i++){
      if( dmax < (dcur = DwVec3.distSqPointSegment(v0, v1, ep[i])) ){
        dmax = dcur;
        v2 = ep[i]; 
      }
    }

    if(v2 == null){
      System.out.println("QuickHull Error: pointcloud degenerate (collinear)");
      return null;
    }
  
    // 4) pyramid: find the most distant point to the triangle.
    final float[] N = DwVec3.normal_new(v0, v1, v2);
    final float   D = DwVec3.dot(N, v0);
    
    for(v3 = null, dmax = 0, i = 0; i < NUM_POINTS; i++){
      if( dmax < (dcur = Math.abs((DwVec3.dot(point_cloud[i], N) - D)))){
        dmax = dcur;
        v3 = point_cloud[ i ]; 
      }
    }
    
    if(v3 == null){
      System.out.println("QuickHull Error: pointcloud degenerate (coplanar)");
      return null;
    }
    

    return new float[][]{ v0,v1,v2,v3 };
  }
 

  
  /**
   * main convex-hull building method.
   * 
   * 1) gets a face from the stack.
   * 2) finds its max distant point.
   * 3) determines all adjacent faces, that are seen from that point (=light-faces)
   * 4) extracts a horizon-ring of all that faces. ( this detaches light-faces )
   * 5) creates new faces from the horizon-edges to the point (pyramid: point=apex, horizon=base)
   * 6) assigns the points of the light-faces to the new created faces.
   * 7) pushes faces, with points on the stack.
   * 
   * @throws Exception (on really bad errors)
   */
  private void quickHullIteration() throws Exception{
    
    // (0) when using FIFO some face's data might be null now -> skip them.
    // !!! only for demo, to avoid any errors, when switching between
    // LIFO/FIFO during building !!!
    while(!STACK.isEmpty() && STACK.peek().data == null){
      STACK.pop();
    }
    if( STACK.isEmpty() ){
      return;
    }
    
    if( DEBUG_VERBOSE ) {
      System.out.println("\n\n");
      System.out.println("___________________              ____________________");
      System.out.println("___________________NEXT ITERATION____________________");
      System.out.println("___________________              ____________________");
      System.out.println("");
      System.out.println("# "+iterations);
    }
    
    iterations++; // for debugging
    
    

    // (1)
    HalfEdge.Face face_cur = STACK.pop();

    // (2)
    // getting max distant point from current face.
    float[] POINT_MAX = getPointMax(face_cur);

    // (3)
    // collecting light-faces (faces, seen from POINT_MAX), starting at current 
    // face (which certainly is a light-face) since POINT_MAX is from its point-set.
    HalfEdge.Edge hor_start = getLightFaces(face_cur, POINT_MAX);

    // (4) and (5)
    // find horizon and create new pyramid using the horizon as the base and 
    // POINT_MAX as the apex. the old (light) faces are detached from the mesh.
    getNewFaces(hor_start, POINT_MAX);
    
    // (6)
    // assign all points of all light-faces to the new created face's point-sets.
    // points inside the hull-volume and on the hull-surface are ignored.
    assignPointSets(POINT_MAX);

    // (7)
    // push the new faces (in case they contain any points now) on the stack.
    pushFacesOnStack(faces_new, NUM_FACES_NEW);
    
    
    // (8) when using FIFO some face's data might be null now -> skip them.
    while(!STACK.isEmpty() && STACK.peek().data == null){
      STACK.pop();
    }
    
    // keep number of mesh-faces up to date
    NUM_FACES_CVH -= NUM_FACES_LIGHT; // removed light faces
    NUM_FACES_CVH += NUM_FACES_NEW;   // but added new faces
  }
  
  
  
  
  /**
   * get max distant point from a face, checking the face's point-set.
   * @param face
   * @return
   */
  private final float[] getPointMax(final HalfEdge.Face face){
    HalfEdge.FaceData data = face.data;
//    if(data == null ){
//      System.out.println("data is null" +data);
//      System.out.println("STACK.LIFO : "+STACK.LIFO);
//    }
    final PointSet<float[]> points = point_sets.get(data.tmp_val_2);
//    if(points == null ){
//      System.out.println("points is null" +data.tmp_val_2);
//    }
    final int num_points = points.size();
    
    float[] cur_pnt, max_pnt = points.get(0);
    float   cur_dis, max_dis = data.distance( max_pnt );
    
    for( int i = 1; i < num_points; i++){
      cur_pnt = points.get(i);
      cur_dis = data.distance( cur_pnt );
      if( max_dis < cur_dis){
        max_dis = cur_dis;
        max_pnt = cur_pnt;
      }
    }
    return max_pnt;
  }
  
  
  
  
  
  
  private HalfEdge.Edge getLightFaces(HalfEdge.Face face_cur, float[] POINT_MAX){
    
    NUM_FACES_LIGHT = 0;
    HalfEdge.Edge hor_start = null;
    
    if( faces_light.length <= NUM_FACES_CVH){
      faces_light = new HalfEdge.Face[ (int)(NUM_FACES_CVH*1.75) ]; // increase buffer
    }
    
    final int VISITED = ++mesh.flag.visited; // current label-id for visited faces.
    int stack_ptr = 0;
    faces_light[stack_ptr++] = face_cur; // faces_light is used as stack and as the list of light-faces.
    face_cur.flag.visited = VISITED;
    face_cur.data.tmp_val = 1; // 1 = "in light", 0 = "not in light"

    // find adjacent light faces.
    while( stack_ptr > NUM_FACES_LIGHT ){
      
      face_cur = faces_light[NUM_FACES_LIGHT++];

      // find adjacent faces
      HalfEdge.Edge edge = face_cur.edge;
      do {   

        // WARNING (TODO)
        // due to precision-errors, some "non-light-faces" may be surrounded by "light-faces"
        // which clearly can't be. this is taken care of when building the horizon ( at least a bit).
        // in such a case, those faces ( and assigned points!!!) are removed from
        // the current convex-hull-mesh.
        // this introduces some problems: 
        //    points that need to be computed are lost.
        //    number of faces differ from NUM_FACES_CVH
        //    horizon during current iteration contains same vertices one or more times.
        //    horizon might not be "convex" (seen from the point 2 plane vector)
        // anyway, the hull continues without any errors, but containing less face.
        // edit: when using FIFI, problem is obsolete
        
        if(  edge.pair == null ){
          throw new NullPointerException("QUICKHULL ERROR, while getting faces: edge.pair == null." );
        }
        HalfEdge.Face adj_face = edge.pair.face;
        if( adj_face.flag == null ){
          throw new NullPointerException("QUICKHULL ERROR, while getting faces: edge.pair.face.flag == null" );
        }
        // if the current face has not been visited yet, check its distance, and
        // push it on the "faces_light"-stack.
        if( adj_face.flag.visited != VISITED ){  
          adj_face.flag.visited = VISITED;
          adj_face.data.tmp_val = (adj_face.data.distance(POINT_MAX) > EPSILON) ? 1 : 0;
          if( adj_face.data.tmp_val == 1){
            faces_light[stack_ptr++] = adj_face;
          }
        } 

        // save horizon-edge if the adjacent faces differ in their label (in light/not in light)
        if( (hor_start == null) && (adj_face.data.tmp_val ^ face_cur.data.tmp_val ) == 1){
          hor_start = edge;
        }
        
      } while( (edge=edge.next) != face_cur.edge);
    }

    return hor_start;
  }
  
  
  

  
  
  private void getNewFaces(HalfEdge.Edge hor_start, float[] POINT_MAX){
    NUM_FACES_NEW = 0;
    {
      // check/swap direction, to assure CCW horizon.
      if( hor_start.face.data.tmp_val == 0 ) 
        hor_start = hor_start.pair;

      hor_start.orig.edge   = hor_start;
      hor_start.face.edge   = hor_start;
      HalfEdge.Edge hor_cur = hor_start;

      // in general all light-faces are adjacent to each other and form
      // a convex horizon (when seen from POINT_MAX). 
      do{
        hor_cur.orig.edge = hor_cur;
        hor_cur.face      = hor_start.face;

        while( (hor_cur.next.pair).face.data.tmp_val == 1 ){
          hor_cur.next = hor_cur.next.pair.next;
        }

        NUM_FACES_NEW++;
      } while ((hor_cur = hor_cur.next) != hor_start);

      // WARNING (TODO)
      // the following is somehow an ugly hack.
      // due to imprecision not all light-faces were recognized, and so
      // the horizon may contains some vertices multiple times -> self intersecting.
      // the following loop, tries repair the horizon, by only using the 
      // outermost edges. of course, some faces inside the new ring, were
      // never labeled as light faces, and therefore, their point-sets are lost!!!
      // this hack, only applies in rare cases: 
      // given a VERY BIG point-cloud with a VERY LITTLE x/y/z extent.
      // though, the problem of having an invalid convex hull remains.
      // edit: when using FIFI, problem is obsolete
      
//      if( STACK.LIFO )
      {
        // ~~~~~~~ HACK START ~~~~~~~
        int num_faces_after_hack = 0;
        do{
          hor_cur.face.edge = hor_cur;
          hor_cur.next = hor_cur.next.orig.edge;
          if( ++num_faces_after_hack > NUM_FACES_NEW ) {
//            System.out.println("ERROR");
            num_faces_after_hack = 1;
            hor_start = hor_cur;
            hor_start.next.face.edge = hor_start;
            hor_cur.face.edge = hor_cur;
            hor_cur = hor_cur.next;
            continue;
          }
  
        } while ((hor_cur = hor_cur.next) != hor_start);
  //      hor_cur.face.edge = hor_cur;
        
        if( NUM_FACES_NEW != num_faces_after_hack ){
          if(DEBUG_VERBOSE | DEBUG_MODE) 
            System.err.println("ERROR: bad horizon");
          NUM_FACES_NEW = num_faces_after_hack;  // number of new faces, may be lower now.
        }
        // ~~~~~~~ HACK END ~~~~~~~
      }
      
      if( faces_new.length < NUM_FACES_NEW ){
        faces_new = new HalfEdge.Face[NUM_FACES_NEW]; 
      }
      
      // hor_start.face is the face (polygon, not planar) of the current horizon.
      NUM_FACES_NEW = HalfEdge.MeshOps.triangulateFace(mesh, hor_start.face, POINT_MAX, faces_new);
//      int NUM_FACES_NEW_tmp = HalfEdge.MeshOps.triangulateFace(mesh, hor_start.face, POINT_MAX, faces_new);
//      if( NUM_FACES_NEW_tmp != NUM_FACES_NEW){
//        System.out.println("__"+NUM_FACES_NEW_tmp+", "+NUM_FACES_NEW);
//      }
    }
  }
 


  
  
  
  
  private void assignPointSets(float[] POINT_MAX){
    if( DEBUG_VERBOSE ){
      System.out.println("\n___________________CLEARING LIGHT FACES ("+NUM_FACES_LIGHT+")____________________");
    }
    // assign all points of all light-faces to the new created face's point-sets.
    // points inside the hull-volume and on the hull-surface are ignored.
    for( int k = 0; k < NUM_FACES_LIGHT; k++)
    {
      HalfEdge.Face L_FACE = faces_light[k];
      HalfEdge.Edge L_EDGE = L_FACE.edge;
      
      // make sure, that nothing else is linked to the light-faces (being removed)
      // doesn't look that nice, and is a TODO
      if( L_EDGE != null ){
        // if light faces edge-face, is pointing to the light-face itself, every connection can be deleted
        if( L_EDGE.face == L_FACE ){
          if( L_EDGE.orig.edge == L_EDGE ){
            L_EDGE.orig.edge = null;
          }
          if( L_EDGE.pair != null ){
            if( L_EDGE.pair.pair == L_EDGE){
              L_FACE.edge.pair.pair = null;
              L_FACE.edge.pair = null;
            }
          }
          L_EDGE.face = null;
          L_EDGE.next = null;
          L_EDGE.orig = null; 
          L_EDGE.pair = null; 
          L_EDGE.data = null;
          L_EDGE.flag = null;
        }
        L_FACE.edge = null;
      }

      // hor_start.face is -1 by default! (to throw an error if used in the wrong place)
      int point_set_idx = L_FACE.data.tmp_val_2;
      int num_points = 0;
      if( point_set_idx != -1 ){
        // get point-set of current light-face
        final PointSet<float[]> point_set = point_sets.get(point_set_idx);
        num_points = point_set.size();
        // assign to the new faces
        assignPointSet(faces_new, NUM_FACES_NEW, point_set, POINT_MAX);
        // delete the old set. this also improves performance.
        point_set.clear();
        point_sets.set(point_set_idx, null);
      }
      if( DEBUG_VERBOSE ){
        System.out.printf("   light-face:  stack[%2d]  num_pts=%5d (point_sets[%3d])\n", L_FACE.data.id, num_points, point_set_idx);
      }
      L_FACE.data = null;
      L_FACE.edge = null;
      L_FACE.flag = null;
    }

  }
  
  
  
  /**
   * used at the beginning, for assigning the initial point-cloud to the initial faces.
   * 
   * @param faces          new created faces, not added to the stack yet
   * @param num_faces      number of new created faces. (!= faces.length)
   * @param points         new points (from initial point-cloud)
   * @param exlude_points  points, that are not considered. ( because they are on the hull, or whatever)
   */
  private void assignPointSet(HalfEdge.Face[] faces, int num_faces, float[][] points, float[] ... exlude_points ){
    int num_points = points.length;
    int est_points = Math.max(num_points/(num_faces), 1);
//    est_points = Math.min(32, est_points);

    for( int i = 0; i < num_points; i++){
      assignPoint(faces, num_faces, points[i], est_points, exlude_points);
    }
  }
  /**
   * 
   * @param faces           new created faces, not added to the stack yet                                
   * @param num_faces       number of new created faces. (!= faces.length)                               
   * @param point_set       new points, from other faces point-sets                                                              
   * @param exlude_points   points, that are not considered. ( because they are on the hull, or whatever)
   */
  private void assignPointSet(HalfEdge.Face[] faces, int num_faces, PointSet<float[]> point_set, float[] ... exlude_points ){
    int num_points = point_set.size();
    int est_points = Math.max(num_points/(num_faces), 1);
    est_points = Math.min(16, est_points);

//    float[][] points = point_set.list;
    for( int i = 0; i < num_points; i++){
//      assignPoint(faces, num_faces, points[i], exlude_points);
      assignPoint(faces, num_faces, point_set.get(i), est_points, exlude_points);
    }
  }


  // some testing for assigning points.
  // in my original version, each point is assigned to the first face in the list.
  // here i made some, test, to assign the point to the face, which it is most distant/close to.
//  private void assignPointSet_UNUSED(HalfEdge.Face[] faces, int num_faces, ArrayList<float[]> points, float[] ... exlude_points ){
//    int num_points = points.size();
//    
//    __NEXT_POINT__:
//    for( int i = 0; i < num_points; i++){
//      float[] point = points.get(i);
//      for(int j = 0; j < exlude_points.length; j++){
//        if( point == exlude_points[j]) 
//          continue __NEXT_POINT__;
//      }
//      
//      // save point to first face it is seen from, and return.
//      HalfEdge.FaceData data = null;
////      float dis_max = Float.MAX_VALUE;
//      float dis_max = (float) EPSILON;
//      for(int j = 0; j < num_faces; j++ ){
//        HalfEdge.FaceData data_cur = faces[j].data;
//        float dis_cur = data_cur.distance(point);
//        if( dis_cur > EPSILON && dis_cur > dis_max ){ 
//          dis_max = dis_cur;
//          data = data_cur;
//        }
//      }
//
//      if( data != null){
//        if( data.tmp_val_2 == -1 ) {
//          data.tmp_val_2 = point_sets.size();
//          point_sets.add( new ArrayList<float[]>() );
//        } 
//        point_sets.get(data.tmp_val_2).add(point);
//      }
//      
//    }
//  }
  
  

  
  /**
   * 
   * @param faces           new created faces, not added to the stack yet
   * @param num_faces       number of new created faces. (!= faces.length) 
   * @param point           point, that is tested for being seen by the faces.
   * @param est_num_pts     number of estimated points for this face, to do a nice allocation of the point-set.
   * @param exlude_points   points, that are not considered. ( because they are on the hull, or whatever)
   */
  private void assignPoint(HalfEdge.Face[] faces, int num_faces, float[] point, int est_num_pts, float[] ... exlude_points){
    // skip hull-points
    for(int j = 0; j < exlude_points.length; j++)
      if( point == exlude_points[j]) 
        return;
    
    // save point to first face it is seen from, and return.
    for(int j = 0; j < num_faces; j++ ){
      HalfEdge.FaceData data = faces[j].data;
      if( data.distance(point) > EPSILON ){ 
        PointSet<float[]> ps;
        if( data.tmp_val_2 == -1 ) {
          data.tmp_val_2 = point_sets.size(); // assign index to point-set
          point_sets.add( ps = new PointSet<float[]>(est_num_pts) );
        } else {
          ps = point_sets.get(data.tmp_val_2);
        }
        ps.add(point);
        return;
      }
    }
  }

  

  
  

  
  /**
   * pushes new created faces (in case they contain points) on the stack.
   * 
   * @param faces
   * @param num_faces
   */
  private final void pushFacesOnStack( final HalfEdge.Face[] faces, final int num_faces){
    
    if( DEBUG_VERBOSE ){
      System.out.println("\n___________________PUSH NEW FACES ON STACK____________________");
    }
    
    STACK.assureSize(num_faces);
    
    
    for( int i = 0; i < num_faces; i++){
      if( faces[i].data.tmp_val_2 != -1){
        faces[i].data.id = STACK.STACK_PTR_LAST; // just for debugging
        STACK.push(faces[i]);
              
        if( DEBUG_VERBOSE ){
          int point_set_idx = faces[i].data.tmp_val_2;
          int num_points = point_sets.get(point_set_idx).size();
          System.out.printf("   new face:  stack[%2d]  num_pts=%5d (point_sets[%3d])\n", faces[i].data.id, num_points, point_set_idx);
        }
        
        // insertion-sort, to have face sorted by their number of assigned points
//        HalfEdge.Face face_new = faces[i];
//        int face_new_num_pts = point_sets.get(face_new.data.tmp_val_2).size();
//        int j;
//        for(j = STACK_PTR_LAST; j > STACK_PTR_FIRST; j--){
//          HalfEdge.Face face_cur = STACK[j-1];
//          if( face_cur.data == null){
//            break;
////            STACK[j] = face_cur;
//          }
//                   
//          int face_cur_num_pts = point_sets.get(face_cur.data.tmp_val_2).size();
//          
//          if( face_cur_num_pts > face_new_num_pts ){
//            STACK[j] = face_cur;
//          } else {
////            STACK[j] = face_new;
//            break;
//          }
//        }
//        STACK[j] = face_new;
//        STACK_PTR_LAST++; 
        
      }
    }
    
    // print current stack, for debugging
//    System.out.println("");
//    for( int i = STACK_PTR_FIRST; i < STACK_PTR_LAST; i++){
//      HalfEdge.Face face = STACK[i];
//      if( face.data == null ) {
//        System.out.println("stack["+i+"] num_points = NULL");
//      } else {
//        ArrayList<float[]> points = point_sets.get(face.data.tmp_val_2);
//        System.out.println("stack["+i+"] num_points = "+points.size());
//      } 
//    }
    
    
  }
  
  
  
  private boolean validInput(float[][] point_cloud){
    if( point_cloud == null ){
      System.err.println( "point-array (float[][]) is null");
      return false;
    }
    if( point_cloud.length < 4 ){
      System.err.println( "point-array (float[][]) has length < 4");
      return false;
    }
    return true;
  }
 
  
  /**
   *  just in case.
   */
  private void cleanUp(boolean clear_pointclouds){
    faces_light = null;
    faces_new   = null;
    point_sets.clear();
    point_sets.trimToSize();
    STACK.reset(0);
    //TODO: reset flags too?
    if( clear_pointclouds )
      clearPointClouds();
  }
  
  public void clearPointClouds(){
    for( int i = 0; i < POINT_CLOUDS.size(); i++) 
      POINT_CLOUDS.set(i, null);
    POINT_CLOUDS.clear();
//    POINT_CLOUDS = null; // ?
  }
  
  public int numPoints(){
    int num_points = 0;
    for( float[][] pc : POINT_CLOUDS ){
      num_points += pc.length;
    }
    return num_points;
  }
  
  public int numFaces(){
    return NUM_FACES_CVH;
  }

  public int numFacesExact(){
    return (NUM_FACES_CVH = HalfEdge.Collect.faces(null, mesh).length);
  }
  
  public int numIterations(){
    return iterations;
  }
  
  public void updateMeshData(){
    HalfEdge.Update.verts(mesh); // always do this before rendering, may takes some time
  }
  public HalfEdge.Mesh getMesh(){
    return mesh;
  }
  
  public boolean finished(){
    return STACK.isEmpty();
  }
  
  //TODO
  public IFS getMeshAsIFS(){
    System.out.println("getMeshAsIFS() is not implemented yet");
    return null;
  }
  

  
  
  
  
  
  
  
  
  
  
  
  
  /**
   * 
   * setting for stack(iteration)-behavior. (IMPORTANT)<br>
   * <br>
   * LIFO = depth first   (LIFO = true )<br>
   * FIFO = breadth first (LIFO = false)<br>
   * <br>
   * FIFO may perform better in general (less iterations,...), but worse in a <br>
   * worst case, when the convex hull contains almost all points of the given <br>
   * point-cloud.<br>
   * e.g. a geosphere as a given point-cloud -> every vertex will be a hull-vertex.<br>
   * <br>
   * FIFO produces a bigger stack and therefore has a higher memory-usage.<br>
   * 
   * <br>
   * default: "LIFO = true" (might change?)<br>
   * 
   */
  public void LIFO(boolean LIFO){
    STACK.LIFO = LIFO;
  }
  public boolean LIFO(){
    return STACK.LIFO;
  }
  
  
  
  
  /**
   * Stack for quick-hull iterations.
   * contains faces with point-sets.
   * traversal: breadth-first/depth-first
   * 
   * @author Thomas Diewald
   *
   */
  private static class DwStack{
    
    private boolean LIFO = true;
    private int STACK_PTR_LAST  = 0;
    private int STACK_PTR_FIRST = 0;
    private HalfEdge.Face[] STACK = new HalfEdge.Face[500];
    //private final ArrayDeque<HalfEdge.Face> STACK = new ArrayDeque<HalfEdge.Face>();
    
    public void reset(int stack_size){
      STACK_PTR_LAST  = 0;
      STACK_PTR_FIRST = 0;
      STACK = new HalfEdge.Face[stack_size];
      //  STACK.clear();
    }
    
    public boolean isEmpty(){
      return (STACK_PTR_FIRST >= STACK_PTR_LAST);
      // return (STACK.isEmpty() );
    }
    
    public void assureSize(int extend){
      if( STACK.length < STACK_PTR_LAST+extend){
//        int new_size = STACK_PTR_LAST+extend;
//        HalfEdge.Face[] STACK_tmp = new HalfEdge.Face[(int)(new_size*1.5f)];
//        System.arraycopy(STACK, 0, STACK_tmp, 0, STACK_PTR_LAST);
//        STACK = STACK_tmp;
        
        STACK = Arrays.copyOf(STACK, (STACK_PTR_LAST+extend)<<1);
  //      System.out.println(stack_count++ +", "+new_size);
      }
    }
    
    public void push(HalfEdge.Face face){
      STACK[STACK_PTR_LAST++] = face;
//      STACK.addLast(face); //STACK.addFirst(face);
    }
    // FIFO performs noticeable better on a more uniform distribution of points.
    // LIFO performs way better on a distribution on for example a sphere surface,
    // where each point is part of the convex hull. in that case, the stack gets really big.
    public HalfEdge.Face pop(){
      return STACK[ LIFO ? --STACK_PTR_LAST : STACK_PTR_FIRST++ ];
//      return LIFO ? STACK.pollLast()/*LIFO*/ : STACK.pollFirst()/*FIFO*/;
 
    }
    public HalfEdge.Face peek(){
      return STACK[ LIFO ? STACK_PTR_LAST-1 : STACK_PTR_FIRST ];
    }
  }
  
  
  
  

  /**
   * dynamic list of points (float[3])
   * could also use an ArrayList, but for testing/debugging, this is more comfortable.
   * 
   * @author thomas
   *
   * @param <T> point-type (e.g. float[3]).
   */
  private static class PointSet<T>{
  //  static int count_reallocs = 0;
  //  final ArrayList< T> list = new ArrayList< T>();
  
    private int end, size;
    private T[] list;
  
    @SuppressWarnings("unchecked")
    public PointSet(int size){
      this.list = (T[]) new Object[this.size=size];
    }
    
    public T get(int pos){
      return list[pos];
  //    return list.get(pos);
    }
    
    public void add(T point){
      if( end >= size) reAlloc( size<<1 );
      list[end++] = point;
  //    list.add(point);
    }
    public void shorten(int new_end){
      end = new_end;
    }
//    public void trimToSize(){
//      reAlloc(end);
//    }
    @SuppressWarnings("unchecked")
    public void clear(){
      shorten(0);
      end = size = 0;
      this.list = (T[]) new Object[0];
    }
    private final void reAlloc(int new_size){
  //    count_reallocs++;
  //    if( count_reallocs %50 == 0) System.out.println(count_reallocs);
//      @SuppressWarnings("unchecked")
//      T[] list_new = (T[]) new Object[size = new_size];
//      end = Math.min(end, size);
//      System.arraycopy(list, 0, list_new, 0, end);
//      list = list_new;
      list = Arrays.copyOf(list, size = new_size);
    }
    public int size(){
      return end;
  //    return list.size();
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  // for debugging only.
  // this method can savely be removed when used within another rendering context 
  // than processing.
  public void draw(PGraphics3D g){
    if( mesh != null ){

      // draw convex hull
      g.noStroke();        g.fill(150);          HalfEdge_Render.faces(g, mesh);
      g.stroke(10);        g.strokeWeight(1.6f); HalfEdge_Render.edges(g, mesh);
      g.stroke(255,200,0); g.strokeWeight(6);    HalfEdge_Render.verts(g, mesh);
//      g.strokeWeight(2); HalfEdge_Render.faceNormals(g, mesh, 150);
////    g.strokeWeight(1); HalfEdge_Render.vertNormals(g, mesh, 150);
//      g.strokeWeight(3); HalfEdge_Render.halfEdges(g, mesh);
//      g.strokeWeight(2); HalfEdge_Render.halfEdgesPairing(g, mesh);
      

      // for demo/debugging, simulate next iteration.
      try {
        NEXT_ITERATION:
        {
          // STUFF need for drawing
          HalfEdge.Face face_cur;
          float[] POINT_MAX;
          ArrayList<HalfEdge.Edge> horizon = new ArrayList<HalfEdge.Edge>();
          
          // simulate next iteration
          {
            if( STACK.isEmpty() )
              break NEXT_ITERATION;

            face_cur = STACK.peek();
   
            POINT_MAX = getPointMax(face_cur);
              
            // get horizon edges:
            HalfEdge.Edge hor_start = getLightFaces(face_cur, POINT_MAX);
            // check/swap direction, to assure CCW horizon.
            if( hor_start.face.data.tmp_val == 0 ) 
              hor_start = hor_start.pair;
    
            HalfEdge.Edge hor_next, hor_cur = hor_start;
            do{
              horizon.add(hor_cur);
              hor_next = hor_cur.next;
              while( (hor_next.pair).face.data.tmp_val == 1 ){
                hor_next = hor_next.pair.next;
              }
            } while ((hor_cur = hor_next) != hor_start);
          } // end simulation

          
          
          // draw all light faces
          g.strokeWeight(1);
          g.noStroke();
          for( int i = 0; i < NUM_FACES_LIGHT; i++){
            if( faces_light[i] == face_cur) continue;
            g.fill(125,0,100);
            HalfEdge_Render.face(g, faces_light[i]);
          }
          
          // draw FIRST light-face
          g.fill(255,0,200);
          HalfEdge_Render.face(g, face_cur);
  
          // draw horizon
          g.noFill();
          g.stroke(0,180,200);
          g.strokeWeight(5);
          g.beginShape();
          for( HalfEdge.Edge edge : horizon ){
            vertex(g, edge.orig.data.v);
          }
          g.endShape(PConstants.CLOSE);
          
          // draw pyramid: horizon vertices -> point
          g.stroke(0,180,200);
          g.strokeWeight(2);
          g.beginShape(PConstants.LINES);
          for( HalfEdge.Edge edge : horizon ){
            vertex(g, POINT_MAX);
            vertex(g, edge.orig.data.v);
          }
          g.endShape();
          
          // draw currents face point-set
          g.stroke(50,50,175,200);
          g.strokeWeight(4);
          int point_set_idx = face_cur.data.tmp_val_2;
          if( point_set_idx != -1 ){
            PointSet<float[]> pointset = point_sets.get(point_set_idx);
            for(int i = 0; i < pointset.size(); i++){
              point(g, pointset.get(i));
            }
          }

          // draw POINT_MAX of FIRST light-face
          g.strokeWeight(10);
          g.stroke(235,0,190);
          point(g, POINT_MAX);
          
        } // end NEXT_ITERATION
      } catch(NullPointerException e){
//        System.errr.println("(draw) ERROR during draw");
      }
        
    } // if( mesh != null ){
      
    // draw (all) pointclouds
    if( POINT_CLOUDS != null ){
      for( float[][] pc : POINT_CLOUDS ){
        drawPointCloud(g, pc);
      }   
    }
    
    
  }
  
  
  private void point(PGraphics3D g, float[] v){
    g.point(v[0], v[1], v[2]);
  }
  
  private void vertex(PGraphics3D g, float[] v){
    g.vertex(v[0], v[1], v[2]);
  }

  private void drawPointCloud(PGraphics3D g, float[][] points){
    int num_points_ptcloud = points.length;
    int num_points_to_draw = num_points_ptcloud;
    int skip = Math.max(1, num_points_ptcloud/num_points_to_draw);
    
    if( num_points_ptcloud > 10000){
      num_points_to_draw = 20000;
      skip = Math.max(1, num_points_ptcloud/num_points_to_draw);
    }
    g.stroke(0,60,125,100); 
    g.strokeWeight(3); 
    g.beginShape(PConstants.POINTS);
    for(int i = 0; i < points.length; i+=skip){
      float[] p = points[i];
      g.vertex(p[0], p[1], p[2]);
    }
    g.endShape();
  }

}
