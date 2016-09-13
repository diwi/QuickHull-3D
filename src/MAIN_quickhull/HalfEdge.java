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
import java.util.HashMap;
import java.util.List;

import processing.core.PGraphics3D;


/**
 * 
 * 
 * 
 * TODO, when deleting mesh elements,the deleted elements also should be explicitly 
 * unlinked, in case old links remain (only in one direction) this may cause problems!!!
 */





/**
 * 
 * @author (c) Thomas Diewald, 2013
 * 
 * 
 *  Half Edge Data Structure
 *  
 *  ____________________________________________________________________________
 *  mesh structure:
 *  
 * 
 *  each "edge" consists of TWO half-edges (pairwise linked to each other),
 *  that point in the exact opposite direction, are of the same length, and share
 *  two adjacent faces.
 *  
 *  each half-edge is linked to ONE pair half-edge on an adjacent face (and vice versa).
 *  each half-edge is linked to ONE next half-edge on the same face (CCW).
 *  each half-edge is linked to ONE face which is on the left side it. (CCW)
 *  each half-edge is linked to ONE vertex that is its origin.
 *  each face      is linked to ONE of its half-edges.
 *  each vertex    is linked to ONE of its half-edges.
 *  
 *  each face is encircled by a single-linked-list of half-edges (CCW).
 *  each edge(= 2 half-edges) is shared by exactly two faces (2-manifold).
 *  so, the mesh can't have any "open" edges (=half-edges that have no pair-link).
 *  this means, a mesh has no open boundary and therefore divides the space R3
 *  into inside and outside.
 *  
 *  valence (number of adjacent ...):
 *  vertex: valence >= 2
 *  edge:   valence == 2
 *  face:   valence == 3
 *  
 *  (average vertex valence is 6, see below)
 *  
 *  ____________________________________________________________________________
 *  mesh creation:
 *  
 *  1 vertex: invalid
 *    edges: 0 
 *    faces: 0 
 *    volume: false
 *    
 *  2 vertices: ~ valid
 *    edges: 1 (= 2 half-edges)
 *    faces: 0
 *    volume: false
 *    
 *  3 vertices: ~ valid
 *    edges: 3 (= 6 half-edges)
 *    faces: 2 (back and front, on the same plane)
 *    volume: false
 *    face-normals: true
 *    vertex-normals: false
 *    shape: triangle
 *    
 *  4 vertices: VALID
 *    edges: 6 (= 12 half-edges)
 *    faces: 4
 *    volume: true
 *    face-normals: true
 *    vertex-normals: true
 *    shape: tetrahedron
 *  
 *  
 *  ____________________________________________________________________________
 *  additional:
 *  
 *  euler formula for closed, 2-manifold meshes: 
 *  
 *  V - E + F = 2*(1-g)
 *
 *  V ... number of vertices
 *  E ... number of edges
 *  F ... number of faces
 *  g ... genus (=number of "holes"/"handles", sphere.genus: 0, torus.genus: 1)
 *  
 * example tetrahedron: V(4) - E(6) + F(4) = 2*( 1-g(0) );
 * 
 * other formulas:
 * 2*E = 3*F  --> E = 1.5*F
 * 2*E = HE ( 1 face per half-edge )
 * 3*F = HE
 * 
 * F = 2*V-4; // when genus = 1
 * 
 * often used estimation: V - E + F ~ 0
 * 
 * --> V - E + F ~ 0 // replacing: F = (2/3)*E
 * --> V - E + (2/3)*E ~ 0
 * --> V - (1/3)*E ~ 0
 * --> V  ~ E/3
 * --> E  ~ 3*V
 * --> HE ~ 6*V (= average valence of 6)
 * 
 * 
 * 
 * links:
 * http://wiki.delphigl.com/index.php/Mesh
 * TODO 
 * 
 * 
 */
public abstract class HalfEdge{

  /**
   *  TODO TODO
   * mesh analysis, 
   * error checking, 
   * link checking
   * hole filling
   * face(s) extraction 
   * 
   * LoD (vertex/edge/face - collapse/split)
   * 
   * 
   */

  
  
  
  public static final double EPSILON = 0.001;
  
  //============================================================================
  //
  // DATA STRUCTURE
  //
  //============================================================================
  
  public static class Mesh {
    public Edge edge;
    public Flag flag;
    public MeshData data;
    public Mesh(){
      flag = new Flag();
      data = new MeshData();
    }
    //TODO: face-count, edge-count, vert-count - keep counters?
  }
  
  public static class Edge {
    public Edge next;
    public Edge pair;
    public Vert orig;
    public Face face;
    public Flag flag;
    public EdgeData data;
    
    public Edge(){
      flag = new Flag();
      data = new EdgeData();
    }
  }
  
  public static class Vert {
    public Edge edge;
    public Flag flag;
    public VertData data;
    
    public Vert(float[] v){
      flag = new Flag();
      data = new VertData();
      data.v = v;
    }
  }
  
  public static class Face{
    public Edge edge;
    public Flag flag;
    public FaceData data;
    
    public Face(){
      flag = new Flag();
      data = new FaceData();
    }
  }
  
  public static class Flag {
    public int selection = 0;
    public int visited   = 0;
  }

  
  public static class VertData{
    public float[] v; // vertex: Float3
    public float[] t; // texture coord: Float2 or null
    public float[] n; // normal: Float3 or null TODO: for hard edges/vertices, the face normal should be used, so some flag might be needed!
    public int id; // TODO: use for conversion to/from ifs !?!
    //TODO: shared, reference(ID,Object)???
    public int tmp_val;  // used for: slicing
    public float d; // use for: slicing[distance to plane.
    public int tmp_val_2 = 0;
  }
  public static class FaceData{
    public float n[], d; // plane: n=normal, d=distance to origin
    public float c[];    // center
    public int tmp_val;  // tmp_val can be used for anything
    public int tmp_val_2 = -1;  // tmp_val can be used for anything
//    int tmp_val_3 = 0;  // tmp_val can be used for anything
//    int tmp_val_4 = 0;  // tmp_val can be used for anything
//    TODO: number of vertices?
//    ArrayList<float[]> points;

    public int id = -1;
    
    
    public float distance(float[] f3_point){
      return (DwVec3.dot(f3_point, n) - d);
    }
  }
  public static class EdgeData{
//    float[] dir = new float[3];
  }
  public static class MeshData{
    String name = "halfedge_mesh_unnamed";
  }
  
  
  
  
  
  
  //============================================================================
  //
  // UPDATE/COMPUTE
  //
  //============================================================================
  
  public static abstract class Update{
    
    public static void faces(HalfEdge.Mesh mesh){
      HalfEdge.Face[] list = Collect.faces(null, mesh);
      for(HalfEdge.Face item : list) Update.face(item);
    }

    /**
     * update/compute face data of a given face.<br>
     * Center (c)<br>
     * Plane: normal (n) and distance (d) to origin<br>
     * 
     * @param face
     */
    public static void face(HalfEdge.Face face){
      List<HalfEdge.Vert> list = new ArrayList<HalfEdge.Vert>();
      Query.adjacentVerts(list, face);

      float[] C = new float[3];
      float[] N = new float[3];
      float D;
      
      // center
      for( HalfEdge.Vert item : list ){
        DwVec3.add_ref_slf(item.data.v, C);
      }
      DwVec3.scale_ref_slf(C, 1f/list.size());
    
      // normal, TODO: for non-planar faces (or faces with more than 3 vertices), 
      // an average normal should be generated
      final float[] E1 = DwVec3.sub_new(list.get(1).data.v, C);
      final float[] E2 = DwVec3.sub_new(list.get(0).data.v, C);
      DwVec3.cross_ref(E1, E2, N); // TODO
      DwVec3.normalize_ref_slf(N);
      
      // distance to origin
      D = DwVec3.dot(C, N);
      
      face.data.c = C;
      face.data.n = N;
      face.data.d = D;
    }
    
    private static final float _1_DIV_3_ = 1f/3f;
    
    
    /**
     * returns true, if the triangle is valid
     * @param face
     * @return
     */
    public static boolean faceTriangle(HalfEdge.Face face){
      // vertices
      float[] v0 = face.edge          .orig.data.v;
      float[] v1 = face.edge.next     .orig.data.v;
      float[] v2 = face.edge.next.next.orig.data.v;
      // edges, TODO orientation, left-handed/right-handed?
      float[] E1 = DwVec3.sub_new(v2, v0);   
      float[] E2 = DwVec3.sub_new(v1, v0);  
      
      ///// TRIANGLE - CENTER
      float[] C = new float[3];
      DwVec3.add_ref_slf(v0, C);
      DwVec3.add_ref_slf(v1, C);
      DwVec3.add_ref_slf(v2, C);
      DwVec3.scale_ref_slf(C, _1_DIV_3_);
      face.data.c = C;

      ///// TRIANGLE - PLANE
      // normal
      float[] N = DwVec3.cross_new(E1, E2); // TODO
      float mag_sq = DwVec3.mag_sq(N);
      if( mag_sq == 0.0 ) {
        System.err.println("faceTriangle(HalfEdge.Face face): triangle is degenerated");
        return false; // triangle is degenerate !!!! -> delete triangle
      }
      DwVec3.scale_ref_slf(N, (float)(1d/Math.sqrt(mag_sq)));

      face.data.n = N;
      face.data.d = DwVec3.dot(C, N); // distance to origin
      
      return true;
    }
    
    /**
     * computes/updates the vertex normals of the whole mesh.<br>
     * mostly, vertex-normals are only need for rendering, or debugging.<br>
     * 
     * @param mesh mesh, the vertices are update.
     */
    public static void verts(HalfEdge.Mesh mesh){
      HalfEdge.Vert[] list = Collect.verts((HalfEdge.Vert[])null, mesh);
//      System.out.println("number of vertices to update: "+list.size());
      for(HalfEdge.Vert item : list) Update.vert(item);
    }
    
    /**
     * computes/updates the given vertex normal of the mesh.<br>
     * mostly, vertex-normals are only need for rendering, or debugging.<br>
     * @param vert vertex to compute.
     */
    public static boolean vert(HalfEdge.Vert vert){

      // get adjacent vertices
      List<HalfEdge.Vert> list = new ArrayList<HalfEdge.Vert>();
      Query.adjacentVerts(list, vert);

      final int num_verts = list.size();
      
      // special case
      // TODO interpolate normal? from neighbors?
      if( num_verts == 2 ){
        
        HalfEdge.Face face_L = vert.edge.face;
        HalfEdge.Face face_R = vert.edge.pair.face;
        
        float[] N = DwVec3.add_new(face_L.data.n, face_R.data.n );
        DwVec3.scale_ref_slf(N, 0.5f);
        vert.data.n = N;
        return true;
      }
      
      // compute edges to adjacent vertices
      final float[]   V = vert.data.v;
      final float[][] E = new float[num_verts][3]; // edges to adjacent vertices
      for( int i = 0; i < num_verts; i++ ){
        DwVec3.sub_ref(list.get(i).data.v, V, E[i]);
        DwVec3.normalize_ref_slf(E[i]);
      }
      if( num_verts == 2){
      System.out.println("_____");
      System.out.println(num_verts);
      }
      
      // compute vertex normal
      // sum up cross-products (normal vectors) of adjacent edges.
      // this takes the size of the face into account (if it is a triangle)
      // TODO: alternatives:
      //       sum up: cross-products of adjacent edges * angle-between-edges
      //       sup up: face-normals(normalized)*face_size
      //       sum up: face-normals(normalized)
      float[] N = new float[3];  // summed up cross-products (final normal)
      float[] CP = new float[3]; // cross-product (normal)
      float[] E1 = E[num_verts-1];
      for( int i = 0; i < num_verts; i++ ){
        
        DwVec3.cross_ref(E1, E[i], CP); // normal vector E1->E[i]

        float mag_sq = DwVec3.mag_sq(CP);
        if( mag_sq < EPSILON) {
          System.err.println("update vertex normal, crossproduct has length 0.0");
        }
        DwVec3.scale_ref_slf(CP, (float)(1d/Math.sqrt(mag_sq))); // normalize it
        
        float dot = DwVec3.dot(E1, E[i]); // cos-angle of E1->E[i], both normalized
        
        float s = (float)(Math.acos(dot)); //TODO?
        if( num_verts == 2){
        System.out.println(dot+", "+s);
        System.out.println("mag_sq = "+mag_sq);
        
        }
//      System.out.println("  dot = "+dot);
//      System.out.printf(Locale.ENGLISH, "    len E1 = %5.2f, E[i] =%5.2f\n",DwVec3.mag(E1), DwVec3.mag(E[i]));
//      System.out.printf(Locale.ENGLISH, "s =%5.2f,  n = [%5.2f, %5.2f, %5.2f], len = %5.2f\n", s, CP[0], CP[1], CP[2], DwVec3.mag(CP));
        DwVec3.scale_ref_slf(CP, s); // scale normal with by angle of E1->E[i]
        DwVec3.add_ref_slf(CP, N);
        E1=E[i];
      }
      DwVec3.normalize_ref_slf(N);
      
      vert.data.n = N;
      return true;
    }
  }
  
  
  
  //============================================================================
  //
  // QUERY
  //
  //============================================================================

  public static abstract class Query{
    /** get adjacent VERTICES of given VERTEX.<br> rotation: CW */
    public static void adjacentVerts(List<HalfEdge.Vert> list, HalfEdge.Vert vert){
      HalfEdge.Edge edge = vert.edge;
      do {   list.add(edge.pair.orig);    } while( (edge=edge.pair.next) != vert.edge );
    }
    /** get 2 adjacent VERTICES of given EDGE.<br> A -> B */
    public static void adjacentVerts(List<HalfEdge.Vert> list, HalfEdge.Edge edge){
      list.add( edge.orig);
      list.add( edge.next.orig);
    }
    /** get adjacent VERTICES of given FACE.<br> rotation: CCW */
    public static void adjacentVerts(List<HalfEdge.Vert> list, HalfEdge.Face face){
      HalfEdge.Edge edge = face.edge;
      do {   list.add(edge.orig);         } while( (edge=edge.next) != face.edge);
    }
    
    /** get adjacent (outgoing) EDGES of given VERTEX.<br> rotation: CW */
    public static void adjacentEdges(List<HalfEdge.Edge> list, HalfEdge.Vert vert){
      HalfEdge.Edge edge = vert.edge;
      do {   list.add(edge);              } while( (edge=edge.pair.next) != vert.edge);
    }
    /** get adjacent EDGES of given EDGE (edges of vertex A + edges of vertex B).<br> rotation: CW */
    public static void adjacentEdges(List<HalfEdge.Edge> list, HalfEdge.Edge edge){
      adjacentEdges(list, edge.orig);
      adjacentEdges(list, edge.next.orig);
    }
    /** get adjacent EDGES of given FACE.<br> rotation: CCW */
    public static void adjacentEdges(List<HalfEdge.Edge> list, HalfEdge.Face face){
      HalfEdge.Edge edge = face.edge;
      do {   list.add(edge);              } while( (edge=edge.next) != face.edge);
    }
    
    /** get adjacent FACES of given VERTEX.<br> rotation: CW */
    public static void adjacentFaces(List<HalfEdge.Face> list, HalfEdge.Vert vert){
      HalfEdge.Edge edge = vert.edge;
      do {   list.add(edge.face);         } while( (edge=edge.pair.next) != vert.edge);
    }
    /** get 2 adjacent (left+right) FACES of given EDGE. */
    public static void adjacentFaces(List<HalfEdge.Face> list, HalfEdge.Edge edge){
      list.add( edge.face);
      list.add( edge.pair.face);
    }
    /** get adjacent FACES (that share an edge) of given FACE.<br> rotation: CCW */
    public static void adjacentFaces(List<HalfEdge.Face> list, HalfEdge.Face face){
      HalfEdge.Edge edge = face.edge;
      do {   list.add(edge.pair.face);    } while( (edge=edge.next) != face.edge);
    }
    public static void adjacentFacesAll(List<HalfEdge.Face> list, HalfEdge.Face face){
      //TODO
      System.err.println("TODO: Query.adjacentFacesAll is not impleemented");
    }

    /**
     * test if a given face is planar.
     * @param face
     * @return true, if planar
     */
    public static boolean isPlanar(HalfEdge.Face face){
      List<HalfEdge.Vert> adj = new ArrayList<HalfEdge.Vert>();
      Query.adjacentVerts(adj, face);
      for(HalfEdge.Vert vert : adj){
        float dis = (DwVec3.dot(vert.data.v, face.data.n) - face.data.d); //TODO: Plane in face-data
        if( Math.abs(dis) > EPSILON ){ //TODO: relative epsilon?
          return false;
        }
      }
      return true;
    }
    /**
     * test all mesh-faces for beeing planar, and return them in a list.
     * 
     * @param mesh half-edge mesh
     * @return list of non-planar faces.
     */
    public static List<HalfEdge.Face> nonPlanarFaces(HalfEdge.Mesh mesh){
      // get all faces
      HalfEdge.Face[] faces = Collect.faces(null, mesh);
      // save non-planar faces
      List<HalfEdge.Face> nonplanar = new ArrayList<HalfEdge.Face>();
      for( HalfEdge.Face face : faces ){
        if( !Query.isPlanar(face) ) nonplanar.add(face);
      }
      return nonplanar;
    }
    
    
    //TODO: non convex faces
    
    
    
    /**
     * tries to find an edge, that shares both, face_A and face_B.
     * this method can also be used, for checking if two faces are adjacent.
     * 
     * @param face_A
     * @param face_B
     * @return the sharing edge (on face_A), or null, of the faces are not adjacent.
     */
    public static HalfEdge.Edge sharingEdge(HalfEdge.Face face_A, HalfEdge.Face face_B){
      HalfEdge.Edge start = face_A.edge;
      HalfEdge.Edge iter  = start;
      do {
        if( iter.pair.face == face_B) 
          return iter;
      } while( (iter=iter.next) != start );
      
      return null;
    }
  }
  
  
  
  
  
  
  
  
  
  //============================================================================
  //
  // COLLECT NODES
  //
  //============================================================================
  
  public static abstract class Collect{

    /**
     * collect all faces of the mesh.<br>
     * stack, is used for iteration, and also for storing the faces.<br>
     * the method will have the best performance if stack is allocated to about <br>
     * the number of faces that are expected.<br>
     * anyways, the returned list contains all faces.
     */
    public static HalfEdge.Face[] faces(HalfEdge.Face[] stack, HalfEdge.Mesh mesh){

      final int flag = ++mesh.flag.visited;
  
      HalfEdge.Face face = mesh.edge.face;
 
      int ptr_beg = 0;
      int ptr_end = 0;
      if( stack == null ) stack = new HalfEdge.Face[100];
      stack[ptr_end++] = mesh.edge.face;

      face.flag.visited = flag;
      while( ptr_beg < ptr_end ){

        face = stack[ptr_beg++];
        HalfEdge.Edge start = face.edge;
        HalfEdge.Edge edge  = start;
        do {     
          face = edge.pair.face;
          if( face.flag.visited != flag ){
            if( ptr_end >= stack.length){
              stack = Arrays.copyOf(stack, ptr_end<<1);
            }
            stack[ptr_end++] = face;
            face.flag.visited = flag;
          }
        } 
        while( (edge=edge.next) != start);
      }
      
      if( stack.length > ptr_end )
        stack = Arrays.copyOf(stack, ptr_end);
      return stack;
    }
    

 
    /**
     * collect all edges of the mesh.<br>
     * stack, is used for iteration, and also for storing the edges.<br>
     * the method will have the best performance if stack is allocated to about <br>
     * the number of edges that are expected.<br>
     * anyways, the returned list contains all edges, this is: only
     * one half-edge per edge, although both are flagged as visited.
     */
    public static HalfEdge.Edge[] edges(HalfEdge.Edge[] stack, HalfEdge.Mesh mesh){
     
      final int flag = ++mesh.flag.visited;
      HalfEdge.Edge edge = mesh.edge;
      
      int ptr_beg = 0;
      int ptr_end = 0;
      if( stack == null ) stack = new HalfEdge.Edge[100];
      stack[ptr_end++] = mesh.edge;
      
      edge.flag.visited = flag;
      while( ptr_beg < ptr_end ){
        
        edge = stack[ptr_beg++];
        
        HalfEdge.Edge start    = edge.pair;
        HalfEdge.Edge edge_adj = start;
        do {   
          if( edge_adj.     flag.visited != flag &&  // only one of the two half-edges
              edge_adj.pair.flag.visited != flag )  
          { 
            if( ptr_end >= stack.length){
              stack = Arrays.copyOf(stack, ptr_end<<1);
            }
            stack[ptr_end++] = edge_adj;
            edge_adj.flag.visited = flag;
          }                   
        } while( (edge_adj=edge_adj.pair.next) != start);

      }
      if( stack.length > ptr_end )
        stack = Arrays.copyOf(stack, ptr_end);
      return stack;
    }
    

    /**
     * collect all verts of the mesh.<br>
     * stack, is used for iteration, and also for storing the verts.<br>
     * the method will have the best performance if stack is allocated to about <br>
     * the number of verts that are expected.<br>
     * anyways, the returned list contains all verts.<br>
     * during collecting, the vertex' ID is set, based on the position in the list.<br>
     * this is useful for converting to an IFS-representation.
     */
    public static HalfEdge.Vert[] verts(HalfEdge.Vert[] stack, HalfEdge.Mesh mesh){

      int flag = ++mesh.flag.visited;
      HalfEdge.Vert vert = mesh.edge.orig;
   
      int ptr_beg = 0;
      int ptr_end = 0;
      if( stack == null ) stack = new HalfEdge.Vert[100];
      vert.flag.visited = flag;
      vert.data.id = ptr_end;
      stack[ptr_end++] = vert;

      while( ptr_beg < ptr_end ){
        vert = stack[ptr_beg++];

        HalfEdge.Edge edge = vert.edge.pair;
        HalfEdge.Edge iter = edge;
        do{
          if( iter.orig.flag.visited != flag ){
            if( ptr_end >= stack.length){
              stack = Arrays.copyOf(stack, ptr_end<<1);
            }
            iter.orig.flag.visited = flag;
            iter.orig.data.id = ptr_end;
            stack[ptr_end++] = vert = iter.orig;
          }
        } while( (iter=iter.next.pair) != edge);
      }

      if( stack.length > ptr_end )
        stack = Arrays.copyOf(stack, ptr_end);
      return stack;
    }
  }
  
  
  

  
  
  
  
  
  
  //============================================================================
  //
  // SELECT
  //
  //============================================================================

  
  
  public static abstract class Selection {
    public final static int DESELECTION_MARK = 0;
    protected int mark = 1;
    protected int size = 0;
    
    public int size(){ 
      return size;
    }
    public void setMark(int mark){
      this.mark = (mark <= DESELECTION_MARK) ? DESELECTION_MARK+1 : mark;
    }
    public int getMark(){
      return mark;
    }
    public void increaseMark(){
      mark++;
    }
    public void decreaseMark(){
      setMark(--mark);
      
    }
    static public boolean isSelected(HalfEdge.Face face){ return (face.flag.selection != DESELECTION_MARK); }
    static public boolean isSelected(HalfEdge.Edge edge){ return (edge.flag.selection != DESELECTION_MARK); }
    static public boolean isSelected(HalfEdge.Vert vert){ return (vert.flag.selection != DESELECTION_MARK); }
    
    abstract public void select(HalfEdge.Mesh mesh);
    abstract public void clear();
    abstract public void grow();
    abstract public void shrink();
    abstract public void display(PGraphics3D g);
  
    
    
    public static class Face extends Selection{
      public ArrayList<HalfEdge.Face> selection = new ArrayList<HalfEdge.Face>();

      public void select(HalfEdge.Face face){
        if( !isSelected(face)){
          selection.add(face);
          face.flag.selection = mark;
          size++;
        }
      }
      
      @Override
      public void clear(){
        for( HalfEdge.Face face : selection) 
          face.flag.selection = DESELECTION_MARK;
        selection.clear();
        size = 0;
        setMark(1);
      }
     
      @Override
      public void grow(){
        increaseMark();
        ArrayList<HalfEdge.Face> adj = new ArrayList<HalfEdge.Face>();
        for( HalfEdge.Face face : selection){
          Query.adjacentFaces(adj, face);
        }
        for( HalfEdge.Face item : adj ) select(item);
      }

      @Override
      public void display(PGraphics3D g) {
        HalfEdge_Render.selection(g, this);
      }

      @Override
      public void shrink() {
//        System.err.println("Selection.Face.shrink() not implemented yet");
        // better to swapping, so that faces that are going to be unselected.
        // are moved to the end. then get a view of the unselected portion and delete.
        ArrayList<HalfEdge.Face> A = new ArrayList<HalfEdge.Face>(); // remains selected
        ArrayList<HalfEdge.Face> B = new ArrayList<HalfEdge.Face>(); // gets unselected
        
        for(int i = 0; i < selection.size(); i++){
          HalfEdge.Face face = selection.get(i);
          ArrayList<HalfEdge.Face> adj = new ArrayList<HalfEdge.Face>();
          Query.adjacentFaces(adj, face);
          boolean deselect = false;
          for( HalfEdge.Face item : adj ) {
            if( !isSelected(item) ){
              deselect = true;
              break;
            }
          }
          if( deselect )
            B.add(face);
           else 
            A.add(face);
        }
        
        if( B.size() > 0 ){
          decreaseMark();
          for( HalfEdge.Face face : B) 
            face.flag.selection = DESELECTION_MARK;
        }
        selection = A;
        size = selection.size();
      }

      @Override
      public void select(Mesh mesh) {
        select(mesh.edge.face);
      }
    }
    
    
    public static class Edge extends Selection {
      public ArrayList<HalfEdge.Edge> selection = new ArrayList<HalfEdge.Edge>();


      static public boolean isSelected(HalfEdge.Edge edge){
        return (edge.flag.selection != DESELECTION_MARK);
      }
      
      public void select(HalfEdge.Edge edge){
        // only one of the two half-edges needs to be selected
        if( !isSelected(edge) && !isSelected(edge.pair)){ 
          selection.add(edge);
          edge.flag.selection = mark;
          size++;
        }
      }
      
      @Override
      public void clear(){
        for( HalfEdge.Edge edge : selection)
          edge.flag.selection = DESELECTION_MARK;
        selection.clear();
        size = 0;
        setMark(1);
      }
     
      @Override
      public void grow(){
        increaseMark();
        ArrayList<HalfEdge.Edge> adj = new ArrayList<HalfEdge.Edge>();
        for(int i = 0; i < selection.size(); i++){
          Query.adjacentEdges(adj, selection.get(i));
        }
        for( HalfEdge.Edge item : adj ) select(item);
      }

      @Override
      public void shrink() {
        System.err.println("Selection.Edge.shrink() not implemented yet");
      }

      @Override
      public void display(PGraphics3D g) {
        HalfEdge_Render.selection(g, this);
      }
      @Override
      public void select(Mesh mesh) {
        select(mesh.edge);
      }
    }
    
    
    
    public static class Vert extends Selection {
      public ArrayList<HalfEdge.Vert> selection = new ArrayList<HalfEdge.Vert>();

      public void select(HalfEdge.Vert vert){
        if( !isSelected(vert)){ 
          selection.add(vert);
          vert.flag.selection = mark;
          size++;
        }
      }
      
      @Override
      public void clear(){
        for( HalfEdge.Vert vertex : selection)
          vertex.flag.selection = DESELECTION_MARK;
        selection.clear();
        size = 0;
        setMark(1);
      }
     
      @Override
      public void grow(){
        int old_size = size;
        ArrayList<HalfEdge.Vert> adj = new ArrayList<HalfEdge.Vert>();
        for(int i = 0; i < selection.size(); i++){
          Query.adjacentVerts(adj, selection.get(i));
        }
        for( HalfEdge.Vert item : adj ) select(item);
        if( size-old_size > 0) increaseMark();
      }
      @Override
      public void shrink() {
        // TODO: optimized

        // better to swapping, so that vertices that are going to be unselected.
        // are moved to the end. then get a view of the unselected portion and delete.
        ArrayList<HalfEdge.Vert> A = new ArrayList<HalfEdge.Vert>(); // remains selected
        ArrayList<HalfEdge.Vert> B = new ArrayList<HalfEdge.Vert>(); // gets unselected
        
        for(int i = 0; i < selection.size(); i++){
          HalfEdge.Vert vert = selection.get(i);
          ArrayList<HalfEdge.Vert> adj = new ArrayList<HalfEdge.Vert>();
          Query.adjacentVerts(adj, selection.get(i));
          boolean deselect = false;
          for( HalfEdge.Vert item : adj ) {
            if( !isSelected(item) ){
              deselect = true;
              break;
            }
          }
          if( deselect )
            B.add(vert);
           else 
            A.add(vert);
        }
        
        if( B.size() > 0 ){
          decreaseMark();
          for( HalfEdge.Vert vert : B) 
            vert.flag.selection = DESELECTION_MARK;
        }
        selection = A;
        size = selection.size();
      }

      @Override
      public void display(PGraphics3D g) {
        HalfEdge_Render.selection(g, this);
      }
      @Override
      public void select(Mesh mesh) {
        select(mesh.edge.orig);
      }
    }

  }
  
  
  
  
  
  
  
  
  
  
  
  
  public static abstract class MeshOps{
    
    /**
     * test all mesh faces for planarity. those being not planar are triangulated.<br>
     * using the face-center as new point.<br>
     * 
     * @param mesh
     * @return number of new new faces
     */
    public static int trianguleNonPlanarFaces(HalfEdge.Mesh mesh){
      int num_new_faces = 0;
      List<HalfEdge.Face> faces_non_planar = HalfEdge.Query.nonPlanarFaces(mesh);
      for( HalfEdge.Face item : faces_non_planar ){
        num_new_faces += HalfEdge.MeshOps.triangulateFace(mesh, item, (HalfEdge.Face[])null);
      }
      HalfEdge.Update.verts(mesh);
      return num_new_faces;
    }

        
    /**
     * triangulate given face using its face-center as a new mesh-vertex. <br>
     * @param mesh
     * @param face
     * @param faces_new
     * @return number of new faces
     */
    public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Face face, List<HalfEdge.Face> faces_new){
      return triangulateFace(mesh, face, face.data.c, faces_new);
    }
    public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Face face, HalfEdge.Face[] faces_new){
      return triangulateFace(mesh, face, face.data.c, faces_new);
    }
    /**
     * triangulate given face using "point" as a new mesh-vertex. <br>
     * 
     * "face" fill be replaced by new faces.<br>
     * <br>
     * mesh.edge.orig = point<br>
     * <br>
     * @param mesh
     * @param face face that gets triangulated, and replaced by new faces.
     * @param point new mesh-point
     * @param faces_new list, to save the new faces to
     * @return number of new faces
     */
    public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Face face, float[] point, List<HalfEdge.Face> faces_new){
      
      HalfEdge.Vert vert__ = new HalfEdge.Vert(point);
      
      HalfEdge.Edge edge_start = face.edge;
      HalfEdge.Edge edge_cur   = edge_start;
      
      edge_cur.orig.edge = edge_cur;

      do{
        HalfEdge.Face face__ = new HalfEdge.Face();
        HalfEdge.Edge edge_2 = new HalfEdge.Edge();
        HalfEdge.Edge edge_1 = new HalfEdge.Edge();
        HalfEdge.Edge edge_0 = edge_cur;

        edge_cur = edge_cur.next;

        // linking: FACE <-> EDGES
        face__.edge = edge_0;
        edge_0.face = face__; 
        edge_1.face = face__;
        edge_2.face = face__;
        
        // linking: EDGE ORIGINS
        // !!!! in case the same vertex appears multiple times in an edge-ring edge-ring
        // which must NOT happen, but can (during convex-hull creation due to precision errors)
        // some created faces will be skipped in the next loop. therefore the saving and of
        // the faces is done in the next loop.
        // edit: during convex hull the horizon is granted to be "nice".
        edge_0.orig.edge  = edge_0;  // -->
        edge_1.orig = edge_cur.orig;
        edge_2.orig = vert__;
             
        // linking: NEXT EDGES
        edge_0.next = edge_1;
        edge_1.next = edge_2;
        edge_2.next = edge_0;
        
        // update new created face
        if( !HalfEdge.Update.faceTriangle(face__) ){
          //TODO
          System.err.println( "TRIANGLE IS DEGENERATE, cant comput normal");
        }
        
        
      } while( edge_cur != edge_start);
      
      int num_new_faces = 0;
      do{

        if(faces_new!=null) faces_new.add(edge_cur.face);
        num_new_faces++;

        HalfEdge.Edge edge_1_cur = edge_cur.next;
        HalfEdge.Edge edge_2_nxt = edge_cur.next.orig.edge.next.next;
        
        edge_1_cur.pair = edge_2_nxt;
        edge_2_nxt.pair = edge_1_cur;
        
//        edge_cur = edge_2_nxt.next;
        edge_cur = edge_cur.next.orig.edge;

      } while( edge_cur != edge_start);

      vert__.edge = edge_start.next.pair;
      mesh.edge = vert__.edge;
      
      //TODO: update everything that was in the neighborhood (vertices, faces)
      return num_new_faces;
    }

   
    /**
     * triangulate given face using "point" as a new mesh-vertex. <br>
     * 
     * "face" fill be replaced by new faces.<br>
     * <br>
     * mesh.edge.orig = point<br>
     * <br>
     * @param mesh
     * @param face face that gets triangulated, and replaced by new faces.
     * @param point new mesh-point
     * @param faces_new list, to save the new faces to
     * @return number of new faces
     */
   public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Face face, float[] point, HalfEdge.Face[] faces_new){
      
      HalfEdge.Vert vert__ = new HalfEdge.Vert(point);
      
      HalfEdge.Edge edge_start = face.edge;
      HalfEdge.Edge edge_cur   = edge_start;
      
      edge_cur.orig.edge = edge_cur;

      do{
        HalfEdge.Face face__ = new HalfEdge.Face();
        HalfEdge.Edge edge_2 = new HalfEdge.Edge();
        HalfEdge.Edge edge_1 = new HalfEdge.Edge();
        HalfEdge.Edge edge_0 = edge_cur;

        edge_cur = edge_cur.next;

        // linking: FACE <-> EDGES
        face__.edge = edge_0;
        edge_0.face = face__; 
        edge_1.face = face__;
        edge_2.face = face__;
        
        // linking: EDGE ORIGINS
        // !!!! in case the same vertex appears multiple times in an edge-ring edge-ring
        // which must NOT happen, but can (during convex-hull creation due to precision errors)
        // some created faces will be skipped in the next loop. therefore the saving and of
        // the faces is done in the next loop.
        // edit: during convex hull the horizon is granted to be "nice".
        edge_0.orig.edge  = edge_0;  // -->
        edge_1.orig = edge_cur.orig;
        edge_2.orig = vert__;
             
        // linking: NEXT EDGES
        edge_0.next = edge_1;
        edge_1.next = edge_2;
        edge_2.next = edge_0;
        
        // update new created face
        if( !HalfEdge.Update.faceTriangle(face__) ){
          //TODO
          System.err.println( "TRIANGLE IS DEGENERATE, cant comput normal");
        }
        
        
      } while( edge_cur != edge_start);
      
      int num_new_faces = 0;
      do{

        if(faces_new!=null) faces_new[num_new_faces] = edge_cur.face;
        num_new_faces++;

        HalfEdge.Edge edge_1_cur = edge_cur.next;
        HalfEdge.Edge edge_2_nxt = edge_cur.next.orig.edge.next.next;
        
        edge_1_cur.pair = edge_2_nxt;
        edge_2_nxt.pair = edge_1_cur;
        
//        edge_cur = edge_2_nxt.next;
        edge_cur = edge_cur.next.orig.edge;

      } while( edge_cur != edge_start);

      vert__.edge = edge_start.next.pair;
      mesh.edge = vert__.edge;
      
      //TODO: update everything that was in the neighborhood (vertices, faces)
      return num_new_faces;
    }
    
    
//    //TODO: delete
//    public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Face face, float[] point, FaceList faces_new){
//      
//      int num_new_faces = 0;
//      
//      HalfEdge.Vert vert__ = new HalfEdge.Vert(point);
//      
//      HalfEdge.Edge edge_start = face.edge;
//      HalfEdge.Edge edge_cur   = edge_start;
//
//      do{
//        HalfEdge.Face face__ = new HalfEdge.Face();
//        HalfEdge.Edge edge_2 = new HalfEdge.Edge();
//        HalfEdge.Edge edge_1 = new HalfEdge.Edge();
//        HalfEdge.Edge edge_0 = edge_cur;
//
//        edge_cur = edge_cur.next;
//
//        // linking: FACE <-> EDGES
//        face__.edge = edge_0;
//        edge_0.face = face__; 
//        edge_1.face = face__;
//        edge_2.face = face__;
//        
//        // linking: EDGE ORIGINS
//        edge_0.orig.edge  = edge_0; // !!!!
//        edge_1.orig = edge_cur.orig;
//        edge_2.orig = vert__;
//             
//        // linking: NEXT EDGES
//        edge_0.next = edge_1;
//        edge_1.next = edge_2;
//        edge_2.next = edge_0;
//        
//        // update new created face
//        if( !HalfEdge.Update.faceTriangle(face__) ){
//          //TODO
//          System.err.println( "TRIANGLE IS DEGENERATE, cant comput normal");
//        }
//
//        // add to face list
//        if( faces_new != null ) faces_new.add(face__);
//        
//        num_new_faces++;
//        
//      } while( edge_cur != edge_start);
//
//      do{
//        HalfEdge.Edge edge_1_cur = edge_cur.next;
//        HalfEdge.Edge edge_2_nxt = edge_cur.next.orig.edge.next.next;
//
//        edge_1_cur.pair = edge_2_nxt;
//        edge_2_nxt.pair = edge_1_cur;
//        edge_cur = edge_2_nxt.next;
//      } while( edge_cur != edge_start);
//      
//      vert__.edge = edge_start.next.pair;
//      mesh.edge = vert__.edge;
//      
//      //TODO: update everything that was in the neighbourhood (vertices, faces)
//      
//      return num_new_faces;
//    }
    
    
    /**
     * deletes a given vertex from a mesh, including its adjacent edges and faces.<br>
     * the resulting face might not be planar!<br>
     * <br>
     * TODO: check old links, check id ok, when deleting multiple edges (creating a ring of faces)<br>
     * 
     * @param mesh  mesh, to delete the vertex from
     * @param vert  vertex to delete
     * @return number of deleted edges ( = number of deleted faces +1)
     */
    public static int deleteVert(HalfEdge.Mesh mesh, HalfEdge.Vert vert){
      int num_removed_edges = 0;
      
      HalfEdge.Edge start = vert.edge.next;
      start.face.edge = start;
      
      HalfEdge.Edge edge = start;

      do{
        edge.face = start.face;
        if( edge.next.pair.orig == vert){
          edge.next = edge.next.pair.next;
          edge.next.orig.edge = edge.next;
          num_removed_edges++;
        }
      }
      while( (edge = edge.next) != start );
      
      HalfEdge.Update.face(start.face); // face is not a triangle
      //TODO: update everything that was in the neighbourhood (vertices, faces)

      mesh.edge = start;
      return num_removed_edges;
    }
    
    
    
    /**
     * deletes a given edge from the mesh -> adjacent faces are merged to 1 face.<br>
     * the resulting face might not be planar.<br>
     * <br>
     * mesh.edge = edge.next<br>
     * 
     * @param mesh
     * @param edge
     * @return number of removed faces, which is always 1
     */
    public static int deleteEdge(HalfEdge.Mesh mesh, HalfEdge.Edge edge){

      HalfEdge.Face face = edge.face; // this face remains part of the mesh
      face.edge = edge.next;
      
      HalfEdge.Edge edge_A = edge;
      HalfEdge.Edge edge_B = edge.pair;
      
      HalfEdge.Vert vert_A = edge_A.orig;
      HalfEdge.Vert vert_B = edge_B.orig;
      
      vert_A.edge = edge_B.next;
      vert_B.edge = edge_A.next;
      
      HalfEdge.Edge prev;
      
      prev = edge_A;
      while((prev=prev.next).pair.orig != vert_A);
      prev.next = edge_B.next;
       
      prev = edge_B;
      while((prev=prev.next).pair.orig != vert_B) prev.face = face;
      prev.face = face;
      prev.next = edge_A.next;
      
      
      // control loop, if everything is linked correctly
//      {
//        HalfEdge.Edge cur, start;
//        cur = start = edge.next;
//        int it = 0;
//  
//        do {
//          if( cur.face != face ){
//            System.out.println("(HalfEdge.MeshOps.delete(Mesh mesh, Edge edge) EDGE NOT LINKED TO RIGHT FACE "+it);
//          }
//          it++;
//        } while( (cur=cur.next) != start);
//        System.out.println("number of edges = "+it);
//      }
      

      HalfEdge.Update.face(face);
      //TODO: update everything that was in the neighborhood (vertices, faces?)

      mesh.edge = edge.next;
      return 1; // return number of removed faces, ...always 1 here
    }
    
    
    public static int delete(HalfEdge.Mesh mesh, HalfEdge.Face edge){
      //TODO
      return 0;
    }
    
    
    /**
     * flips and edge, that is shared by two TRIANGLES!
     * 
     * @param mesh
     * @param edge
     */
    
    
    
    /*
 
     
     CCW
     
     |           |          |           |  
   --o-----------o--      --o-----------o-- 
     |         / |          | \         |  
     |  A    /   |          |   \   B   |  
     |     /     |    ->    |     \     |  
     |   /   B   |          |  A    \   |  
     | /         |          |         \ |  
   --o-----------o--      --o-----------o-- 
     |           |          |           | 
     */
    /**
     * 
     * flips an edge with two adjacent triangles.
     * 
     * 
     * @param mesh
     * @param edge
     */
    public static void edgeFlip(HalfEdge.Mesh mesh, HalfEdge.Edge edge){
  
      // copy old values
      final HalfEdge.Edge edge_A_0 = edge;
      final HalfEdge.Edge edge_A_1 = edge_A_0.next;
      final HalfEdge.Edge edge_A_2 = edge_A_1.next;
      
      final HalfEdge.Edge edge_B_0 = edge.pair;
      final HalfEdge.Edge edge_B_1 = edge_B_0.next;
      final HalfEdge.Edge edge_B_2 = edge_B_1.next;
      
      // check, if faces are triangles
      if( edge_A_2.next != edge_A_0) { System.out.println("(edgeFlip) 1. face is not a triangle."); }
      if( edge_B_2.next != edge_B_0) { System.out.println("(edgeFlip) 2. face is not a triangle."); }
            
      // EDGES
      // link next-edges for the flip-edge
      edge_A_0.next = edge_A_2;
      edge_B_0.next = edge_B_2;
      
      // link next.next-edges for the flip-edge
      edge_A_0.next.next = edge_B_1;
      edge_B_0.next.next = edge_A_1;
      
      // link next.next.next-edges for the flip-edge ( to close the loop )
      edge_A_0.next.next.next = edge_A_0;
      edge_B_0.next.next.next = edge_B_0;
      
      // FACES
      // assure that adjacent faces are linked to the flip-edge.
      edge_A_0.face.edge = edge_A_0;
      edge_B_0.face.edge = edge_B_0;
      
      // link prev-edge to current face
      edge_A_0.next.next.face = edge_A_0.face;
      edge_B_0.next.next.face = edge_B_0.face;
      
      // VERTS
      // link new origins of flip-edge
      edge_A_0.orig = edge_B_2.orig;
      edge_B_0.orig = edge_A_2.orig;
      
      // assure that flip-edge origins point to flip-edge
      edge_A_0.orig.edge = edge_A_0;
      edge_B_0.orig.edge = edge_B_0;
      
      // assure that opposite origins link to next-next
      edge_A_0.next.next.orig.edge = edge_A_0.next.next;
      edge_B_0.next.next.orig.edge = edge_B_0.next.next;
      
      
      // UPDATE FACES
//      HalfEdge.Update.faceTriangle(edge_A_0.face);
//      HalfEdge.Update.faceTriangle(edge_B_0.face);
      if( !HalfEdge.Update.faceTriangle(edge_A_0.face) ){
        //TODO
        System.err.println( "TRIANGLE IS DEGENERATE, cant comput normal");
      }
      if( !HalfEdge.Update.faceTriangle(edge_B_0.face) ){
        //TODO
        System.err.println( "TRIANGLE IS DEGENERATE, cant comput normal");
      }
      
      
      mesh.edge = edge;
      
      
      // error checking - TODO only for debugging purposes
      int count;
      HalfEdge.Edge start;
      HalfEdge.Edge iter;
     
      if( edge_A_0.orig.edge != edge_A_0 ) System.err.println("(edgeFlip) error: vert-edge not linked 1");
      if( edge_B_0.orig.edge != edge_B_0 ) System.err.println("(edgeFlip) error: vert-edge not linked 2 ");
      
      if( edge_A_0.next.next.orig.edge != edge_A_0.next.next ) System.err.println("(edgeFlip) error: vert-edge not linked 3");
      if( edge_B_0.next.next.orig.edge != edge_B_0.next.next ) System.err.println("(edgeFlip) error: vert-edge not linked 4");
      
      start = iter = edge_A_0;
      count = 0;
      do
      {
        if( iter.face != start.face ){
          System.err.println("(edgeFlip) error: edge not linked to face "+count);
        }
        count++;
      } while( (iter = iter.next) != start);
      
      if( count != 3 ) System.err.println("(edgeFlip) error: to many edges");
      
      start = iter = edge_B_0;
      count = 0;
      do
      {
        if( iter.face != start.face ){
          System.err.println("(edgeFlip) error: edge not linked to face "+count);
        }
        count++;
      } while( (iter = iter.next) != start);
      
      if( count != 3 ) System.err.println("(edgeFlip) error: to many edges");

    }
    
    /**
     * inserts a new edge, connecting the two vertices of the given edges.<br>
     * make sure, the two edges are part of the same face (and are not connected already).<br>
     * the result is a new edge, that splits the original face.<br>
     * The two edges of the two vertices (vert_A.edge and vert_B.edge) must point
     * to the same face.<br>
     * After inserting the edge, the two vertice-edges point to the new edge.<br>
     * 
     * new created: edge_A, edge_B, face_B<br>
     * 
     * mesh.edge = edge_Av
     * vert_A.edge      = edge_A<br>
     * vert_A.edge.face = face_A<br>
     * vert_B.edge      = edge_B<br>
     * vert_B.edge.face = edge_B<br>
     * 
     * @param mesh
     * @param vert_A
     * @param vert_B
     */
    public static boolean insertEdge(HalfEdge.Mesh mesh, HalfEdge.Edge edge_A, HalfEdge.Edge edge_B){
      edge_A.orig.edge = edge_A;
      edge_B.orig.edge = edge_B;
      return insertEdge(mesh, edge_A.orig, edge_B.orig);
    }
    /**
     * inserts a new edge, connecting the two given vertices.<br>
     * make sure, the two vertices are part of the same face (and are not connected already).<br>
     * the result is a new edge, that splits the original face.<br>
     * The two edges of the two vertices (vert_A.edge and vert_B.edge) must point
     * to the same face.<br>
     * After inserting the edge, the two vertice-edges point to the new edge.<br>
     * 
     * new created: edge_A, edge_B, face_B<br>
     * 
     * mesh.edge = edge_Av
     * vert_A.edge      = edge_A<br>
     * vert_A.edge.face = face_A<br>
     * vert_B.edge      = edge_B<br>
     * vert_B.edge.face = edge_B<br>
     * 
     * @param mesh
     * @param vert_A
     * @param vert_B
     */
    public static boolean insertEdge(HalfEdge.Mesh mesh, HalfEdge.Vert vert_A, HalfEdge.Vert vert_B){
      // check if the two vertice edges point to the same face.
      if( vert_A.edge.face != vert_B.edge.face ){
        System.err.println("(insertEdge) vert_A.edge.face != vert_B.edge.face" );
        return false;
      }
      
      HalfEdge.Edge edge_A_old = vert_A.edge;
      HalfEdge.Edge edge_B_old = vert_B.edge;
      // check if the two vertice are already connected
      if( edge_A_old.next == edge_B_old ||edge_B_old.next == edge_A_old ){
        System.err.println("(insertEdge) edge_A.next == edge_B ||edge_B.next == edge_A" );
        return false;
      }
        
      HalfEdge.Face face_A     = edge_A_old.face;
      HalfEdge.Face face_B     = new HalfEdge.Face();
 
      HalfEdge.Edge edge_A_ins = new HalfEdge.Edge();
      HalfEdge.Edge edge_B_ins = new HalfEdge.Edge();
      
      // pairs
      edge_A_ins.pair = edge_B_ins;
      edge_B_ins.pair = edge_A_ins;
      
      // faces
      edge_A_ins.face = face_A;  face_A.edge = edge_A_ins;
      edge_B_ins.face = face_B;  face_B.edge = edge_B_ins;
      
      // origins
      edge_A_ins.orig = vert_A;   vert_A.edge = edge_A_ins;
      edge_B_ins.orig = vert_B;   vert_B.edge = edge_B_ins;

      HalfEdge.Edge iter;
      
      // Face B: next-links
      for(iter = edge_A_old; iter.next != edge_B_old; iter = iter.next){
        iter.face = face_B;
      }
      iter.face = face_B;
      iter.next = edge_B_ins;
      edge_B_ins.next = edge_A_old;
      
      // Face A (original face): next-links
      for(iter = edge_B_old; iter.next != edge_A_old; iter = iter.next){
//        iter.face = face_A; // already linked
      }
      iter.face = face_A;
      iter.next = edge_A_ins;
      edge_A_ins.next = edge_B_old;
      
      mesh.edge = edge_A_ins;

      HalfEdge.Update.face(face_A);
      HalfEdge.Update.face(face_B);
      return true;
    }
    
    
    /**
     * inserts an edge vertex / splits the given edge at param t.<br>
     * <br>
     * mesh.edge           = edge<br>
     * mesh.edge.next.orig = new vertex (given edge points towards the new vertex)<br>
     * 
     * @param mesh
     * @param edge the edge to split
     * @param where to split. t > 0 && t < 1
     */
    public static void insertVertex(HalfEdge.Mesh mesh, HalfEdge.Edge edge, float t){
      if( t <= 0 || t >= 1) return;
      float[] point =  DwVec3.lerp_new(edge.orig.data.v, edge.next.orig.data.v, t);
      insertVertex(mesh, edge, point);
    }
    
    /**
     * inserts an edge vertex / splits the given edge in half, using the edges'
     * mitpoint.<br>
     * <br>
     * mesh.edge           = edge<br>
     * mesh.edge.next.orig = edge-mitpoint (given edge points towards the new vertex)<br>
     * 
     * @param mesh
     * @param edge the edge to split
     */
    public static void insertVertex(HalfEdge.Mesh mesh, HalfEdge.Edge edge){
      float[] mid_point = DwVec3.add_new(edge.orig.data.v, edge.next.orig.data.v);
      DwVec3.scale_ref_slf(mid_point, 0.5f);
      insertVertex(mesh, edge, mid_point);
    }
    
    
    
    /**
     * inserts an edge vertex / splits the given edge in half, using the given 
     * point as the new vertex.<br>
     * 
     * mesh.edge           = edge<br>
     * mesh.edge.next.orig = edge-mitpoint (given edge points towards the new vertex)<br>
     * 
     * @param mesh
     * @param edge
     * @param point
     */
    public static void insertVertex(HalfEdge.Mesh mesh, HalfEdge.Edge edge, float[] point){

      
      HalfEdge.Vert vert_new = new HalfEdge.Vert(point);
      
      HalfEdge.Edge edge_L = edge;
      HalfEdge.Edge edge_R = edge.pair;
      
      HalfEdge.Edge edge_L_next = new HalfEdge.Edge();
      edge_L_next.face = edge_L.face;
      edge_L_next.orig = vert_new;
      edge_L_next.pair = edge_R;
      edge_L_next.next = edge_L.next;
      edge_R.pair = edge_L_next;
      edge_L.next = edge_L_next;
     
      HalfEdge.Edge edge_R_next = new HalfEdge.Edge();
      edge_R_next.face = edge_R.face;
      edge_R_next.orig = vert_new;
      edge_R_next.pair = edge_L;
      edge_R_next.next = edge_R.next;
      edge_L.pair = edge_R_next;
      edge_R.next = edge_R_next;
      
      vert_new.edge = edge_L_next;
      
      mesh.edge = edge;
    }
    
    /**
     * splits an edge of two adjacent faces (polygons, triangles, ...), using the given point as the new
     * vertex position.
     * the adjacent faces are triangulated from the new vertex, to every
     * other vertex, the polygon contains
     * 
     * @param mesh
     * @param point
     * @param edge
     * @return
     */
//    public static int edgeSplit(HalfEdge.Mesh mesh, HalfEdge.Edge edge, float[] point){
//      int num_new_faces = 0;
//      
//      HalfEdge.Vert vert_new = new HalfEdge.Vert(point);
//      
//      
//      // edge
//      vert_new.edge = edge;
//      
//      return num_new_faces;
//    }
    
    /**
     * triangulate a face, given by its adjacent edge.<br>
     * The triangulation is done from the edges origin, to all other face-verts.<br>
     * except the edge's preceding and next vertex (of course).<br>
     * in other words, the polygon becomes a triangle-fan.<br>
     * 
     * mesh.edge remains unchanged.<br>
     * 
     * <br>
     * in case the face has 3 verts, the face stays the same.<br>
     * in case the face has 4 verts, 1 edge is added.<br>
     * in case the face has 5 verts, 2 edges are added.<br>
     * 
     * 
     * @param edge
     * @return
     */
    
    public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Edge edge, List<HalfEdge.Face> faces_new){

      int num_new_faces = 0;
      
      final HalfEdge.Edge edge__ = edge;
      final HalfEdge.Face face__ = edge.face;
      final HalfEdge.Vert vert__ = edge.orig;
      
      HalfEdge.Edge edge_0 = edge__;
      HalfEdge.Edge edge_1 = edge_0.next;
      HalfEdge.Edge edge_2 = edge_1.next;
           
//      int counter = 0;
      while( edge_2.next.orig != vert__ ){

        HalfEdge.Face face_N = new HalfEdge.Face();  
        if( faces_new != null) faces_new.add(face_N); 
        num_new_faces++;
        
        HalfEdge.Edge edge_N = newEdge(edge_2.orig, face_N, vert__, face__);
        
//        edge_0.next = edge_1;
        edge_1.next = edge_N;
        edge_N.next = edge_0;
        
        edge_0.face = face_N;
        edge_1.face = face_N;
//        edge_N.face = face_N;

        if( !HalfEdge.Update.faceTriangle(face_N) ){
          //TODO
          System.err.println( "TRIANGLE IS DEGENERATE, cant comput normal");
        }
        edge_0 = edge_N.pair; 
        edge_1 = (edge_0.next = edge_2);
        edge_2 = edge_1.next;

//        counter++;
      }
      
//      System.out.println("counter = "+counter);
      if( num_new_faces > 0 ){
        edge_0.next.next.next = edge_0; // close triangle loop
        face__.edge = edge_0.next.next;
        
        if( !HalfEdge.Update.faceTriangle(face__) ){
          //TODO
          System.err.println( "TRIANGLE IS DEGENERATE, cant comput normal");
        }
      }
      
//      vert__.edge  = edge__;
//      edge__.orig = vert__;
//      mesh.edge = edge__;
      return num_new_faces;
    }
    
    
//    public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Edge edge, List<HalfEdge.Face> faces_new){
//
//      int num_new_faces = 0;
//      HalfEdge.Edge edge__ = edge;
//      HalfEdge.Face face__ = edge.face;
//      HalfEdge.Vert vert__ = edge.orig;
//           
//      int counter = 0;
//      while( edge.next.next.next.orig != vert__ ){
//        HalfEdge.Edge edge_0 = edge;
//        HalfEdge.Edge edge_1 = edge.next;
//        HalfEdge.Edge edge_2 = edge_1.next;
//  
//        HalfEdge.Face face_N = new HalfEdge.Face(); faces_new.add(face_N); num_new_faces++;
//        HalfEdge.Edge edge_N = newEdge(edge_2.orig, face_N, vert__, face__);
//        edge_1.next = edge_N;
//        edge_N.next = edge;
//        
//        edge_0.face = face_N;
//        edge_1.face = face_N;
//
//        HalfEdge.Update.face(face_N);
//        
//        edge = edge_N.pair;
//        edge.next = edge_2;
//
//        counter++;
//      }
//      
//      System.out.println("counter = "+counter);
//      if( num_new_faces > 0 ){
//        edge.next.next.next = edge; // close triangle loop
//        face__.edge = edge.next.next;
//        HalfEdge.Update.face(face__);
//      }
//      
//      mesh.edge = edge__;
//      return num_new_faces;
//    }
    
   
    private static HalfEdge.Edge newEdge(
        HalfEdge.Vert vert_L, HalfEdge.Face face_L,
        HalfEdge.Vert vert_R, HalfEdge.Face face_R){
      
      HalfEdge.Edge edge_L = new HalfEdge.Edge();
      HalfEdge.Edge edge_R = new HalfEdge.Edge();
      
      edge_L.pair = edge_R;
      edge_R.pair = edge_L;
      
      edge_L.face = face_L; face_L.edge = edge_L;
      edge_R.face = face_R; face_R.edge = edge_R;
      
      edge_L.orig = vert_L; vert_L.edge = edge_L;
      edge_R.orig = vert_R; vert_R.edge = edge_R;
      
      return edge_L;
    }
    
    
    
    
    
    
//   public static int triangulateFace(HalfEdge.Mesh mesh, HalfEdge.Edge edge, List<HalfEdge.Face> faces_new){
//      
//      if( edge.next.next.next == edge ){
//        System.out.println("(triangulateFace(HalfEdge.Edge edge) - face is triangle");
//        return 0; // face is a triangle, so no triangulation makes sense
//      }
//      
//      int num_new_faces = 0;
//      
//      // keep start values
////      HalfEdge.Edge edge_start = edge;
//      HalfEdge.Vert vert__ = edge.orig;
//      HalfEdge.Face face_A = edge.face;
//      
//      // get last edge
//      HalfEdge.Edge last = edge;
//      while( last.next != edge ) last = last.next;
//      
//      HalfEdge.Edge edge_0 = edge;
//      HalfEdge.Edge edge_1 = edge.next;
//      HalfEdge.Edge edge_2 = edge.next.next;
//      
//      // new face
//      HalfEdge.Face face_B = new HalfEdge.Face();
//      faces_new.add(face_B);
//      num_new_faces++;
//
//      // new edge
//      HalfEdge.Edge edge_2_new = new HalfEdge.Edge();
//      edge_2_new.next = edge_0;
//      edge_2_new.orig = edge_2.orig;
//      
//      face_B.edge = edge_0;
//      edge_0.face = face_B;
//      edge_1.face = face_B;
//      edge_2_new.face = face_B;
//      
//      edge_0.next = edge_1;
//      edge_1.next = edge_2_new;
//      edge_2_new.next = edge_0;
//     
//      // new pair edge == start edge for next iteration
//      HalfEdge.Edge edge_2_new_pair = new HalfEdge.Edge();
//      edge_2_new_pair.next = edge_2;
//      edge_2_new_pair.face = face_A;
//      edge_2_new_pair.orig = vert__;
//      face_A.edge = edge_2_new_pair;
//      
//      edge_2_new.pair = edge_2_new_pair;
//      edge_2_new_pair.pair = edge_2_new;
//      
//      last.next = edge_2_new_pair;
//      
//      //TODO: loop
//      
//      HalfEdge.Update.face(face_B);
//      
//      HalfEdge.Update.face(face_A);
//      
//      mesh.edge = edge_2_new_pair;
//   
//      return num_new_faces;
//    }
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
    
    
    
    
  
  
  
    
    
  
  
  
    
    
  //============================================================================
  //
  // CREATE
  //
  //============================================================================
    
    
  public abstract static class Creator{
    
    public static boolean ERROR_CHECKING = true;
    
    
    static public IFS convertToIFS(HalfEdge.Mesh mesh){

      // collecting half-edge stuff
      final HalfEdge.Vert[] verts = HalfEdge.Collect.verts(null, mesh);
      final HalfEdge.Face[] faces = HalfEdge.Collect.faces(null, mesh);

      final int num_verts = verts.length;
      final int num_faces = faces.length;
      
      // creating IFS-structure
      final IFS         ifs   = new IFS();
      final IFS.Group[] ifs_g = { new IFS.Group(ifs, "group_halfedge_mesh") };
      final IFS.Mesh[]  ifs_m = { new IFS.Mesh (ifs_g[0], mesh.data.name) };
      final IFS.Face[]  ifs_f = new IFS.Face[num_faces];
      final float[][]   ifs_v = new float[num_verts][];
      final float[][]   ifs_t = null; // TODO
      final float[][]   ifs_n = null; // TODO
     
      // vertices
      for(int i = 0; i < num_verts; i++){
        ifs_v[i] = verts[i].data.v;
      }
      
      // faces
      int[] idx_v_buffer = new int[3];
      for(int i = 0; i < num_faces; i++){
        HalfEdge.Face face = faces[i];

        // get face vertices
        int v_pos = 0;
        HalfEdge.Edge edge = face.edge;
        do {   
          if( v_pos >= idx_v_buffer.length){
            idx_v_buffer = Arrays.copyOf(idx_v_buffer, v_pos<<1);
          }
          idx_v_buffer[v_pos++] = edge.orig.data.id;
        } while( (edge=edge.next) != face.edge);
        
        final int[] idx_v = Arrays.copyOf(idx_v_buffer, v_pos);
        final int[] idx_n = null; // TODO
        final int[] idx_t = null; // TODO
         
        ifs_f[i] = new IFS.Face(ifs_m[0], idx_v, idx_n, idx_t);
      }
      
      ifs_g[0].m.add(ifs_m[0]);
      ifs_m[0].f.addAll(Arrays.asList(ifs_f));
      
      return ifs.set( ifs_g, ifs_m ,ifs_f, ifs_v, ifs_t, ifs_n );
    }
    
    
    
    static public HalfEdge.Mesh convertFromIFS(IFS ifs){
      //TODO handle incorrect ifs-input
      int num_verts = ifs.v.length;
      int num_faces = ifs.f.length;
      int num_edges = 0; // depends of the number of vertices per face.
      
      HalfEdge.Vert[] he_verts;
      HalfEdge.Face[] he_faces;
      HalfEdge.Edge[] he_edges;
      
      
      //------------------------------------------------------------------------
      // 1) INIT VERTICES
      he_verts = new HalfEdge.Vert[num_verts];
      for(int i = 0; i < num_verts; i++){
        he_verts[i]         = new HalfEdge.Vert(ifs.v[i]);
//        he_verts[i].data.id = i;
//        he_verts[i].edge    = null; // done in step 3
      }
      
      
      //------------------------------------------------------------------------
      // 2) INIT FACES ... and compute number of edges
      he_faces = new HalfEdge.Face[num_faces];
      for(int i = 0; i < num_faces; i++){
        he_faces[i]      = new HalfEdge.Face();
//        he_faces[i].edge = null; // done in step 3
        num_edges += ifs.f[i].IDX_V.length;
      }

      
      //------------------------------------------------------------------------
      // 3) INIT EDGES ... link to face and edge.next
      he_edges = new HalfEdge.Edge[num_edges];
      for(int id_f = 0, id_e = 0; id_f < num_faces; id_f++){
        
        int[] vert_indices = ifs.f[id_f].IDX_V;

        int num_face_edges = vert_indices.length;
        int id_e_cur = id_e;                  // first edge index of current face
        int id_e_nxt = id_e+num_face_edges;   // first edge index of next face
        
        // init edges for current face
        for(int id_v = 0; id_v < num_face_edges; id_v++, id_e++){
          he_edges[id_e]           = new HalfEdge.Edge();
          he_edges[id_e].orig      = he_verts[vert_indices[id_v]];
          he_edges[id_e].face      = he_faces[id_f];
          he_edges[id_e].orig.edge = he_edges[id_e]; // link vertex with edge
          he_edges[id_e].face.edge = he_edges[id_e]; // link face with edge
          he_edges[id_e].next      = null; // done in next step
          he_edges[id_e].pair      = null; // done, after all edges/faces are created
        }
        
        // link edges (edge.next) for current face in CCW-Direction
        he_edges[id_e_nxt-1].next = he_edges[id_e_cur];
        while( id_e_cur < id_e_nxt-1) {
          he_edges[id_e_cur].next = he_edges[++id_e_cur];
        }
      }
      
      

      //------------------------------------------------------------------------
      // 4) COMPLETE HE-STRUCTURE (pair-connection)
      
//      long timer;
      
//      timer = System.currentTimeMillis();
//      HalfEdge.Edge[] e = he_edges;
//      HalfEdge.Edge ei, ej;
//      for( int i = 0; i < num_edges; i++ ){
//        if( (ei=e[i]).pair != null ) continue;
//        
//        HalfEdge.Vert ei_v0 = ei.orig;
//        HalfEdge.Vert ei_v1 = ei.next.orig;
//        
//        for( int j = i+1; j < num_edges ; j++ ){
//          if( (ej=e[j]).pair != null ) continue;
//
//          if( ei_v0 == ej.next.orig && ei_v1 == ej.orig) {
//            ei.pair = ej;
//            ej.pair = ei;
//            break;
//          }
//        }
//      }
//      timer = System.currentTimeMillis()-timer;
//      System.out.println("  linking pairs: "+timer+" ms");

      
//      timer = System.currentTimeMillis();
      // setup a hash-map of all edges, using the edge's vertices (A->B) as key.
      HashMap<EdgeVerts, HalfEdge.Edge> map = new HashMap<EdgeVerts, HalfEdge.Edge>(num_edges);
      for( int i = 0; i < num_edges; i++ ){
        if( map.put( EdgeVerts.AB(he_edges[i]), he_edges[i] ) != null ){
          System.err.println("HashMap-Error during map creation (hashcode collision)");
        }
      }
      
      // to find an edge-pair, the map is queried by using the reversed vertex order (B->A) as key.
      // this is really fast!
      for( int i = 0; i < num_edges; i++ ){
        if( he_edges[i].pair == null ) {
          HalfEdge.Edge e1 = he_edges[i];                // find pair for e1
          HalfEdge.Edge e2 = map.get( EdgeVerts.BA(e1)); // search map, key = swapped verts
          if( e2 != null ){ // on closed meshes, edge2 is never null!!
            e1.pair = e2;
            e2.pair = e1; 
          } else {
            System.err.println("HashMap-Error during pair-linking (hashcode collision)");
          }
        }
      }
//      timer = System.currentTimeMillis()-timer;
//      System.out.println("  linking pairs: "+timer+" ms");

      //------------------------------------------------------------------------
      // 5) DEBUGGING: CHECK PAIRS
      for( int i = 0; ERROR_CHECKING && i < num_edges; i++ ){
        HalfEdge.Edge edge = he_edges[i];
        
        if( edge.pair == null ) {
          System.err.println("WARNING: edge.pair == null  "+i);
          continue;
        }
        if( edge == edge.pair ){
          System.err.println("ERROR: edge.pair == edge  "+i);
          continue;
        }
        if( edge != edge.pair.pair) {
          System.err.println("ERROR: edge.pair.pair != edge  "+i);
          continue;
        }
        Vert edge_a = edge.orig;
        Vert edge_b = edge.next.orig;
        Vert pair_a = edge.pair.orig;
        Vert pair_b = edge.pair.next.orig;
        if( edge_a != pair_b || edge_b != pair_a){
          System.err.println("ERROR:  edge_a != pair_b || edge_b != pair_a  "+i);
          continue;
        }
      }
      
      //------------------------------------------------------------------------
      // 6) INIT MESH
      HalfEdge.Mesh he_mesh = new HalfEdge.Mesh();
      he_mesh.edge = he_edges[0];
      
      return he_mesh;
    }
    
    
    
    
    
    /**
     * save the 2 vertices of an edge.<br>
     * is used during pairwise linking of half-edges, as the key for the hashmap,<br>
     * to find an edge's pair.<br>
     * 
     * http://stackoverflow.com/questions/7032961/java-how-to-use-a-pair-of-keys-for-hashmap
     * 
     * @author thomas diewald
     */
    private static final class EdgeVerts {
      
      private final Vert a, b;
      private final int hashcode;
      
      public static final EdgeVerts AB(final Edge e){  return new EdgeVerts( e.orig, e.next.orig ); }
      public static final EdgeVerts BA(final Edge e){  return new EdgeVerts( e.next.orig, e.orig ); }

      private EdgeVerts(final Vert a, final Vert b) { 
        this.a = a; 
        this.b = b; 
        this.hashcode = a.hashCode() + b.hashCode()*31;
      }
      
      @Override  
      public final int hashCode() { return hashcode; }
      
      @Override
      public final boolean equals(Object o) {
        return equals( (EdgeVerts)o );
      }
      public final boolean equals(final EdgeVerts e){
        return  (a.equals(e.a) && b.equals(e.b) );
      }
    }

    
    
    
    
    
    
    
    
    
    
    
    
    
//    /**
//     * creates a new mesh (triangle, flat!!).
//     * @param p0
//     * @param p1
//     * @param p2
//     * @return new created mesh
//     */
//    public static HalfEdge.Mesh newTriangle(float[] p0, float[] p1, float[] p2){
//      HalfEdge.Face face = newFace(p0, p1, p2);
//      HalfEdge.Face oppo = addOpposite(face);
//      HalfEdge.Mesh mesh = new HalfEdge.Mesh();
//      mesh.edge = face.edge;
//      HalfEdge.Update.verts(mesh);
//      // TODO: add mesh face count
//      return mesh;
//    }
    
    

    

    /**
     * <br>
     * creates a new mesh from a given number of base points and one apex.<br>
     * <br>
     * base.length  = 2 --> triangle<br>
     * base.length  = 3 --> tetrahedron<br>
     * base.length  = 4 --> square pyramid (in case base is planar)<br>
     * base.length  = 5 --> pentagonal pyramid (in case base is planar)<br>
     * (http://en.wikipedia.org/wiki/Pyramid_%28geometry%29)<br>
     * <br>
     * points are not tested for being co-planar/linear/... .<br>
     * <br>
     * mesh.edge.orig = apex<br>
     * </p>
     * @param base  number of n base points:  float[n][3]
     * @param apex  apex
     * @return new created mesh
     */
    public static HalfEdge.Mesh newPyramid(List<HalfEdge.Face> faces_new, float[] apex, float[] ... base){
      HalfEdge.Mesh mesh = new HalfEdge.Mesh();
      HalfEdge.Face face = newFace(base);
      HalfEdge.Face oppo = addOpposite(face);
      
      // make sure, the face, pointing to apex, gets triangulated.
      if( face.data.distance(apex) < 0 ) 
        face = face.edge.pair.face; //opposite face
      
      if( faces_new != null ) faces_new.add(face.edge.pair.face);
      HalfEdge.MeshOps.triangulateFace(mesh, face, apex, faces_new);
//
      HalfEdge.Update.verts(mesh);
      // TODO: add mesh face count
      return mesh;
    }
    
    
    /**
     * creates a new face (polygon) from the given list of points.<br>
     * the face-vertices are in the same order as the given points are.<br>
     * so, the resulting face might not be planar.<br>
     * to make a valid mesh (pair-links) the opposite face has to be generated<br>
     * separately by calling HalfEdge.Creator.createOpposite(...)<br>
     * 
     * @param points
     * @return
     */
    public static HalfEdge.Face newFace(float[] ... points){
      
      HalfEdge.Face face = new HalfEdge.Face();
     
      HalfEdge.Edge iter = new HalfEdge.Edge();
      iter.orig          = new HalfEdge.Vert( points[0] );
      iter.orig.edge     = iter;
      iter.face          = face;
      iter.face.edge     = iter; 
      iter.pair          = null;

      HalfEdge.Edge last = iter;

      for( int i = 1; i < points.length; i++, iter=iter.next){
        (iter.next)           = new HalfEdge.Edge();
        (iter.next).orig      = new HalfEdge.Vert( points[i] );
        (iter.next).orig.edge = iter.next;
        (iter.next).face      = face; 
        (iter.next).pair      = null;
        (iter.next).next      = last; // close edge-loop
      }
      HalfEdge.Update.face(face);
        
      return face;
    }
    
    
    /**
     * create the opposite/back face of a given face.<br>
     * be aware, all possible former links to other faces will be destroyed.<br>
     * <br>
     * so this can also be used to disconnect a polygon/face from a mesh.<br>
     * the old mesh will be broken then.<br>
     * <br>
     * TODO: fill holes, repair mesh<br>
     * 
     * 
     * @param face 
     * @return the opposite face
     */
    public static HalfEdge.Face addOpposite(HalfEdge.Face face){
      
      HalfEdge.Face oppo = new HalfEdge.Face();
      
      HalfEdge.Edge iter = face.edge;

      // create pair-edges, set their origins and face, and make pairs-links
//      int counter = 0;
      do
      {
        HalfEdge.Edge pair = new HalfEdge.Edge();
        pair.face = oppo;
        oppo.edge = pair;
        pair.orig = iter.next.orig;
        pair.pair = iter;
        iter.pair = pair;
//        counter++;
      } while( (iter=iter.next) != face.edge);
//      System.out.println(counter);
      
      
      // link new edges amongst each other, make next-links
      do
      {
        iter.next.pair.next = iter.pair;
      } while( (iter=iter.next) != face.edge);
      
      
      HalfEdge.Update.face(oppo);
      return oppo;
    }
    

    
    
  }
  
  
  
  
  
  
  
  

  
  
}