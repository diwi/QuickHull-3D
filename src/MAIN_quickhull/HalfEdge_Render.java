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
import java.util.List;

import MAIN_quickhull.HalfEdge.Query;
//import MAIN_quickhull.HalfEdge.Selection;

import processing.core.PConstants;
import processing.core.PGraphics3D;



public class HalfEdge_Render {
  //============================================================================
  //
  // DISPLAY/RENDER
  //
  //============================================================================
  

//    public static IFS ifs;
    
//    static public void faces(PGraphics3D g, HalfEdge.Mesh mesh){
//      List<HalfEdge.Face> list = HalfEdge.Collect.faces_List(null, mesh);
//      List<HalfEdge.Face> list = HalfEdge.Collect.faces_old(null, mesh);
//      for( HalfEdge.Face item : list ) HalfEdge_Render.face(g, item);
//    } 
    
    // for faster collecting
    private static HalfEdge.Face[] face_list;
    private static HalfEdge.Edge[] edge_list;
    private static HalfEdge.Vert[] vert_list;
    
    static public void faces(final PGraphics3D g, final HalfEdge.Mesh mesh){
      face_list = HalfEdge.Collect.faces(face_list, mesh);
      for( HalfEdge.Face face : face_list ){ 
        HalfEdge.Edge edge = face.edge;
        if( edge == edge.next.next.next ) 
          HalfEdge_Render.faceTriangle(g, face);
        else if( edge == edge.next.next.next.next ) 
          HalfEdge_Render.faceQuad(g, face);
        else
          HalfEdge_Render.face(g, face);
      }
    } 
    static public void edges(final PGraphics3D g, final HalfEdge.Mesh mesh){
      edge_list = HalfEdge.Collect.edges(edge_list, mesh);
      for( HalfEdge.Edge edge : edge_list ) HalfEdge_Render.edge(g, edge);
    }
    static public void verts(final PGraphics3D g, final HalfEdge.Mesh mesh){
      vert_list = HalfEdge.Collect.verts(vert_list, mesh);
      for( HalfEdge.Vert vert : vert_list ) HalfEdge_Render.vert(g, vert);
    }
    
    
    
    
    static public void faceNormals(final PGraphics3D g, final HalfEdge.Mesh mesh, float scale){
      face_list = HalfEdge.Collect.faces(face_list, mesh);
      for( HalfEdge.Face face : face_list ) HalfEdge_Render.faceNormal(g, face, scale);
    } 
    static public void faceNormal(final PGraphics3D g, final HalfEdge.Face face, float scale){
      g.stroke(255,0,0); new DwRay3D().start(face.data.c).end(face.edge.orig.data.v     ).draw(g, .50f);  // X
      g.stroke(0,255,0); new DwRay3D().start(face.data.c).end(face.edge.next.orig.data.v).draw(g, .50f);  // Y
      g.stroke(0,0,255); new DwRay3D().start(face.data.c).dir(face.data.n               ).draw(g, scale); // Z = normal
    } 
    
    static public void vertNormals(final PGraphics3D g, final HalfEdge.Mesh mesh, float scale){
      vert_list = HalfEdge.Collect.verts(vert_list, mesh);
      for( HalfEdge.Vert vert : vert_list ) HalfEdge_Render.vertNormal(g, vert, scale);
    } 
    static public void vertNormal(final PGraphics3D g, final HalfEdge.Vert vert, float scale){
      g.stroke(0,0,255); new DwRay3D().start(vert.data.v).dir(vert.data.n).draw(g, scale); // Z = normal
    } 
    
    static public void halfEdges(final PGraphics3D g, final HalfEdge.Mesh mesh){
      face_list = HalfEdge.Collect.faces(face_list, mesh);
      for( HalfEdge.Face face : face_list ) HalfEdge_Render.halfEdges(g, face);
    } 
    static public void halfEdges(final PGraphics3D g, final HalfEdge.Face face){
      List<HalfEdge.Edge> list = new ArrayList<HalfEdge.Edge>();
      Query.adjacentEdges(list, face);

      DwRay3D ray_he = new DwRay3D();
      DwRay3D ray_v0 = new DwRay3D();
      DwRay3D ray_v1 = new DwRay3D();
      float[] v0 = new float[3]; // start of current edge
      float[] v1 = new float[3]; // start of next edge
      float[] v2 = new float[3]; // end of current edge
      
      g.beginShape(PConstants.LINES);
      for(HalfEdge.Edge item : list){
        ray_v0.start(face.data.c).end(item.orig.data.v     ).getPoint(.9f, v0);
        ray_v1.start(face.data.c).end(item.next.orig.data.v).getPoint(.9f, v1);
        ray_he.start(v0).end(v1).getPoint(.9f, v2);
//        g.stroke(0);
        g.stroke(255,0,0);
        vertexArr(g, v0);
//        g.stroke(255,100,0);
        g.stroke(0,255,0);
        vertexArr(g, v2);
      }
      g.endShape();
    } 
    
    
    
    
    static public void halfEdgesPairing(final PGraphics3D g, final HalfEdge.Mesh mesh){
      // face-edge links are always valid, so therefore:
      // 1) collect all faces ...
      face_list = HalfEdge.Collect.faces(face_list, mesh);
      // 2) for each face ...
      for( HalfEdge.Face face : face_list ) {
        // 3) get it's edges ...
        List<HalfEdge.Edge> edges = new ArrayList<HalfEdge.Edge>();
        Query.adjacentEdges(edges, face);
        // 4) for each edge ...
        for( HalfEdge.Edge edge : edges ) {
          // 5) display parent edge.
          HalfEdge_Render.halfEdgesPairing(g, edge);
        }
      }
    } 
    static public void halfEdgesPairing(final PGraphics3D g, final HalfEdge.Edge edge){
      if( edge.pair == null ){
        System.err.println("(HalfEdge.Display.halfEdgesPairing()): ERROR -> edge got pair");
        return;
      }
      
      float[] face_center = edge.face.data.c;
      float[] edge_center = DwVec3.add_new(edge     .orig.data.v, edge.next     .orig.data.v); 
      float[] pair_center = DwVec3.add_new(edge.pair.orig.data.v, edge.pair.next.orig.data.v); 
      DwVec3.scale_ref_slf(pair_center, .5f);
      DwVec3.scale_ref_slf(edge_center, .5f);
      
      DwRay3D ray_c_2_e = new DwRay3D().start(face_center             ).end(edge_center);
      DwRay3D ray_e_2_p = new DwRay3D().start(ray_c_2_e.getPoint(0.5f)).end(pair_center);
      
      g.beginShape(PConstants.LINE);
      {
        g.stroke(0); vertexArr(g,  ray_c_2_e.o);
        g.stroke(0); vertexArr(g,  ray_c_2_e.getPoint(0.66f));
        
        g.stroke(0);   vertexArr(g,  ray_c_2_e.getPoint(0.66f));
        g.stroke(255); vertexArr(g,  ray_e_2_p.getPoint(1));
      }
      g.endShape();
    } 
    
    
   
    static public void selection(PGraphics3D g, HalfEdge.Selection.Face selection) {
      g.noStroke();

      float col_step = 255f/selection.mark;
      for( HalfEdge.Face item : selection.selection ){
        g.fill(0,0, col_step*item.flag.selection);
        HalfEdge_Render.face(g, item);
      }
    }
    
    static public void selection(PGraphics3D g,  HalfEdge.Selection.Edge selection) {
//      g.strokeWeight(2);
      float col_step = 255f/selection.mark;
      for( HalfEdge.Edge item : selection.selection ){
        g.stroke(col_step*item.flag.selection,0,0);
        HalfEdge_Render.edge(g, item);
      }
    }
    
    static public void selection(PGraphics3D g,  HalfEdge.Selection.Vert selection) {
//      g.strokeWeight(5);
      float col_step = 255f/selection.mark;
      for( HalfEdge.Vert item : selection.selection ){
        g.stroke(0,col_step*item.flag.selection,0);
        HalfEdge_Render.vert(g, item);
      }
    }
    
    static public void faceTriangle(final PGraphics3D g, final HalfEdge.Face face){
      HalfEdge.Edge e = face.edge;
      g.beginShape(PConstants.TRIANGLES);
      {
//        normalArr(g, e.orig.data.n);
        vertexArr(g, e.orig.data.v); e=e.next;
//        normalArr(g, e.orig.data.n);
        vertexArr(g, e.orig.data.v); e=e.next;
//        normalArr(g, e.orig.data.n);
        vertexArr(g, e.orig.data.v);
      }
      g.endShape();
    }
    
    static public void faceQuad(final PGraphics3D g, final HalfEdge.Face face){
      HalfEdge.Edge e = face.edge;
      g.beginShape(PConstants.QUADS);
      {
        vertexArr(g, e.orig.data.v); e=e.next;
        vertexArr(g, e.orig.data.v); e=e.next;
        vertexArr(g, e.orig.data.v); e=e.next;
        vertexArr(g, e.orig.data.v);
      }
      g.endShape();
    }
    static public void face(final PGraphics3D g, final HalfEdge.Face face){
//      g.beginShape(PConstants.POLYGON);
//      {
//        HalfEdge.Edge e = face.edge;
//        do{
//          vertexArr(g, e.orig.data.v);
//        } while( (e = e.next) != face.edge);
//      }
//      g.endShape();
      
      g.beginShape(PConstants.TRIANGLE_FAN);
      {
        HalfEdge.Edge e = face.edge;
        vertexArr(g, e.orig.data.v);
        do{
          vertexArr(g, e.orig.data.v);
        } while( (e = e.next) != face.edge);
      }
      g.endShape();
    }
    
    static public void edge(final PGraphics3D g, final HalfEdge.Edge edge){
      g.beginShape(PConstants.LINES);
      {
        vertexArr(g, edge.orig.data.v);
        vertexArr(g, edge.next.orig.data.v);
      }
      g.endShape();
    }
    
    static public void vert(final PGraphics3D g, final HalfEdge.Vert vert){
      g.beginShape(PConstants.POINTS);
      {
        vertexArr(g, vert.data.v);
      }
      g.endShape();
    }
    
    static public void vertexArr(final PGraphics3D g, final float[] v){
      g.vertex(v[0], v[1], v[2]);
    }
    static public void normalArr(final PGraphics3D g, final float[] v){
      g.normal(v[0], v[1], v[2]);
    }

}
