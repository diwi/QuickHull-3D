/**
 * 
 *   author: (c) Thomas Diewald, http://thomasdiewald.com
 *   date: 21.01.2013
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



import processing.core.PConstants;
import processing.core.PGraphics3D;


public class IFS_Render {
  static private PGraphics3D g;
  static int id = 0;
  
  static public void g(PGraphics3D g){
    IFS_Render.g = g;
  }
  
  static public void render(IFS ifs){
    for(IFS.Face f : ifs.f){
      render(ifs, f);
    }
  }
  static public void render(IFS.Mesh ifs_mesh){
    IFS ifs = ifs_mesh.parent.parent;
    for(IFS.Face f : ifs_mesh.f){
      render(ifs, f);
    }
  }
  static public void render(IFS.Group ifs_group){
    IFS ifs = ifs_group.parent;
    for(IFS.Mesh m : ifs_group.m ){
      for(IFS.Face f : m.f){
        render(ifs, f);
      }
    }
  }
  
  static public void renderVerts(IFS ifs){
    float[][] verts = ifs.v;
    g.beginShape(PConstants.POINTS);
      for(float[] v : verts )
        vertex(v);
    g.endShape();

  }

  

  // TODO: currrently not optimized for speed!!!!
  static public void render(IFS ifs, IFS.Face f){
    int num_verts = f.IDX_V.length;
   
    if( num_verts == 2 ) { renderLine    (ifs, f); return;  };
    if( num_verts == 3 ) { renderTriangle(ifs, f); return;  };
    if( num_verts == 4 ) { renderQuad    (ifs, f); return;  };
    if( num_verts >  4 ) { renderPolygon (ifs, f); return;  };
  }
  
  
  static public void renderLine(IFS ifs, IFS.Face f){
    g.beginShape(PConstants.LINE);
    {
      vertex(ifs.v[ f.IDX_V[0] ]);
      vertex(ifs.v[ f.IDX_V[1] ]);
    }
    g.endShape();
  }
  
  static public void renderTriangle(IFS ifs, IFS.Face f){
    g.beginShape(PConstants.TRIANGLES);
    {
      vertex(ifs.v[ f.IDX_V[0] ]);
      vertex(ifs.v[ f.IDX_V[1] ]);
      vertex(ifs.v[ f.IDX_V[2] ]);
    }
    g.endShape();
  }
  static public void renderQuad(IFS ifs, IFS.Face f){
    g.beginShape(PConstants.QUADS);
    {
      vertex(ifs.v[ f.IDX_V[0] ]);
      vertex(ifs.v[ f.IDX_V[1] ]);
      vertex(ifs.v[ f.IDX_V[2] ]);
      vertex(ifs.v[ f.IDX_V[3] ]);        
    }
    g.endShape();
  }
  static public void renderPolygon(IFS ifs, IFS.Face f){
//    g.beginShape(PConstants.POLYGON);
//    {
//      for( int i = 0; i < f.IDX_V.length; i++){
//        vertex(ifs.v[ f.IDX_V[i] ]);
//      }      
//    }
//    g.endShape(PConstants.CLOSE);
    
    g.beginShape(PConstants.TRIANGLE_FAN);
    {
      vertex(ifs.v[0]);
      for(int i : f.IDX_V){
        vertex(ifs.v[i]);
      }
    }
    g.endShape();
  }
  
  static public void vertex(float[] v){
    g.vertex(v[0], v[1], v[2]);
  }
}
