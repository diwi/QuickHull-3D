/**
 * 
 *   author: (c) Thomas Diewald, http://thomasdiewald.com
 *   date: 23.04.2012
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

import processing.core.PGraphics3D;


public class DwRay3D {
  public float[] o, d, d_rec;

  public DwRay3D(){
  }
  
  public DwRay3D(float[] o, float[] d){
    start(o).dir(d);
  }
  
  public DwRay3D start(float[] startpoint){
    o = startpoint;
    return this;
  }
  
  public DwRay3D end(float[] endpoint){
    d = DwVec3.sub_new(endpoint, o);
    d_rec = new float[]{1f/d[0], 1f/d[1], 1f/d[2]};
    return this;
  }
  
  public DwRay3D dir(float[] direction){
    this.d = direction;
    this.d_rec = new float[]{1f/d[0], 1f/d[1], 1f/d[2]};
    return this;
  }
 
  public DwRay3D copy(){
    return new DwRay3D(DwVec3.copy_new(o), DwVec3.copy_new(d));
  }

  public DwRay3D normalize(){
    DwVec3.normalize_ref_slf(d);
    dir(d);
    return this;
  }
  
  public float[] getPoint(float t){
    return DwVec3.add_new(o, DwVec3.scale_new(d, t) );
  }
  
  public DwRay3D getPoint(float t, float[] dst){
    DwVec3.add_ref(o, DwVec3.scale_new(d, t), dst);
    return this;
  }

  public DwRay3D draw(PGraphics3D g){
    float[] e = getPoint(1);
    g.line(o[0], o[1], o[2], e[0], e[1], e[2]);
    return this;
  }
  
  public DwRay3D draw(PGraphics3D g, float t){
    float[] e = getPoint(t);
    g.line(o[0], o[1], o[2], e[0], e[1], e[2]);
    return this;
  }
}
