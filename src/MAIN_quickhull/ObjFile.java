/**"+
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


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;



public class ObjFile {

  public static boolean DEBUG = false;

  /**
   * obj-file parser:
   * parses only vertices, texture-coordinated, normals and faces
   * TODO: materials, ...
   * 
   * http://en.wikipedia.org/wiki/Wavefront_.obj_file
   * 
   * @param lines obj-file content (line separated)
   * @return Index Face Set (IFS)
   */
  public static IFS parse(String[] lines){
    IFS ifs = new IFS();
    
    // tmp data buffers
    ArrayList<float[]>  tmp_v = new ArrayList<float[]>();
    ArrayList<float[]>  tmp_n = new ArrayList<float[]>();
    ArrayList<float[]>  tmp_t = new ArrayList<float[]>();
    ArrayList<IFS.Face> tmp_f = new ArrayList<IFS.Face>();
    
    // TODO: might better to keep the HashMaps than the array in the IFS class!
    HashMap<String, IFS.Group>tmp_g = new HashMap<String, IFS.Group>();
    HashMap<String, IFS.Mesh> tmp_m = new HashMap<String, IFS.Mesh>();
    
    IFS.Group g_curr = null;
    IFS.Mesh  m_curr = null;
    
    // buffer to save the (usually) four split-tokens per line.
    String[] token = new String[4];
    
    
    for(int i = 0; i < lines.length; i++){
//      System.out.println("["+i+"] "+lines[i]);
      String line = lines[i]; 
      
      // skip line-comment
      if( line.indexOf('#') == 0){
//        System.out.println(line);
        continue;
      }

      // rhino export breaks too long lines, so this tries to rebuild the line as it should be.
      while( line.endsWith("\\") ){
        line = line.substring(0, line.length()-1) +" "+ lines[++i]; // delete '\' from current line, and add next line
      }
      
      
      // split token (slow version), and also skips multiple white-spaces, so no empty tokens are generated
//      token = line.split("\\s+"); 
//      int NUM_TOKENS = token.length; 
      
      
      // tokenize line, skip Whitespaces
      int NUM_TOKENS = 0;
      int a = -1, b = -1, end = line.length();
      while( ++a != end ){
        if( ( b = line.indexOf(' ', a))  == -1) {
          token[NUM_TOKENS++] = line.substring(a);
          break;
        } else if( a != b){
          token[NUM_TOKENS++] = line.substring(a, b);
          a = b;
        } 
      }

      if( NUM_TOKENS == 0 ) continue;
      
      String ID = token[0];
      
      // GROUP/LAYER
      if( ID.equals("g") ){
        g_curr = tmp_g.get(token[1]);
        if( g_curr == null){
          g_curr = new IFS.Group(ifs, token[1]);
          tmp_g.put(g_curr.name, g_curr );
        }
      }
      // OBJECT/MESH
      else if( ID.equals("o") ){
        m_curr = tmp_m.get(token[1]);
        
        if( m_curr == null ){
          if( g_curr == null ){
            g_curr = new IFS.Group(ifs, "default_group");
            tmp_g.put(g_curr.name, g_curr );
          }
          m_curr = new IFS.Mesh(g_curr, token[1]);
          g_curr.m.add(m_curr);
          tmp_m.put(m_curr.name, m_curr );
        }
      }
      // VERTICES (float[3])
      else if( ID.equals("v") ){
        float vx = Float.parseFloat( token[1] );
        float vy = Float.parseFloat( token[2] );
        float vz = Float.parseFloat( token[3] );
        tmp_v.add( new float[]{vx, vy, vz} );
      }
      // TEXTURE COORDINATES (float[2])
      else if( ID.equals("vt")){
        float u = Float.parseFloat( token[1] );
        float v = Float.parseFloat( token[2] );
        tmp_t.add( new float[]{u, v} );
      }
      // VERTEX NORMALS (float[3])
      else if( ID.equals("vn") ){
        float nx = Float.parseFloat( token[1] );
        float ny = Float.parseFloat( token[2] );
        float nz = Float.parseFloat( token[3] );
        tmp_n.add( new float[]{nx, ny, nz} );
      }
      // FACES (indices for v/t/n ... only v is mandatory ... indices start at 1!!!)
      else if( ID.equals("f") ){
        IFS.Face f = new IFS.Face( NUM_TOKENS-1 );

        for(int j = 1, id = 0; j < NUM_TOKENS; j++, id++){
          // obj-indices start at 1, so i initialize them with -1 and add the parsed value.
          // -1 also indicate: no data present
          f.IDX_V[id] = f.IDX_T[id] = f.IDX_N[id] = -1; 
          
          // face format: v   v/t   v//n    v/t/n
          // triangle example: f v/t/n v/t/n v/t/n

//          String[] vtn = token[j].split("/");
//          try{  f.IDX_V[id] += Integer.parseInt(vtn[0]);  } catch (Exception e){} 
//          try{  f.IDX_T[id] += Integer.parseInt(vtn[1]);  } catch (Exception e){} 
//          try{  f.IDX_N[id] += Integer.parseInt(vtn[2]);  } catch (Exception e){}
          
//          String[] vtn = token[j].split("/");
//          if( vtn.length >= 1 && !vtn[0].isEmpty() )   f.IDX_V[id] += Integer.parseInt(vtn[0]);  else   continue;
//          if( vtn.length >= 2 && !vtn[1].isEmpty() )   f.IDX_T[id] += Integer.parseInt(vtn[1]);  else   continue;
//          if( vtn.length >= 3 && !vtn[2].isEmpty() )   f.IDX_N[id] += Integer.parseInt(vtn[2]);  else   continue;

          
          // fastest way to extract available index data
          String vtn = token[j];
          int pos_v = vtn.indexOf('/', 0);
          if( pos_v == -1 ){
            f.IDX_V[id] += Integer.parseInt(vtn);
            continue;
          } else 
            f.IDX_V[id] += Integer.parseInt(vtn.substring(0, pos_v));
          

          int pos_t = vtn.indexOf('/', ++pos_v);
          if( pos_t == -1 ){
            f.IDX_T[id] += Integer.parseInt(vtn.substring(pos_v));
            continue;
          } else if( pos_v != pos_t)
              f.IDX_T[id] += Integer.parseInt(vtn.substring(pos_v, pos_t));
  
          f.IDX_N[id] += Integer.parseInt(vtn.substring(++pos_t));
         
        }
        
        if( m_curr == null ){ //TODO
          if( g_curr == null ){
            g_curr = new IFS.Group(ifs, "default_group");
            tmp_g.put(g_curr.name, g_curr );
          }
          m_curr = new IFS.Mesh(g_curr, "default_mesh");
          g_curr.m.add(m_curr);
          tmp_m.put(m_curr.name, m_curr );
        }
        
        tmp_f.add( f );
        m_curr.f.add(f);
        f.parent = m_curr;
      }
      
    }
    
    
    float[][]  v  = tmp_v.toArray( new float[tmp_v.size()][3] ); 
    float[][]  t  = tmp_t.toArray( new float[tmp_t.size()][2] ); 
    float[][]  n  = tmp_n.toArray( new float[tmp_n.size()][3] ); 
    
    IFS.Group[] g = tmp_g.values().toArray( new IFS.Group[tmp_g.size()] ); 
    IFS.Mesh [] m = tmp_m.values().toArray( new IFS.Mesh [tmp_m.size()] );
    IFS.Face [] f = tmp_f.toArray( new IFS.Face [tmp_f.size()] ); 
      
    if( DEBUG ){
      System.out.println("(ObjFile)-------------");
      System.out.println("    v = "+ v.length);
      System.out.println("    t = "+ t.length);
      System.out.println("    n = "+ n.length);
      System.out.println("    g = "+ g.length);
      System.out.println("    m = "+ m.length);
      System.out.println("    f = "+ f.length);
      System.out.println("----------------------");
    }
    
    return ifs.set(g, m, f, v, t, n);
  }
  
  
  
  
  
  private static final SimpleDateFormat sdf = new SimpleDateFormat("dd.MM.yyyy - HH:mm:ss");
  private static final String NL = System.getProperty("line.separator");
  
  private static int NUM_LINES = 0;
  private static String NL(){
    NUM_LINES++;
    return NL;
  }
  

  
  
  //TODO: boolean for overwriting file, boolean for adding comments
  /**
   * creates an obj-file of the given IFS. Any existing file will be overwritten.<br>
   * Faces may be triangles, quads or polygons with any number of vertices.<br>
   * <br>
   * Use the mask VTNF to indicate which data gets exported.<br>
   * 0x1111 vertices, texcoords, normals and faces are exported.<br>
   * 0x1000 vertices only --> point-cloud.<br>
   * 0x1001 vertices and faces.
   * 
   * 
   * @param filename full path, e.g. System.getProperty("user.dir")+"/data/obj_export/test.obj"
   * @param ifs      IFS to export
   * @param VTNF     mask
   * 
   * @return true on success.
   */
  public static boolean create(String filename, IFS ifs, int VTNF, boolean replace_file){ 
    
    
    
    
    File file = new File(filename);
    
    // check file
    {
      String name = file.getName();
      String obj_name = "";
      String obj_ext  = "";
      int i = name.lastIndexOf('.');
      if (i > 0) {
        obj_name = name.substring(0,i);
        obj_ext  = name.substring(i+1);
      }
      if( !obj_ext.equals("obj")){
        System.err.println("(ObjFile.create) file has wrong extension \""+obj_ext+"\", not \"obj\"");
        return false;
      }
      
//      System.out.println("obj_name = "+obj_name);
//      System.out.println("obj_ext  = "+obj_ext);
  
      if( !replace_file && file.exists()){
  
        String path = file.getParentFile().getAbsolutePath();
        int file_nr = 0;
        do{
          file = new File(path, String.format("%s_%03d.%s", obj_name, file_nr++, obj_ext ));
        } while( file.exists() );
        
      }
    }
//    System.out.println(file.getName());
    
    
    
    boolean CREATED = true;
    BufferedWriter bw = null;

    try{
      bw = new BufferedWriter(new FileWriter(file));
      
      NUM_LINES = 1;
      
      // prepare/repair mask, necessary?
      int VTNF_tmp = VTNF;
      VTNF &= 0xFFFF0000; // clear bits
      if( (VTNF_tmp & 0xF000) != 0 ) VTNF |= 0x1000;
      if( (VTNF_tmp & 0x0F00) != 0 ) VTNF |= 0x0100;
      if( (VTNF_tmp & 0x00F0) != 0 ) VTNF |= 0x0010;
      if( (VTNF_tmp & 0x000F) != 0 ) VTNF |= 0x0001;
      
      VTNF &= (ifs.v != null) & (ifs.f[0].IDX_V!=null) & (ifs.f[0].IDX_V[0]!=-1) ? 0xFFFF : 0x0FFF;
      VTNF &= (ifs.t != null) & (ifs.f[0].IDX_T!=null) & (ifs.f[0].IDX_T[0]!=-1) ? 0xFFFF : 0xF0FF;
      VTNF &= (ifs.n != null) & (ifs.f[0].IDX_N!=null) & (ifs.f[0].IDX_N[0]!=-1) ? 0xFFFF : 0xFF0F;
      
      // file comments
      bw.write("#"+NL());
      bw.write("# obj export: Thomas Diewald - 2013."+NL());
      bw.write("# date: "+sdf.format(Calendar.getInstance().getTime())+NL());
      bw.write("# "+file.getName() +NL());
      bw.write("#"+NL());
      bw.write("# groups     = "+((ifs.g==null)?0:ifs.g.length) +NL());
      bw.write("# objects    = "+((ifs.m==null)?0:ifs.m.length) +NL());
      bw.write("# faces      = "+(((VTNF & 0x000F)!=0)?ifs.f.length:0) +NL());
      bw.write("# vertices   = "+(((VTNF & 0xF000)!=0)?ifs.v.length:0) +NL());
      bw.write("# tex-coords = "+(((VTNF & 0x0F00)!=0)?ifs.t.length:0) +NL());   
      bw.write("# normals    = "+(((VTNF & 0x00F0)!=0)?ifs.n.length:0) +NL());  
      bw.write("#"+NL());

      
//      final int VTN = (V?0x100:0) | (T?0x010:0) | (N?0x001:0);
      
      if( (VTNF & 0xFFF0) != 0 ) // v, t, n
      {
        float[][] v = ifs.v, t = ifs.t, n = ifs.n;

        if( (VTNF & 0xF000) != 0 ){ // vertices
          bw.write(NL()+"# vertices L:"+(NUM_LINES+1)+"-"+(NUM_LINES+1+v.length)+NL()); 
          for( int i=0; i<v.length; i++) bw.write("v " +v[i][0]+" "+v[i][1]+" "+v[i][2]+NL());
        }
        if( (VTNF & 0x0F00) != 0 ){ // tex-coords
          bw.write(NL()+"# texcoords L:"+(NUM_LINES+1)+"-"+(NUM_LINES+1+t.length)+NL()); 
          for( int i=0; i<t.length; i++) bw.write("vt "+t[i][0]+" "+t[i][1]+NL());
        }
        if( (VTNF & 0x00F0) != 0 ){ // normals
          bw.write(NL()+"# normals L:"+(NUM_LINES+1)+"-"+(NUM_LINES+1+n.length)+NL()); 
          for( int i=0; i<n.length; i++)bw.write("vn "+n[i][0]+" "+n[i][1]+" "+n[i][2]+NL());
        }
      }
      
      if( (VTNF & 0x000F) != 0 ) // faces
      {
        bw.write(NL());
        IFS.Group[] ifs_groups = ifs.g;
        int num_groups = ifs_groups.length;
        for(int g_id = 0; g_id < num_groups; g_id++){
          IFS.Group ifs_group = ifs_groups[g_id];
          
          // Group
          bw.write("g "+ifs_group.name+NL()); 
          
          ArrayList<IFS.Mesh> ifs_meshes = ifs_group.m;
          int num_meshes = ifs_meshes.size();
          for(int m_id = 0; m_id < num_meshes; m_id++){
            IFS.Mesh ifs_mesh = ifs_meshes.get(m_id);

            // Object/Mesh
            bw.write("o "+ifs_mesh.name+NL()); 
            
            ArrayList<IFS.Face> ifs_faces = ifs_mesh.f;
            int num_faces = ifs_faces.size();
            int num_triangles = 0, num_quads = 0, num_polygons = 0;
            for( int f_id = 0; f_id < num_faces; f_id++){
              IFS.Face ifs_face = ifs_faces.get(f_id);

              int[] v = ifs_face.IDX_V; 
              int[] t = ifs_face.IDX_T; 
              int[] n = ifs_face.IDX_N; 

              // Face
              bw.write("f");
              int i = 0, j = v.length; 
              if( j == 3 ) num_triangles++;
              if( j == 4 ) num_quads++;
              if( j >= 5 ) num_polygons++;
              switch(VTNF){
                case 0x1001: for(;i<j;i++) bw.write(" "+ (v[i]+1)                               ); break; // f v v v ...            
                case 0x1101: for(;i<j;i++) bw.write(" "+ (v[i]+1) +"/"+ (t[i]+1)                ); break; // f v/t v/t v/t ...      
                case 0x1011: for(;i<j;i++) bw.write(" "+ (v[i]+1) +"/"+           "/"+ (n[i]+1) ); break; // f v//n v//n v//n ...   
                case 0x1111: for(;i<j;i++) bw.write(" "+ (v[i]+1) +"/"+ (t[i]+1) +"/"+ (n[i]+1) ); break; // f v/t/n v/t/n v/t/n ...
                default: System.out.println("(ObjFile.create) no face-vertices?"); 
              }
              
              bw.write(NL());
            }
            bw.write("# "
              +num_faces    +" faces ("
              +num_triangles+" triangles, "
              +num_quads    +" quads, "
              +num_polygons +" polygons)"+NL()); 
            bw.write(NL());
          }
          bw.write(NL());
        }
      }// faces
      bw.write("# EOF L:"+NUM_LINES);

      System.out.println("created obj-file: "+file);

    }catch ( IOException e){
      e.printStackTrace();
      CREATED = false;
    }finally{
      try{
        if ( bw != null ){
          bw.flush();
          bw.close();
        }
      }catch ( final IOException e){
        e.printStackTrace();
      }
    }
    
    return CREATED;
  }



}

