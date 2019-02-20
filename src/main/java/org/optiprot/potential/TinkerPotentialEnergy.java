/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 * Use tinker program (analize) to obtain the potential energy of a protein
 *
 * @author victor
 */
public class TinkerPotentialEnergy {

    private static String tinkerPath="";
    private static String forceField="";
    
    public static final String program="analyze";
    public static final String idLineEnergy="Total Potential Energy";


    /**
     * @return the tinkerPath
     */
    static public String getTinkerPath() {
        return tinkerPath;
    }

    /**
     * @param tinkerPath the tinkerPath to set
     */
    static public void setTinkerPath(String tinkerPath) {
        TinkerPotentialEnergy.tinkerPath = tinkerPath;
    }

    /**
     * @return the forceField
     */
    static public String getForceField() {
        return forceField;
    }

    /**
     * @param forceField the forceField to set
     */
    static public void setForceField(String forceField) {
        TinkerPotentialEnergy.forceField = forceField;
    }

    
    /**
     * launchs a external process (tinker - analyze) to calculate de potential
     * energy of the structure
     * @param xyzFile
     * @return
     */
    static public Double calcEnergy( String xyzFile )
    {
       
       Runtime run=Runtime.getRuntime();
       
       String path=TinkerPotentialEnergy.getTinkerPath()+File.separator;
       
       
       
       String command=path+"bin"+File.separator+TinkerPotentialEnergy.program+
               " "+xyzFile+" "+path+"params"+File.separator+TinkerPotentialEnergy.getForceField()+" E";

       try{
           Process pro=run.exec( command );
           pro.waitFor();

           BufferedReader r = new BufferedReader( new InputStreamReader(pro.getInputStream()) );

           String l_line=r.readLine();

           while( l_line!=null )
           {
               l_line=l_line.trim();

               if( l_line.startsWith(TinkerPotentialEnergy.idLineEnergy) )
               {
                    StringTokenizer l_tokenizer = new StringTokenizer(l_line, " \t");
                    String l_curtoken = "";

                    while (l_tokenizer.hasMoreTokens())
                    {
                        l_curtoken = l_tokenizer.nextToken();

                        try{
                            Double.parseDouble(l_curtoken);
                            return Double.valueOf( l_curtoken );
                        }
                        catch( NumberFormatException ex){

                        }
                    }
               }

               l_line=r.readLine();
           }
           
       }
       catch( Exception e )
       {
            return null;
       }

       return null;
    }

}
