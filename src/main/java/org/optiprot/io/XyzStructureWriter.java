/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import org.biojava.bio.structure.*;

/**
 *
 * @author victor
 */
public class XyzStructureWriter {

    private static String tinkerPath="";
    private static String forceField="";

    public static final String program="pdbxyz";

    public XyzStructureWriter()
    {

    }

    /**
     * @return the tinkerPath
     */
    public static String getTinkerPath() {
        return tinkerPath;
    }

    /**
     * @param tinkerPath the tinkerPath to set
     */
    public static void setTinkerPath(String tinkerPath) {
        XyzStructureWriter.tinkerPath = tinkerPath;
    }

    /**
     * @return the forceField
     */
    public static String getForceField() {
        return forceField;
    }

    /**
     * @param forceField the forceField to set
     */
    public static void setForceField(String forceField) {
        XyzStructureWriter.forceField = forceField;
    }

    public static void writeStructurePDB( Chain chain , String filepdb ) throws Exception
    {
        Structure struc=new StructureImpl(chain);

        writeString2File( struc.toPDB() , filepdb );
    }

    /**
     *
     * @param str : string to write
     * @param filename : path to file
     * @throws java.lang.Exception
     */
    public static void writeString2File( String str, String filename ) throws Exception
    {
        
        BufferedWriter out = new BufferedWriter(new FileWriter(filename));

        out.write( str );

        out.close();
    }

    /**
     * write structure to disk (to PDB file)
     *
     * @param p_chain
     * @param name
     * @param path
     * @throws java.lang.Exception
     */
    public static void writeStructure( Chain p_chain, String name, String path)
            throws Exception{

        XyzStructureWriter.writeStructurePDB( p_chain,
                path+File.separator+name+".pdb");
    }

    


    public static void writeStructureXyz( String filepdb ) throws Exception
    {

         Runtime run=Runtime.getRuntime();

         String path=XyzStructureWriter.getTinkerPath()+File.separator;

         String command=path+"bin"+File.separator+XyzStructureWriter.program
                 +" "+filepdb+" "+path+"params"+File.separator+XyzStructureWriter.getForceField();

         Process pro=run.exec( command );
         pro.waitFor();

    }
}
