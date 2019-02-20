package org.optiprot.configtests;


import java.io.FileNotFoundException;
import java.io.IOException;
import org.optiprot.jgap.OptiProtConfiguration;
import org.optiprot.*;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Test;
import static org.junit.Assert.*;
import org.biojava.bio.structure.*;
import org.jgap.InvalidConfigurationException;
import org.optiprot.io.BioJavaStructureReader;
import org.optiprot.io.CharmmTopParReader;
import org.optiprot.io.Mol2Reader;
import org.optiprot.io.SNNSReader;
import org.optiprot.io.XyzStructureWriter;
import org.optiprot.io.mol2.Mol2Structure;
import org.optiprot.neural.NNPattern;
import org.optiprot.potential.CharmmForceField;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;
import org.optiprot.potential.element.CharmmElement;
import org.optiprot.potential.element.CharmmResidue;
import org.optiprot.rotamer.RotamerLibrary;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author victor
 */
public class ConfigOptiProtTest {

    private static OptiProtParameters parameters=null;
    private static OptiProtConfiguration configuration=null;

    public static String astex_path="/home/victor/software/ccdc_astex";

    private static Chain chain1=null;
    private static Chain chain2=null;

    static{

        parameters=new OptiProtParameters();
        
        getParameters().setTinkerPath("/home/victor/software/tinker");
        getParameters().setTinkerForceField("charmm22.prm");
        getParameters().setWorkDir("/home/victor/work_bio/optiprot");
        getParameters().setPDBDir( "/home/victor/wwwpdb" );

        try {
            configuration = new OptiProtConfiguration(getParameters());
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(ConfigOptiProtTest.class.getName()).log(Level.SEVERE, null, ex);
        }

        Structure struct=null;
        try {
            struct = BioJavaStructureReader.readStructure(
                    getParameters().getWorkDir() + File.separator + "eif6_swiss.pdb");
        } catch (IOException ex) {
            Logger.getLogger(ConfigOptiProtTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        chain1=struct.getChain(0);

        //load topology
        List<CharmmResidue> topology=null;
        ArrayList<CharmmElement> elements=new ArrayList<CharmmElement>();
        try {
            topology = CharmmTopParReader.parseTopFile(
                    "org/optiprot/data/top_all22_prot.inp",
                    elements);
        } 
        catch (Exception ex) {
            topology=null;
        }

        getParameters().setTopology(topology);

        //load force field
        CharmmForceField ffield=new CharmmForceField(elements);

        try {
            CharmmTopParReader.parseParFile(
                    "org/optiprot/data/par_all22_prot.inp",
                    ffield, elements);
        }
        catch (Exception ex) {
            Logger.getLogger(ConfigOptiProtTest.class.getName()).log(Level.SEVERE, null, ex);
            ffield=null;
        }
        getParameters().setForceField(ffield);
    }

    
    public ConfigOptiProtTest() {
    }

    static public void setName(String name){
        getParameters().setName(name);
    }

    /**
     * @return the parameters
     */
    public static OptiProtParameters getParameters() {
        return parameters;
    }

    /**
     * @return the configuration
     */
    public static OptiProtConfiguration getConfiguration() {
        return configuration;
    }

    /**
     * @return the chain1
     */
    public static Chain getChain1() {
        return chain1;
    }

    /**
     * @return the chain2
     */
    public static Chain getChain2() {
        return chain2;
    }

    /**
     * @param aChain2 the chain2 to set
     */
    public static void setChain2(Chain aChain2) {
        chain2 = aChain2;
    }

    /**
     * read a chain froma pdb file in the workdir
     *
     * @param pdbfile
     * @param num
     * @return
     */
    public static Chain readChain(String pdbfile, int [] numchains){

        Structure struct;
        try {
            struct = BioJavaStructureReader.readStructure(getParameters().getWorkDir() + File.separator + pdbfile);
        } catch (IOException ex) {
            return null;
        }

        return getChains( struct, numchains, true );
        
    }
    
    public static Chain readChainPDB(String pdbname, int [] numchains){

        Structure struct;
        try {
            struct = BioJavaStructureReader.readStructure(
                    getParameters().getPDBDir() + File.separator +
                    pdbname.substring(1, 3).toLowerCase() + File.separator +
                    "pdb" + pdbname.toLowerCase() + ".ent.gz");
        } catch (IOException ex) {
           return null;
        }

        return getChains( struct, numchains, true );
        
    }

    /**
     *
     * @param struct
     * @param numchains
     * @param onlyAA : if true gets only AminoAcids
     * @return
     */
    private static Chain getChains( Structure struct, int [] numchains, boolean onlyAA ){

        Chain ch2=new ChainImpl();

        for( int i=0;i<numchains.length;i++){
            Chain ch=struct.getChain(numchains[i]);

            for( Group group : ch.getAtomGroups() ){
                if( onlyAA ){
                    if( group instanceof AminoAcid ){
                        ch2.addGroup(group);
                    }
                }
                else{
                    ch2.addGroup(group);
                }
            }
        }

        return ch2;
    }


    /**
     * read a Hetatom froma pdb file in the workdir
     *
     * @param pdbfile
     * @param hetname
     * @param nchain
     * @return
     */
    public static HetatomImpl readHetatom(String pdbfile, String hetname, int nchain){

        Structure struct;
        try {
            struct = BioJavaStructureReader.readStructure(getParameters().getWorkDir() + File.separator + pdbfile);
        } catch (IOException ex) {
            return null;
        }

        return getHetatom(struct.getChain(nchain), hetname);
    }

    public static HetatomImpl readHetatomPDB(String pdbname, String hetname, int nchain){

        Structure struct;
        try {
            struct = BioJavaStructureReader.readStructure(
                    getParameters().getPDBDir() +
                    File.separator + pdbname.substring(1, 3).toLowerCase() +
                    File.separator + "pdb" + pdbname.toLowerCase() + ".ent.gz");
        } catch (IOException ex) {
            return null;
        }

        return getHetatom(struct.getChain(nchain), hetname);
    }

    
    private static HetatomImpl getHetatom( Chain ch, String hetname){

        for( Group group : ch.getAtomGroups() ){

            if( group instanceof HetatomImpl ){
                if( group.getPDBName().equals(hetname) )
                    return (HetatomImpl) group;
            }
        }

        return null;
    }

    /**
     * loads the rotamer library
     */
    public static void loadRotamerLibrary(){

        //create and load rotamer library
        RotamerLibrary lib=new RotamerLibrary();
        try {
            lib.loadLibrary("/home/victor/work_bio/rotPDBs", "allH");
        } catch (Exception ex) {
            Logger.getLogger(ConfigOptiProtTest.class.getName()).log(Level.SEVERE, null, ex);
            lib=null;
        }

        getParameters().setRotLib(lib);
    }

    /**
     * write structure to disk (to PDB file)
     *
     * @param p_chain
     * @param name
     * @throws java.lang.Exception
     */
    public static void writeStructure( Chain p_chain, String name) throws Exception{

        XyzStructureWriter.writeStructurePDB( p_chain,
                parameters.getWorkDir()+File.separator+name+".pdb");
    }


    /**
     * read a ligand from the astex validation set
     *
     * @param pdbcode
     * @param minimized : true = read the minimized ligand,
     * false = read cristal structure ligand
     * 
     * @return
     */
    public static DockingLigand readAstexLigand( String pdbcode, boolean minimized ){

        String path=getAstexLigPath( pdbcode, minimized);

        Mol2Structure struct=null;

        try {
            struct = Mol2Reader.parseMol2File(path);
        } catch (Exception ex) {
            return null;
        }

        DockingLigand ligand=new DockingLigand( struct.getChain(0) );
        ////ligand.readCDKMol2(path);
        ligand.calcObProperties(path);

        return ligand;

    }

    public static String getAstexLigPath(String pdbcode, boolean minimized) {

        String path=astex_path + File.separator + "test"+
                pdbcode.toLowerCase()+File.separator+"ligand_reference.mol2";

        if( minimized ){
            path=astex_path + File.separator + "test"+
                pdbcode.toLowerCase()+File.separator+"ligand_reference_min.mol2";
        }

        return path;
    }



    /**
     *
     * @param pdbcode
     * @param numchains, if nulls reads all
     * @param onlyAA if true gets only AminoAcids
     * @return
     */
    public static DockingProtein readAstexProtein( String pdbcode, int [] numchains, boolean onlyAA ){

        String path=astex_path + File.separator + "test"+
                pdbcode.toLowerCase()+File.separator+"protein.mol2";

        Mol2Structure struct=null;

        try {
            struct = Mol2Reader.parseMol2File(path);
        } catch (Exception ex) {
            return null;
        }

        DockingProtein prot=null;

        if( numchains!=null ){
            prot=new DockingProtein(
                getChains(struct, numchains, onlyAA) );
        }
        else{
            prot=new DockingProtein(
                BioJavaStructureReader.getAllChains(struct, null, onlyAA) );
        }

        return prot;

    }
    
    /**
     * reads a SNNS pat file from workdir
     * 
     * @param filename
     * @return
     */
    public static ArrayList<NNPattern> readSNNSFile( String filename ){

        try {
            return SNNSReader.parseSNNSFile(
                    getParameters().getWorkDir() + File.separator + filename );
        } catch (FileNotFoundException ex) {
            return null;
        } catch (IOException ex) {
            return null;
        }
    }

    @Test
    public void testParameters() {

        assertTrue(parameters!=null);
        assertTrue(parameters.getForceField()!=null);
        assertTrue(parameters.getTopology()!=null);
    }
}
