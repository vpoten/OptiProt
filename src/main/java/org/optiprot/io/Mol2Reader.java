/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.StringTokenizer;
import org.optiprot.io.mol2.Mol2Atom;
import org.optiprot.io.mol2.Mol2Bond;
import org.optiprot.io.mol2.Mol2Structure;
import org.optiprot.io.mol2.Mol2SubStructure;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;

/**
 * Tripos Mol2 file reader
 * 
 * @author victor
 */
public class Mol2Reader {

    //mol2 file fields
    final static private String MOLECULE="@<TRIPOS>MOLECULE";
    final static public String ATOM="@<TRIPOS>ATOM";
    final static private String BOND="@<TRIPOS>BOND";
    final static private String SUBSTRUCTURE="@<TRIPOS>SUBSTRUCTURE";
    final static private String COMMENT="#";
    final static private String DICTIONARY="@<TRIPOS>DICT";
    final static private String SET="@<TRIPOS>SET";

    //OpenBabel tools tokens
    private static String OB_NUMROTBOND="Number of rotatable bonds:";
    

    /**
     * parse a Tripos Mol2 file
     *
     * @param path to file
     * @return
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     */
    static public Mol2Structure parseMol2File( String path )
            throws FileNotFoundException, IOException {

        BufferedReader r = new BufferedReader( new FileReader(path) );
        return parseMol2File(r);
    }

    static private Mol2Structure parseMol2File( BufferedReader r )
            throws FileNotFoundException, IOException {

        String l_line=r.readLine();

        Mol2Structure struct=null;

        while( l_line!=null )
        {
            l_line=l_line.trim();

            if( l_line.isEmpty() ){
                //nothing to do
            }
            else if( l_line.startsWith(COMMENT) ){
               // a comment
            }
            else if( l_line.startsWith(MOLECULE) ){
                struct = readMolecule(r);
            }


            l_line=r.readLine();
        }

        r.close();

        return struct;
    }

    /**
     *
     * @param line : line to parse
     * @param delimit : delimiters
     * @param num : number of tokens to obtain
     * @param tokens : array of tokens I/O
     * @return number of readed tokens
     */
    static private int getTokens( String line, String delimit, int num, String [] tokens ){

        StringTokenizer l_tokenizer = new StringTokenizer(line, delimit);

        int i=0;

        while (l_tokenizer.hasMoreTokens())
        {
            tokens[i++] = l_tokenizer.nextToken();

            if( i==(num-1) )
                break;
        }

        return i;
    }

    /**
     * read the atoms associated with RTI Atom
     * @param r
     * @param natoms
     * @param list where to append readed atoms
     * @throws java.io.IOException
     */
    private static void readAtoms(BufferedReader r, int natoms, List<Mol2Atom> list)
            throws IOException {

        String [] tokens = new String [16];

        for(int i=0;i<natoms;i++){
            String l_line=r.readLine();
            int ntokens = getTokens( l_line, " \t", 10, tokens);

            Mol2Atom atom=new Mol2Atom();

            atom.setName( tokens[1] );

            atom.setX( Double.parseDouble(tokens[2]) );
            atom.setY( Double.parseDouble(tokens[3]) );
            atom.setZ( Double.parseDouble(tokens[4]) );

            atom.setType( tokens[5] );

            if( ntokens>6 )
                atom.setSubstr( Integer.parseInt(tokens[6]) );

            if( ntokens>8 )
                atom.setCharge( Float.parseFloat(tokens[8]) );


            list.add(atom);
            
        }
    }

    /**
     * read the bonds associated with RTI Bond
     * @param r
     * @param nbonds
     * @param struct
     * @throws java.io.IOException
     */
    private static void readBonds(BufferedReader r, int nbonds, Mol2Structure struct) 
            throws IOException {

        String [] tokens = new String [16];

        for(int i=0;i<nbonds;i++){
            String l_line=r.readLine();
            int ntokens = getTokens( l_line, " \t", 5, tokens);

            Mol2Bond bond=new Mol2Bond();

            bond.setIdOrigin( Integer.parseInt(tokens[1]) );
            bond.setIdTarget( Integer.parseInt(tokens[2]) );
            bond.setType( tokens[3] );

            struct.getBonds().add( bond );
        }
    }

    /**
     * read the substructures associated with RTI Substructure
     * @param r
     * @param nsubstruct
     * @param struct
     * @throws java.io.IOException
     */
    private static void readSubstruct(BufferedReader r, int nsubstruct, Mol2Structure struct) 
            throws IOException {

        String [] tokens = new String [16];

        for(int i=0;i<nsubstruct;i++){
            String l_line=r.readLine();
            int ntokens = getTokens( l_line, " \t", 10, tokens);

            Mol2SubStructure substr=new Mol2SubStructure();

            substr.setName( tokens[1] );
            substr.setRootAtom( Integer.parseInt(tokens[2]) );

            if( ntokens>3 )
                substr.setType( tokens[3] );

            if( ntokens>5 )
                substr.setChain( tokens[5] );

            if( ntokens>6 )
                substr.setSubType( tokens[6] );

            if( ntokens>7)
                substr.setIntBonds( Integer.parseInt(tokens[7]) );

            struct.getSubStructs().add(substr);
        }
    }

    /**
     * 
     * @param r
     * @return
     * @throws java.io.IOException
     */
    private static Mol2Structure readMolecule( BufferedReader r )
            throws IOException {

        String [] tokens = new String [16];
        String l_line="";
        int ntokens=0;

        Mol2Structure struct=new Mol2Structure();

        //read the molecule name
        struct.setName( r.readLine() );

        //read the number of atoms, bonds, ...
        l_line = r.readLine();
        ntokens = getTokens( l_line, " \t", 5, tokens);

        int natoms=Integer.parseInt(tokens[0]);

        int nbonds=0;
        int nsubstruct=0;
        int nfeat=0;
        int nsets=0;

        if( ntokens>1 )
            nbonds=Integer.parseInt(tokens[1]);

        if( ntokens>2 )
            nsubstruct=Integer.parseInt(tokens[2]);

        if( ntokens>3 )
            nfeat=Integer.parseInt(tokens[3]);

        if( ntokens>4 )
            nsets=Integer.parseInt(tokens[4]);

        struct.setMolType( r.readLine() );
        struct.setChargeType( r.readLine() );
        r.readLine();//status bits
        struct.setComment( r.readLine() );

        int nregister=0;

        l_line = r.readLine();

        while( nregister<3 && l_line!=null ){
            

            if( l_line.isEmpty() ){
                //nothing to do
            }
            else if( l_line.startsWith(COMMENT) ){
               // a comment
            }
            else if( l_line.startsWith(ATOM) ){
                readAtoms(r, natoms, struct.getAtoms());
                nregister++;
            }
            else if( l_line.startsWith(BOND) ){
                readBonds(r, nbonds, struct);
                nregister++;
            }
            else if( l_line.startsWith(SUBSTRUCTURE) ){
                readSubstruct(r, nsubstruct, struct);
                nregister++;
            }

            l_line = r.readLine();
        }

        struct.arrange();
        
        return struct;
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
    public static DockingLigand readAstexLigand( String astex_path, String pdbcode, boolean minimized ){

        String path=getAstexLigPath( astex_path, pdbcode, minimized,"");
        
        Mol2Structure struct=null;

        try {
            struct = Mol2Reader.parseMol2File(path);
        }
        catch( FileNotFoundException ex){
            
            path=getAstexLigPath( astex_path, pdbcode, minimized,"1");
            
            try {
                struct = Mol2Reader.parseMol2File(path);
            } catch (Exception ex1) {
                return null;
            } 
        }
        catch (Exception ex) {
            return null;
        }

        DockingLigand ligand=new DockingLigand(
                BioJavaStructureReader.getAllChains(struct, null, false) );
        ////ligand.readCDKMol2(path);
        ligand.calcObProperties(path);

        return ligand;

    }

    private static String getAstexLigPath( String astex_path, String pdbcode, 
            boolean minimized, String suffix) {

        String path=astex_path + File.separator + "test"+
                pdbcode.toLowerCase()+File.separator+"ligand_reference"+suffix+".mol2";

        if( minimized ){
            path=astex_path + File.separator + "test"+
                pdbcode.toLowerCase()+File.separator+"ligand_reference"+suffix+"_min.mol2";
        }

        return path;
    }

    /**
     *
     * @param pdbcode
     * @param numchains : index of chains to read, if nulls reads all
     * @param onlyAA if true gets only AminoAcids
     * @return
     */
    public static DockingProtein readAstexProtein( String astex_path,
            String pdbcode, int [] numchains, boolean onlyAA ){

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
                BioJavaStructureReader.getChains(struct, numchains, onlyAA) );
        }
        else{
            prot=new DockingProtein(
                BioJavaStructureReader.getAllChains(struct, null, onlyAA) );
        }

        return prot;

    }

    /**
     * generates a conformer using OpenBabel obconformer command tool
     *
     * @param filepath path to mol2 file
     * @param numtest num of conformers to test
     * @param optsteps optimization steps for the best conformer
     * @return
     */
    public static DockingLigand generObConformer( String filepath, int numtest, int optsteps ){

        Runtime run=Runtime.getRuntime();

        String command="obconformer "+numtest+" "+optsteps+" "+filepath;

        Mol2Structure struct=null;
        Double energy=null;
        int numRotableBonds=0;

        try {
            Process pro = run.exec(command);
            pro.waitFor();

            BufferedReader r = new BufferedReader( new InputStreamReader(pro.getInputStream()) );
            
            struct=parseMol2File( r );

            r = new BufferedReader( new InputStreamReader(pro.getErrorStream()) );

            //get the last line of stderr to obtain the energy of conformer
            String l_line=r.readLine();
            String last="";

            while( l_line!=null ){
                last=l_line;
                l_line=l_line.trim();

                if( l_line.startsWith(OB_NUMROTBOND.toUpperCase()) ){
                    numRotableBonds =
                        Integer.parseInt( l_line.substring(OB_NUMROTBOND.length()).trim() );
                }

                l_line=r.readLine();
            }

            r.close();

            String [] tokens = new String [3];
            int ntokens = getTokens( last, " \t", tokens.length, tokens);
            energy=Double.parseDouble( tokens[1] );

        } catch (IOException ex) {
            return null;
        } catch (InterruptedException ex) {
            return null;
        }
        
        DockingLigand ligand=new DockingLigand(
                BioJavaStructureReader.getAllChains(struct, null, false) );

        ligand.setInternalEnergy(energy);
        ligand.setNumRotableBonds(numRotableBonds);

        return ligand;
    }

    
}
