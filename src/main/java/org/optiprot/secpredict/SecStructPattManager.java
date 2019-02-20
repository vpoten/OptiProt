/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.optiprot.io.BioJavaStructureReader;
import org.optiprot.neural.NNPattern;
import org.optiprot.neural.NeuralNetwork;

/**
 * Class that read and manages patterns for secondary structures
 *
 * @author victor
 */
public class SecStructPattManager {

    private static String [] AminoAcids = {
                            "ala", "arg", "asn", "asp", "cys",
                            "gln", "glu", "gly", "his", "ile",
                            "leu", "lys", "met", "phe", "pro",
                            "ser", "thr", "trp", "tyr", "val",
                            "sec" };

    //table that associates the amino acid code with its index
    private static HashMap<String, Integer> tableAAindex =
            new HashMap<String, Integer>();

    static{
        for(int i=0; i<AminoAcids.length; i++)
            tableAAindex.put( AminoAcids[i], i);
    }

    //types of secondary structures
    public static final int TYPE_HELIX = 0;
    public static final int TYPE_STRAND = 1;
    public static final int TYPE_COIL = 2;


    private List<String> pdbCodes = new ArrayList<String>();
    private List<List<NNPattern>> listsPatterns = new ArrayList<List<NNPattern>>();

    private HashMap<String, List<SecStructSubPattern>>  tableChains = new
            HashMap<String, List<SecStructSubPattern>>();

    private int numResidues=0;

    private int windowSize=0;
    private int resPredPos=0;
    private ISecStructPatternGen patGen=null;



    /**
     *
     * @param pathtofile path to file containing pdb codes + chain letter (CCCCX)
     * @param pdbdir local pdb directory
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     */
    public SecStructPattManager( String pathtofile, String pdbdir )
            throws FileNotFoundException, IOException {

        List<String> codes=new ArrayList<String>();

        readPDBCodes( pathtofile, codes );

        readChains( pdbdir, codes );

        //fill the list of pdb codes
        Iterator<String> it=tableChains.keySet().iterator();

        while(it.hasNext()){
            pdbCodes.add( it.next() );
            listsPatterns.add(new ArrayList<NNPattern>());
        }
    }


    /**
     * private constructor, used for cloning operation
     * @param other
     */
    private SecStructPattManager(SecStructPattManager other) {
        
        this.setResPredPos(other.getResPredPos());
        this.setWindowSize(other.getWindowSize());

        if( other.getPatGen()!=null )
            this.setPatGen(other.getPatGen().clone());

        this.numResidues=other.getNumResidues();

        for( String code : other.pdbCodes ){
            this.pdbCodes.add(code);
            this.listsPatterns.add(new ArrayList<NNPattern>());

            List<SecStructSubPattern> list=other.tableChains.get(code);

            List<SecStructSubPattern> list2=new ArrayList<SecStructSubPattern>();

            for( SecStructSubPattern subpat : list )
                list2.add(subpat.clone());

            this.tableChains.put(code, list2);
        }

        this.genPatterns( this.getWindowSize(), this.getResPredPos(), this.getPatGen() );

    }


    @Override
    public SecStructPattManager clone(){
        return new SecStructPattManager(this);
    }

    /**
     * returns the number of readed chains
     * 
     * @return
     */
    int getNumChains() {
        return tableChains.size();
    }


    /**
     * Reads a list of PDB codes from a file (one pdb code per line)
     *
     * @param pathtofile
     * @param codes I/O list to fill with the pdbcodes
     *
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     */
    private static void readPDBCodes( String pathtofile, List<String> codes )
            throws FileNotFoundException, IOException {

        BufferedReader r = new BufferedReader( new FileReader(pathtofile) );

        String l_line=r.readLine();

        while( l_line!=null ){

           l_line=l_line.trim();

           if( !l_line.isEmpty() ){
                codes.add(l_line);
           }

           l_line=r.readLine();
        }

        r.close();
    }


    /**
     * 
     * @param pdbdir
     * @param listcodes
     */
    private void readChains( String pdbdir, List<String> listcodes ){

        for( String code : listcodes ){
            try {
                Structure struct =
                        BioJavaStructureReader.readStructurePDB(pdbdir, code.substring(0, 4));

                Chain chAux=null;

                for( Chain ch : struct.getChains() ){
                    if( ch.getName().toUpperCase().equals(code.substring(4)) ){
                        chAux = ch;
                        break;
                    }
                }

                if( chAux==null )
                    continue;

                List<SecStructSubPattern> listSubpatt=new ArrayList<SecStructSubPattern>();

                getSubPatternChain( chAux, listSubpatt );
                this.numResidues+=listSubpatt.size();//count the total of residues

                tableChains.put(code, listSubpatt);

            } catch (IOException ex) {
                continue;
            } 

        }

    }

    /**
     * Extracts subpatterns from chain
     *
     * @param chain
     * @param listSubpatt
     */
    private void getSubPatternChain( Chain chain, List<SecStructSubPattern> listSubpatt ){
        List<Group> chainGrps=chain.getAtomGroups();

        int lim=chainGrps.size();

        for( int i=0; i<lim; i++  ){

            Group group=chainGrps.get(i);

            if( !(group instanceof AminoAcid) )
                continue;

            AminoAcid aa=(AminoAcid) group;

            AminoAcid prev=null;
            AminoAcid next=null;

            if( i>0 ){
                group=chainGrps.get(i-1);
                if( group instanceof AminoAcid )
                    prev=(AminoAcid) group;
            }

            if( i<lim-1 ){
                group=chainGrps.get(i+1);
                if( group instanceof AminoAcid )
                    next=(AminoAcid) group;
            }


            double phi=-1;
            double psi=-1;

            try{
                if( prev!=null )
                    phi=Calc.torsionAngle(prev.getC(),aa.getN(),aa.getCA(),aa.getC());
                if( next!=null )
                    psi=Calc.torsionAngle(aa.getN(),aa.getCA(),aa.getC(),next.getN());
            } catch (StructureException ex) {
                phi=psi=-1;
            }

            Map<String,String> map=aa.getSecStruc();

            int type=TYPE_COIL;

            if( map.containsValue("HELIX") ){
                type = TYPE_HELIX;
            }
            else if( map.containsValue("STRAND") ){
                type = TYPE_STRAND;
            }

            int aaindex=tableAAindex.get( aa.getPDBName().toLowerCase() );

            listSubpatt.add( new SecStructSubPattern( aaindex, type, -1, phi, psi) );
        }
    }


    /**
     * returns the aminoacid 3-letter code
     *
     * @param index : AA index
     * @return
     */
    public static String getAAName(int index){
        return AminoAcids[index];
    }


    public void genPatterns( int windowSize, int resPredPos, ISecStructPatternGen patgen ){
        setWindowSize(windowSize);
        setResPredPos(resPredPos);
        setPatGen(patgen);

        for( int i=0; i<this.pdbCodes.size(); i++ ){

            Iterator<NNPattern> it=new SecStructPattIterator( tableChains.get(pdbCodes.get(i)),
                getWindowSize(), getResPredPos(), getPatGen(), false );

            List<NNPattern> listPat=this.listsPatterns.get(i);

            while( it.hasNext() ){
                listPat.add(it.next());
            }

        }
    }


    public Iterator<NNPattern> getIterator4Chain(int idx){
        return this.listsPatterns.get(idx).iterator();
    }


    /**
     * @return the windowSize
     */
    public int getWindowSize() {
        return windowSize;
    }

    /**
     * @param windowSize the windowSize to set
     */
    protected void setWindowSize(int windowSize) {
        this.windowSize = windowSize;
    }

    /**
     * @return the resPredPos
     */
    public int getResPredPos() {
        return resPredPos;
    }

    /**
     * @param resPredPos the resPredPos to set
     */
    protected void setResPredPos(int resPredPos) {
        this.resPredPos = resPredPos;
    }

    /**
     * @return the patGen
     */
    public ISecStructPatternGen getPatGen() {
        return patGen;
    }

    /**
     * @param patGen the patGen to set
     */
    protected void setPatGen(ISecStructPatternGen patGen) {
        this.patGen = patGen;
    }

    /**
     * @return the numResidues
     */
    public int getNumResidues() {
        return numResidues;
    }

    /**
     * Classifies the internal patterns using a neural network
     *
     * @param network a network with 3 outputs
     * @return : the contingency matrix (double [3][3])
     */
    public double[][] classify( NeuralNetwork network ){

        double [][]mat=new double[3][3];

        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                mat[i][j]=0;

        int lim=this.getNumChains();

        for( int i=0; i<lim; i++){

            Iterator<NNPattern> it=this.getIterator4Chain(i);
            network.init();

            while(it.hasNext()){
                NNPattern pat=it.next();
                network.evaluate( pat.getInputs() );

                int exp=pat.getMaxOutput();
                int out=network.getMaxOutput();

                mat[exp][out]++;
            }
        }

        return mat;
    }


    /**
     * Classifies the patterns using a neural network
     *
     * @param pathtofile file that contains the pdbcodes to classify
     * @param pdbdir
     * @param network
     * @return
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     */
    public double[][] classify( String pathtofile, String pdbdir, NeuralNetwork network )
            throws FileNotFoundException, IOException{

        List<String> listCodes=new ArrayList<String>();

        //contingency matrix
        double [][]mat=new double[3][3];

        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                mat[i][j]=0;

        readPDBCodes( pathtofile, listCodes );

        for( String code : listCodes ){

            try {
                Structure struct =
                        BioJavaStructureReader.readStructurePDB(pdbdir, code.substring(0, 4));

                Chain chAux=null;

                for( Chain ch : struct.getChains() ){
                    if( ch.getName().toUpperCase().equals(code.substring(4)) ){
                        chAux = ch;
                        break;
                    }
                }

                if( chAux==null )
                    continue;

                List<SecStructSubPattern> listSubpatt=new ArrayList<SecStructSubPattern>();

                getSubPatternChain( chAux, listSubpatt );

                Iterator<NNPattern> it=new SecStructPattIterator( listSubpatt,
                    getWindowSize(), getResPredPos(), getPatGen(), true );

                while( it.hasNext() ){
                    NNPattern pat=it.next();
                    network.evaluate( pat.getInputs() );

                    int exp=pat.getMaxOutput();
                    int out=network.getMaxOutput();

                    mat[exp][out]++;
                }


            } catch (IOException ex) {
                continue;
            }
        }
        
        return mat;
    }
    
}
