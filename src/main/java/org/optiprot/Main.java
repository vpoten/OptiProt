/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot;

import org.jumpmind.symmetric.csv.CsvReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.optiprot.io.BioJavaStructureReader;
import org.optiprot.io.XyzStructureWriter;
import org.optiprot.maths.CalcAlignment;
import org.optiprot.metaheuristic.de.DifferentialEvolution;
import org.optiprot.metaheuristic.de.IDESolutionEval;
import org.optiprot.metaheuristic.de.impl.BinCrossover;
import org.optiprot.metaheuristic.de.impl.RandMutation;
import org.optiprot.metaheuristic.fans.FANSDouble;
import org.optiprot.metaheuristic.fans.impl.FANSNeighborBase;
import org.optiprot.neural.ElmanNetwork;
import org.optiprot.neural.JordanNetwork;
import org.optiprot.neural.NNPattern;
import org.optiprot.neural.NeuralNetwork;
import org.optiprot.neural.ProdUnitNetwork;
import org.optiprot.pockets.CalcPockets;
import org.optiprot.pockets.GridProbeCluster;
import org.optiprot.secpredict.SecStructDeEvaluator;
import org.optiprot.secpredict.SecStructFANSOper;
import org.optiprot.secpredict.SecStructPattGen1;
import org.optiprot.secpredict.SecStructPattManager;

/**
 *
 * @author victor
 */
public class Main {
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        /**
         * Example calls:
         *
         * pocket /home/victor/wwwpdb /home/victor/complexes.csv
         * fans /home/victor/work_bio/optiprot /home/victor/work_bio/romo1.fasta
         * astex /home/victor/software/ccdc_astex /home/victor/astex_pdbs.csv
         * setrain /home/victor/wwwpdb /home/victor/recentpdb_codes25.txt /home/victor/work_bio/optiprot
         */
     
        if( args.length<3 ){
            System.out.println("OptiProt FANS <workdir> <sequence_file>");
            System.out.println("OptiProt POCKET <pdbdir> <complexes (csv)>");
            System.out.println("OptiProt ASTEX <astexdir> <pdbcodes (csv)>");
            System.out.println("OptiProt SETRAIN <pdbdir> <pdbcodes (txt)> <workdir> [<init_pop>]");
            System.exit(0);
        }

        if( !args[0].toUpperCase().equals("FANS") &&
                !args[0].toUpperCase().equals("POCKET") &&
                !args[0].toUpperCase().equals("ASTEX") &&
                !args[0].toUpperCase().equals("SETRAIN") ){
            System.err.println( args[0]+" operation, not valid.");
            System.err.println( "Valid operations : FANS, POCKET, ASTEX, SETRAIN");
            System.exit(-1);
        }


        if( args[0].toUpperCase().equals("FANS") ){
            operationFANS( args );
        }
        else if( args[0].toUpperCase().equals("POCKET") ){
            operationPocket( args );
        }
        else if( args[0].toUpperCase().equals("ASTEX") ){
            OptiProtDockRun.operationAstex( args );
        }
        else if( args[0].toUpperCase().equals("SETRAIN") ){
            operationSETrain( args );
        }

        
    }



    /**
     *
     * 
     * @param args
     */
    private static void operationFANS( String[] args ){

        OptiProtParameters parameters=null;

        try {
            parameters = OptiProtParameters.createParameters(args[1], "/home/victor/work_bio/rotPDBs");
        } catch (Exception ex) {
            System.err.println( "Cannot create parameters (force field/rotamer library).");
            System.exit(-1);
        }

        System.out.println("Loaded force field parameters...");
        System.out.println("Loaded rotamer library...");



        String sequence="";

        try{
            sequence=BioJavaStructureReader.readSequence( args[2] );
        }
        catch ( Exception ex) {
          //not in fasta format or wrong alphabet
           System.err.println( args[2]+" bad file, not in fasta format or wrong alphabet." );
           System.exit(-1);
        }


        System.out.println("Readed sequence...");

        Calendar cal=Calendar.getInstance();
        String outfile=String.format("optiprot.%1$tY%1$tm%1$td_%1$tH%1$tM%1$tS", cal);

        System.out.println("Starting optimization...");

        OptiProtResults results=null;

        parameters.setCalcMMechanics(false);
        parameters.setCalcGBorn(false);

        try {
            results=OptiProtFANSRun.run(sequence, null, parameters, 300, 0);
        } catch (StructureException ex) {
            System.err.println( "Error in FANS.");
            System.exit(-1);
        }

        System.out.println( "\nIntermediate results 1 (SASA):" );
        printResults(results);

        parameters.setCalcGBorn(true);
        parameters.setCalcSASA(false);
        parameters.setMutAngleStepMax(6);

        try {
            results =
                    OptiProtFANSRun.run(sequence, results.getBestSolution(), parameters, 600, 0);
        } catch (StructureException ex) {
            System.err.println( "Error in FANS.");
            System.exit(-1);
        }

        System.out.println( "\nIntermediate results 2 (GBorn):" );
        printResults(results);

        parameters.setCalcMMechanics(true);
        parameters.setMutAngleStepMax(2);

        try {
            results =
                    OptiProtFANSRun.run(sequence, results.getBestSolution(), parameters, 1200, 0);
        } catch (StructureException ex) {
            System.err.println( "Error in FANS.");
            System.exit(-1);
        }

        System.out.println( "\nEnd search." );
        System.out.println( "------------------" );

        printResults(results);


        System.out.println("\nWriting file : "+outfile );
        try {
            XyzStructureWriter.writeStructure( results.getBestSolution(), outfile, parameters.getWorkDir());
        } catch (Exception ex) {
            System.err.println( "Cannot write output structure.");
            System.exit(-1);
        }
    }


    /**
     * 
     * @param results
     */
    private static void printResults( OptiProtResults results ){

        System.out.println( "Num. iterations = "+results.getIterations() );
        System.out.println( "Init. fitness = "+results.getInitFitness() );
        System.out.println( "Best fitness = "+results.getBestFitness() );
        System.out.println( "Start time = "+String.format( "%1$tY-%1$tm-%1$td, %1$tH:%1$tM:%1$tS", results.getStartTime()));
        System.out.println( "End time = "+String.format( "%1$tY-%1$tm-%1$td, %1$tH:%1$tM:%1$tS", results.getEndTime()));
    }
    


    /**
     *
     *
     * @param args
     */
    private static void operationPocket(String[] args) {

        OptiProtParameters parameters=new OptiProtParameters();

        parameters.setPDBDir( args[1] );
        
        CsvReader reader=null;

        try {
            reader = new CsvReader(args[2]);
        } catch (FileNotFoundException ex) {
           System.err.println( args[2]+" file not exists." );
           System.exit(-1);
        }

        int npockets = 3;
        double distLim=4.0;
        int total=0;
        int numHitsComplex=0;
        int numHitsApo=0;
        int numHitsComplexTop=0;
        int numHitsApoTop=0;

        try {
            reader.readHeaders();
            
            while (reader.readRecord() ){

                total++;

                String pdbcomplex=reader.get(0);
                String pdbapo=reader.get(1);
                String protein=reader.get(2);
                String ligand=reader.get(3).trim();
                String otherligand=reader.get(4).trim();
                String remcomplex=reader.get(5).trim();
                String remapo=reader.get(6).trim();
                

                Chain chainComplex = null;
                Chain chainApo = null;
                Chain chligand = null;

                //get the chains and ligands
                try{
                    Structure strcomplex = BioJavaStructureReader.readStructurePDB(
                            parameters.getPDBDir(), pdbcomplex);

                    if( ligand.length()==3 ){
                        // ligand = a single heteroatom

                        chainComplex =
                                BioJavaStructureReader.getAllChains(strcomplex, remcomplex);

                        chligand=new ChainImpl();
                        chligand.addGroup( 
                                BioJavaStructureReader.getHetatom( 
                                strcomplex, ligand, remcomplex)
                                );
                    }
                    else if( ligand.length()==1 ){
                        // ligand = a polipeptide

                        ArrayList<Integer> list= new  ArrayList<Integer>();
                        int i=0;

                        for( Chain ch : strcomplex.getChains() ){
                            if( ch.getName().equals(ligand) ){
                                chligand=ch;
                            }
                            else{
                                list.add(i);
                            }

                            i++;
                        }

                        int [] array=new int [list.size()];

                        for(i=0;i<array.length;i++){
                            array[i]=list.get(i);
                        }

                        chainComplex = 
                                BioJavaStructureReader.getChains(strcomplex,
                                array );
                    }
                    else{
                        // ligand = multiple molecules [XXX-YYY-...]

                        chligand=new ChainImpl();

                        for(int i=1;i<ligand.length();i+=4){

                            chligand.addGroup(
                                    BioJavaStructureReader.getHetatom(
                                    strcomplex,
                                    ligand.substring(i,i+3), remcomplex)
                                    );
                        }

                        chainComplex =
                                BioJavaStructureReader.getAllChains(strcomplex, remcomplex);
                    }

                    
                    Structure strapo = BioJavaStructureReader.readStructurePDB(
                            parameters.getPDBDir(), pdbapo);

                   
                    chainApo=BioJavaStructureReader.getAllChains(strapo, remapo);

                }
                catch( IOException e1 ){
                    System.out.println( "Error reading "+pdbcomplex+", "+pdbapo );
                    continue;
                }
                catch( StructureException e2 ){
                    System.out.println( ligand+" not found" );
                    continue;
                }


                GridProbeCluster[] result=null;
                int hit=-1;

                try {
                    // calc pockets for bound form
                    result =
                            CalcPockets.calc( chainComplex, parameters.getGrid(), npockets);

                    hit=checkHit( result, chligand, distLim);

                } catch (StructureException ex) {

                    System.out.println( "Error while calculating pockets of "+pdbcomplex );
                }
                finally{

                    if( hit>=0 ){
                        System.out.println( "Hit for "+pdbcomplex+" (complex) on pocket "+(hit+1) );

                        if( hit==0 )
                            numHitsComplex++;

                        numHitsComplexTop++;
                    }
                    else{
                        System.out.println( "No hit for "+pdbcomplex+" (complex)" );
                    }
                }

                try {
                    // calc pockets for apo form
                    hit=-1;
                    
                    //transform chainApo to meet chainComplex (structure alignment)
                    CalcAlignment.align(chainComplex, chainApo);

                    result =
                            CalcPockets.calc( chainApo, parameters.getGrid(), npockets);

                    hit=checkHit( result, chligand, distLim);

                    
                } catch (StructureException ex) {

                    System.out.println( "Error while calculating pockets of "+pdbapo );
                }
                finally{

                    if( hit>=0 ){
                        System.out.println( "Hit for "+pdbapo+" (apo) on pocket "+(hit+1) );

                        if( hit==0 )
                            numHitsApo++;

                        numHitsApoTop++;
                    }
                    else{
                        System.out.println( "No hit for "+pdbapo+" (apo)" );
                    }
                }

            }//end loop
        }
        catch( IOException ex){
            System.err.println( "error reading csv file." );
            System.exit(-1);
        }
        
        
        reader.close();

        System.out.println( "-----------------------\n" );
        System.out.println( total+" groups analyzed." );

        System.out.println( numHitsComplex + " hits for complex forms: " +
                100.0*(numHitsComplex/(double)total)+"%" );

        System.out.println( numHitsApo + " hits for apo forms: " +
                100.0*(numHitsApo/(double)total)+"%" );

        System.out.println( numHitsComplexTop + " top n hits for complex forms: "+
                100.0*(numHitsComplexTop/(double)total)+"%" );

        System.out.println( numHitsApoTop + " top n hits for apo forms: "+
                100.0*(numHitsApoTop/(double)total)+"%" );
    }

    
    private static double getMinDistance( GridProbeCluster cluster, Chain ligand )
            throws StructureException{

        Atom center=cluster.getCentroid();

        double minDist=java.lang.Double.MAX_VALUE;

        for( Group gr : ligand.getAtomGroups() ){
            for( Atom at : gr.getAtoms() ){
                double dist=Calc.getDistance(at, center);

                if( dist<minDist ){
                    minDist=dist;
                }
            }
        }

        return minDist;
    }

    
    private static int checkHit( GridProbeCluster [] clusters, Chain ligand,
            double distLimit ) throws StructureException{

        if( ligand==null || clusters==null )
            return -1;

        for( int i=0;i<clusters.length; i++ ){
            if( getMinDistance( clusters[i],ligand)<=distLimit ){
                return i;
            }
        }

        return -1;
    }


    /**
     * secondary structure NN training
     * 
     * @param args
     */
    private static void operationSETrain( String[] args ){

        int windowSize=11;
        int resPredPos=5;
        
        int nhid=125;//number of neurons in hidden layer

        double CR=0.9;// DE CR
        double F=0.8;//DE F
        int NP = 100;// DE NP (number of chromosomes)

        int nworkers=8;//number of threads
        double min_mse=0.05;
        int writeCicle=15;

        String workdir=args[3];


        SecStructPattManager patman=null;

        try {
            patman = new SecStructPattManager(args[2], args[1]);
        } catch (Exception ex) {
            System.err.println( "Error reading files: "+args[2]);
            System.exit(-1);
        }


        SecStructPattGen1 patgen=new SecStructPattGen1();
        patgen.setMode( SecStructPattGen1.MODE_COD3 );

        patman.genPatterns(windowSize, resPredPos, patgen );

        NNPattern pattern=patman.getPatGen().createPattern(windowSize);

        int nin=pattern.getNInputs();
        int nout=pattern.getNOutputs();
        
        NeuralNetwork network = new ElmanNetwork(nin, nout, nhid);
        ///NeuralNetwork network = new NeuralNetwork(nin, nout, nhid);
        ///NeuralNetwork network = new JordanNetwork(nin, nout, nhid);
        ///NeuralNetwork network = new ProdUnitNetwork(nin, nout, nhid);


        //prepare DE
        int dim = network.getNInputs()*network.getNHidden() +
                network.getNOutputs()*network.getNHidden();
        

        SecStructDeEvaluator evaluator= new SecStructDeEvaluator(network, patman);

        DifferentialEvolution devol=null;

        if( args.length<5 ){
            devol = new DifferentialEvolution();
            devol.createRandomPopulation( dim, evaluator, -1, 1 );
        }
        else{
            try {
                BufferedReader r = new BufferedReader(
                    new FileReader(workdir+args[4]) );

                devol = new DifferentialEvolution(r, evaluator);
            } catch (IOException ex) {
                System.err.println( "Error reading file: "+args[4]);
                System.exit(-1);
            }
        }

        devol.setCR( CR );
        devol.setF( F );
        devol.setNP( NP );
        devol.setCrossover( new BinCrossover() );
        devol.setMutation( new RandMutation() );

        //devol.setModifyF(true);
        devol.setHybridFCM(true);

        //add FANS local search to DE
        //devol.setFansFCM(true, null, 0.05);
        //devol.setFcmPeriod(25);
        //devol.setFansIter(50);
        ///

        ArrayList<IDESolutionEval> listEvals=new ArrayList<IDESolutionEval>();

        
        listEvals.add( evaluator );

        for(int i=1;i<nworkers;i++){
            listEvals.add( new SecStructDeEvaluator( new ElmanNetwork(nin, nout, nhid), patman.clone()) );
            ///listEvals.add( new SecStructDeEvaluator( new NeuralNetwork(nin, nout, nhid), patman.clone()) );
            ///listEvals.add( new SecStructDeEvaluator( new JordanNetwork(nin, nout, nhid), patman.clone()) );
            ///listEvals.add( new SecStructDeEvaluator( new ProdUnitNetwork(nin, nout, nhid), patman.clone()) );
        }

        devol.setWorkers( listEvals, 350 );


        double fitness=java.lang.Double.MAX_VALUE;
        int cicles=0;

        Calendar cal=Calendar.getInstance();

        String filename=workdir+String.format("depop.%1$tY%1$tm%1$td_%1$tH%1$tM%1$tS", cal);

        while( fitness > min_mse  ){
            devol.setIterLimit( devol.getIterLimit()+10 );
            devol.doSearch();
            fitness=devol.getBestSolution().getFitness();
            fitness /= (double)patman.getNumResidues();
            cicles++;

            if( (cicles%writeCicle)==0 ){
                System.out.println( "Iteration = "+devol.getIterations()+", mse = "+fitness );

                try {
                    devol.toFile( filename + "_" + cicles + ".txt");
                } catch (Exception ex) {
                    Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        int neval=evaluator.getNevaluations();

        try {
            devol.toFile( filename + "_final.txt" );
        } catch (Exception ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }

        System.out.println( "----------------------------");
        System.out.println( "Num. evaluations = "+neval*nworkers );
        System.out.println( "Num. DE iterations = "+devol.getIterations() );
        System.out.println( "Final mse = "+fitness );
        
    }
   

    
}
