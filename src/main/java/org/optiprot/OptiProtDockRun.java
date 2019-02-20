/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot;

import org.jumpmind.symmetric.csv.CsvReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.biojava.bio.structure.Atom;
import org.optiprot.io.Mol2Reader;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcRmsd;
import org.optiprot.metaheuristic.MetaheuristicBase;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.de.DifferentialEvolution;
import org.optiprot.metaheuristic.de.impl.BinCrossover;
import org.optiprot.metaheuristic.de.impl.DockScoreEval;
import org.optiprot.metaheuristic.de.impl.RandMutation;
import org.optiprot.metaheuristic.meme.Memetic;
import org.optiprot.metaheuristic.meme.impl.MemeOperatorDock;
import org.optiprot.metaheuristic.tabu.TabuSearch;
import org.optiprot.metaheuristic.tabu.impl.DockConfSolution;
import org.optiprot.metaheuristic.tabu.impl.DockTabuNeighbor;
import org.optiprot.potential.docking.Docked;
import org.optiprot.potential.docking.SimpleDockScore;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;

/**
 *
 * @author victor
 */
public class OptiProtDockRun implements Runnable {

    private static int NWORK_THREADS=2;
    private static int TIME_PER_THREAD=70;//seconds

    private static int MODE_DE=1;
    private static int MODE_MEME=2;
    private static int MODE_TABU=3;

    private static int mode=MODE_DE;

    //parameters for memetic
    private static int NP = 5;
    private static int NMUT = 3;
    private static int NITER = 40;

    //parameters for tabu
    private static int TTENURE=10;

    //parameters for DE
    private static int NPDE = 50;
    private static double DECR = 0.9;
    private static double DEF = 0.8;


    /////////////////


    private ArrayList<String> pdbcodes=null;
    private String astexDir=null;
    private AtomGrid grid=new AtomGrid(30,30,30);

    private DockScoreEval evaluator = new DockScoreEval( new SimpleDockScore() );
    private MemeOperatorDock operator= null;
    private Memetic memetic = null;

    private TabuSearch tabu = null;

    private DifferentialEvolution devol = null;

    
    private int hits=0;

    private static double RMSD_MAX = 2.0;
    private int timeLimit = 60;//seconds
    private static int FINAL_OPT_STEPS = 40;

    static final private int OPT_MEME = 1;
    static final private int OPT_TABU = 2;
    static final private int OPT_DE = 3;

    private int modeOpt = OPT_MEME;
    
    /**
     * constructor for memetic optimization
     *
     * @param nmut : mutants introduced per generation
     * @param np : size of population
     * @param numiter : num of initial iterations
     */
    public OptiProtDockRun(int nmut, int np, int numiter) {

        operator=new MemeOperatorDock();
        memetic = new Memetic();

        operator.setNMutants(nmut);
        operator.setInitIterations(numiter);
        memetic.setOperators(operator);
        memetic.setNP(np);
        this.setModeOpt(OPT_MEME);
    }

    /**
     *
     * @param tenure
     */
    public OptiProtDockRun( int tenure ) {

        tabu = new TabuSearch();
        tabu.setNeighManager( new DockTabuNeighbor() );
        tabu.setTenure(tenure);
        this.setModeOpt(OPT_TABU);
    }

    /**
     * 
     * @param NPDE : population size
     * @param DECR : cross rate
     * @param DEF : mutation factor
     */
    OptiProtDockRun(int NPDE, double DECR, double DEF) {
        devol = new DifferentialEvolution();
        devol.setCR(DECR);
        devol.setNP(NPDE);
        devol.setF(DEF);
        devol.setCrossover( new BinCrossover() );
        devol.setMutation( new RandMutation() );
        this.setModeOpt(OPT_DE);
    }


    
    public void run() {

        for( String pdbcode : getPdbcodes() ){
            try{
                doDocking( pdbcode );
            }
            catch(Exception e){
                System.out.println("Error while calculating "+pdbcode );
            }
        }
    }

    /**
     * @return the pdbcodes
     */
    public ArrayList<String> getPdbcodes() {
        return pdbcodes;
    }

    /**
     * @param pdbcodes the pdbcodes to set
     */
    public void setPdbcodes(ArrayList<String> pdbcodes) {
        this.pdbcodes = pdbcodes;
    }

    /**
     * @return the astexDir
     */
    public String getAstexDir() {
        return astexDir;
    }

    /**
     * @param astexDir the astexDir to set
     */
    public void setAstexDir(String astexDir) {
        this.astexDir = astexDir;
    }

    /**
     * @return the hits
     */
    public int getHits() {
        return hits;
    }

    
    public void doDocking(String pdbcode){
        
        DockingProtein protein=Mol2Reader.readAstexProtein( getAstexDir(),pdbcode, null, false);
        DockingLigand ligand=Mol2Reader.readAstexLigand( getAstexDir(),pdbcode, false );

        if( protein==null || ligand==null ){
            System.out.println("Error reading "+pdbcode );
            return;
        }


        Atom actSiteCoM=ligand.getCoM();

        //read reference ligand to compare
        DockingLigand ligand_ref=Mol2Reader.readAstexLigand( getAstexDir(),pdbcode, false );

        //translates ligand to origin
        ligand.shift( CalcGeom.product(actSiteCoM, -1));

        Docked docked=new Docked( protein, ligand, actSiteCoM, grid );

        if( !docked.isAbleAssignActivePoints() ){
            throw new RuntimeException("cannot optimize "+pdbcode);
        }

        //set the evaluator to docked
        evaluator.setDocked( docked );
        docked.setEvaluator( evaluator );

        MetaheuristicBase metah=null;

        if( this.getModeOpt()==OPT_MEME ){
            memetic.reset();
            operator.setDocked(docked);
            metah=memetic;
        }
        else if( this.getModeOpt()==OPT_TABU ){
            tabu.reset();
            DockConfSolution initSol = new DockConfSolution( docked,
                DockConfSolution.MODE_OPT_POSE, FINAL_OPT_STEPS );
            tabu.setSolution(initSol);
            metah=tabu;
        }
        else if( this.getModeOpt()==OPT_DE ){
            devol.reset();
            double boxLength=ligand.calcRadii()*0.5;

            double [] lowers=new double [] { -1, -1, -1, -1,
                actSiteCoM.getX()-boxLength, actSiteCoM.getY()-boxLength, actSiteCoM.getZ()-boxLength };

            double [] uppers=new double [] { 1, 1, 1, 1,
                actSiteCoM.getX()+boxLength, actSiteCoM.getY()+boxLength, actSiteCoM.getZ()+boxLength};
            
            devol.setBounds(lowers, uppers);
            devol.createRandomPopulation( DockScoreEval.NUM_PARAMETERS, evaluator);

            metah=devol;
        }

        
        metah.setTimeLimit(getTimeLimit());

        metah.doSearch();

        if( metah.getBestSolution() instanceof DockConfSolution ){

            DockConfSolution bestSol=(DockConfSolution) metah.getBestSolution();
            bestSol.optimize(FINAL_OPT_STEPS);

            docked.assignAtomToActivePoint( bestSol.getLigActiveAtom(), bestSol.getActivePoint());
            docked.setTransform( bestSol.getMTrans() );
        }
        else{
            DESolution bestsol=(DESolution) metah.getBestSolution();
            docked.setTransform( evaluator.getMTrans(bestsol) );
        }

        double rmsd=CalcRmsd.bruteRmsd(ligand_ref.getAtoms(), docked.getLigandAtoms());

        if( rmsd<RMSD_MAX ){
            System.out.println("Hit with "+pdbcode+", rmsd="+rmsd );
            hits++;
        }
        else{
            System.out.println("No hit with "+pdbcode+", rmsd="+rmsd );
        }
        
    }

    /**
     * @return the timeLimit
     */
    public int getTimeLimit() {
        return timeLimit;
    }

    /**
     * @param timeLimit the timeLimit to set
     */
    public void setTimeLimit(int timeLimit) {
        this.timeLimit = timeLimit;
    }

    /**
     * @return the modeOpt
     */
    public int getModeOpt() {
        return modeOpt;
    }

    /**
     * @param modeOpt the modeOpt to set
     */
    public void setModeOpt(int modeOpt) {
        this.modeOpt = modeOpt;
    }


    /**
     *
     * @param args
     */
    public static void operationAstex( String[] args ){

        String astex_dir=args[1];

        ArrayList<String> pdbcodes=new ArrayList<String>();

        //read pdbcodes
        CsvReader reader=null;

        try {
            reader = new CsvReader(args[2]);
        } catch (FileNotFoundException ex) {
           System.err.println( args[2]+" file not exists." );
           System.exit(-1);
        }

        try {

            while (reader.readRecord() ){
                for(int i=0;i<reader.getColumnCount();i++){
                    String code=reader.get(i).trim();

                    if( !code.isEmpty() )
                        pdbcodes.add(code);
                }
            }

            reader.close();
        }
        catch( IOException ex){
            System.err.println( "error reading csv file." );
            System.exit(-1);
        }


        Calendar cal=Calendar.getInstance();
        System.out.println("Start time "+String.format("%1$tY-%1$tm-%1$td, %1$tH:%1$tM:%1$tS", cal) );
        System.out.println("Threads="+NWORK_THREADS+", Time="+TIME_PER_THREAD);

        if( mode==MODE_MEME ){
            System.out.println("Memetic");
            System.out.println("NP="+NP+", NMUT="+NMUT+", NITER="+NITER);
        }
        else if( mode==MODE_TABU ){
            System.out.println("Tabu");
            System.out.println("TTENURE="+TTENURE );
        }
        else if( mode==MODE_DE ){
            System.out.println("DE");
            System.out.println("NP="+NPDE+", CR="+DECR+", F="+DEF );
        }


        // divide the work in threads
        ArrayList<OptiProtDockRun> listDocks=new ArrayList<OptiProtDockRun>();

        int nthreads=NWORK_THREADS;//num of working threads
        int size=pdbcodes.size()/nthreads;//docks per threads

        ExecutorService es =  Executors.newFixedThreadPool(nthreads+1);


        for( int i=0;i<nthreads;i++){
            OptiProtDockRun dockRun=null;

            if( mode==MODE_MEME )
                dockRun=new OptiProtDockRun( NMUT, NP, NITER );
            else if( mode==MODE_TABU )
                dockRun=new OptiProtDockRun( TTENURE );
            else if( mode==MODE_DE )
                dockRun=new OptiProtDockRun( NPDE, DECR, DEF );

            listDocks.add(dockRun);
            dockRun.setAstexDir(astex_dir);

            ArrayList<String> lcodes=new ArrayList<String>();

            for( int j=i*size; j<(i+1)*size; j++){
                lcodes.add( pdbcodes.get(j) );
            }

            if( i==(nthreads-1) ){
                for( int j=(i+1)*size; j<pdbcodes.size(); j++){
                    lcodes.add( pdbcodes.get(j) );
                }
            }

            dockRun.setPdbcodes( lcodes );
            dockRun.setTimeLimit(TIME_PER_THREAD);

            es.execute( dockRun );
        }


        es.shutdown();

        int sleeptime = 1000*TIME_PER_THREAD/3;//milliseconds

        //wait all threads to finish
        while( !es.isTerminated() ){

            try {
                Thread.sleep(sleeptime);
            } catch (InterruptedException ex) {
                System.err.println( "Interrupted exception in thread." );
            }
        }


        int hits=0;

        for( OptiProtDockRun dockRun : listDocks){
            hits+=dockRun.getHits();
        }



        //end
        System.out.println( "-----------------------\n" );
        cal=Calendar.getInstance();
        System.out.println("End time "+String.format("%1$tY-%1$tm-%1$td, %1$tH:%1$tM:%1$tS", cal) );

        System.out.println( pdbcodes.size()+" complexes analyzed." );

        System.out.println( hits + " hits: " +
                100.0*(hits/(double)pdbcodes.size())+"%" );


    }

}
