/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.LinkedList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * class that implements co-operative strategy for optimization
 *
 * @author victor
 */
public class MetaCooperative {

    protected static int MEMORY_LENGTH=50;

    static private int cores = 2;
    
    private ArrayList<IMetaheuristic> solvers=new ArrayList<IMetaheuristic>();
    private ArrayList<SolverInfo> solversInfo = new ArrayList<SolverInfo>();
    protected LinkedList<IMetaheuristicSolution> latestSol=new LinkedList<IMetaheuristicSolution>();
    protected LinkedList<Double> latestImpr = new LinkedList<Double>();

    //arraylist of sorted costs and improve rates
    private ArrayList<IMetaheuristicSolution> solMemory =
            new ArrayList<IMetaheuristicSolution>(MEMORY_LENGTH);

    private ArrayList<Double> improveMemory = new ArrayList<Double>(MEMORY_LENGTH);

    protected IMetaheuristicSolution bestSol=null;
   
    //constants for fuzzy rules
    ////
    static protected double BETA = 20.0;//percentile limit 1
    static protected double ALPHA = 10.0;//percentile limit 2
    static protected double LAMBDA = 0.7;//acceptability
    static protected double DELTA_LAMBDA = 0.05;

    

    

    /**
     * class that stores information about the solvers
     * 
     */
    protected class SolverInfo {

        //time at t ant t-1
        public Date time=null;
        public Date time_1=null;

        //solutions at t and t-1
        public IMetaheuristicSolution solution=null;
        public IMetaheuristicSolution solution_1=null;

        public SolverInfo() {
        }

        public void update( IMetaheuristicSolution p_solution, Date p_time){
            time_1=time;
            solution_1=solution;

            time=p_time;
            solution=p_solution;
        }

        public Double getImprove(){
            if( time==null || time_1==null )
                return null;

            double improve=solution.getFitness() - solution_1.getFitness();

            if( solution.isBetter(solution_1) ){
                improve=Math.abs(improve);
            }
            else{
                improve=-Math.abs(improve);
            }

            return improve / ( time.getTime() - time_1.getTime() );
        }
    }

    
    /////////////////////

    public MetaCooperative() {
        for( int i=0;i<MEMORY_LENGTH;i++ ){
            this.getImproveMemory().add( -java.lang.Double.MAX_VALUE );
            this.getSolMemory().add(null);
            this.getLatestSol().add(null);
            this.getLatestImpr().add( -java.lang.Double.MAX_VALUE );
        }
    }

    /**
     * @return the solMemory
     */
    protected ArrayList<IMetaheuristicSolution> getSolMemory() {
        return solMemory;
    }

    /**
     * @return the improveMemory
     */
    protected ArrayList<Double> getImproveMemory() {
        return improveMemory;
    }

    /**
     * @return the latestSol
     */
    protected LinkedList<IMetaheuristicSolution> getLatestSol() {
        return latestSol;
    }

    /**
     * @return the latestImpr
     */
    protected LinkedList<Double> getLatestImpr() {
        return latestImpr;
    }

    /**
     * @return the solvers
     */
    protected ArrayList<IMetaheuristic> getSolvers() {
        return solvers;
    }

    /**
     * @return the solversInfo
     */
    protected ArrayList<SolverInfo> getSolversInfo() {
        return solversInfo;
    }

    /**
     * @return the cores
     */
    public static int getCores() {
        return cores;
    }

    /**
     * @param cores the cores to set
     */
    public static void setCores(int ncores) {
        cores = ncores;
    }



    /**
     * adds a solver to the coordinator
     * 
     * @param solver
     */
    public void addSolver( IMetaheuristic solver ){

        this.getSolvers().add(solver);
        this.getSolversInfo().add( new SolverInfo() );
        
        solver.setCoordinator(this);
        solver.setId( this.getSolvers().size()-1 );
    }

    /**
     * stores the new values in memories and keep sorted.
     * keep the N latest reports
     *
     * @param sol
     * @param improve
     */
    protected void storeInMemory( IMetaheuristicSolution sol, Double improve ){

        if( this.getBestSolution()!=null ){
            if( sol.isBetter(this.getBestSolution()) )
                this.setBestSolution(sol);
        }
        else{
            this.setBestSolution(sol);
        }

        this.getLatestImpr().poll();
        this.getLatestSol().poll();

        if( improve==null )
            improve=0.0;

        this.getLatestImpr().add(improve);
        this.getLatestSol().add(sol);
        
        //copy the lists of latest reports (improve & solution) to the
        // sorted arraylists
        for( int i=0;i<MEMORY_LENGTH;i++ ){
            this.getSolMemory().set(i, getLatestSol().get(i) );
            this.getImproveMemory().set(i, getLatestImpr().get(i) );
        }

        Collections.sort( this.getSolMemory(), new SolutionComparator() );
        Collections.sort( this.getImproveMemory() );
        
    }

    /**
     * get the percentile associated with the improve based in
     * improve memory
     * 
     * @param improve
     * @return
     */
    protected double getPercentile( double improve ){

        for(int i=MEMORY_LENGTH-1;i>=0;i--){
            if( improve > this.getImproveMemory().get(i) ){
                return 100.0*(i+1)/(double)MEMORY_LENGTH;
            }
        }

        return 0;
    }

    /**
     * get the percentile associated with the solution based in
     * solutions memory
     * @param sol
     * @return
     */
    protected double getPercentile( IMetaheuristicSolution sol ){

        for(int i=MEMORY_LENGTH-1;i>=0;i--){

            IMetaheuristicSolution curr=this.getSolMemory().get(i);

            if( curr==null ){
                return 100.0*(i+1)/(double)MEMORY_LENGTH;
            }

            if( !curr.isBetter(sol) ){
                return 100.0*(i+1)/(double)MEMORY_LENGTH;
            }
        }

        return 0;
    }

    /**
     * evaluate the function for fuzzy set "BAD solution/improve"
     *
     * @param sol
     * @return
     */
    protected double isBad( double value ){

        if( value>=BETA )
            return 0.0;

        if( value<BETA && value>ALPHA )
            return (BETA-value)/(BETA-ALPHA);

        return 1.0;
    }


    /**
     * Evaluate this antecedent:
     * IF the quality of the solution reported by solver_i is bad 
     * AND the improvement rate of solver_i is bad
     *
     * @param sol
     * @param improve
     * @param lambda
     * @return
     */
    protected boolean evalFuzzyAntecedent( IMetaheuristicSolution sol,
            double improve, double lambda ){

        double val1 = isBad( getPercentile(sol) );
        double val2 = isBad( getPercentile(improve) );

        if( Math.min(val1, val2) > lambda )
            return true;

        return false;
    }

    /**
     * gets the index of the current best solver
     * 
     * @return
     */
    protected Integer getCurrBestSolver(){

        if( getSolversInfo().size()==0 )
            return null;
        
        int index=0;
        IMetaheuristicSolution bestCurrSol = getSolversInfo().get(0).solution;

        for(int i=1; i<getSolversInfo().size(); i++){

            IMetaheuristicSolution sol=getSolversInfo().get(i).solution;

            if( bestCurrSol==null && sol==null )
                continue;

            if( bestCurrSol==null && sol!=null ){
                bestCurrSol=sol;
                index=i;
            }
            else if( sol.isBetter(bestCurrSol) ){
                bestCurrSol=sol;
                index=i;
            }
        }

        if( bestCurrSol==null )
            return null;

        return index;
    }

    /**
     * gets the best solution
     *
     * @return
     */
    public IMetaheuristicSolution getBestSolution() {
        return this.bestSol;
    }

    /**
     * @param bestSol the bestSol to set
     */
    protected void setBestSolution(IMetaheuristicSolution bestSol) {
        this.bestSol = bestSol;
    }


    /**
     * method used for communication between solvers and coordinator
     *
     * @param id
     * @param sol
     * @param time
     * @return true if a command is sent to the solver, false otherwise
     */
    public synchronized boolean sendReport( int id, IMetaheuristicSolution sol, Date time ){

        //update solver's info
        this.getSolversInfo().get(id).update(sol, time);

        //calculate improve
        Double improve = this.getSolversInfo().get(id).getImprove();

        // store in memories the new solution and improve
        this.storeInMemory(sol, improve);       

        if( improve==null )
            return false;

        boolean val=this.evalFuzzyAntecedent(sol, improve, LAMBDA );

        if( val ){
            //if the antecedent is true build a command to modify the behaviour
            // of the solver
            
            //send the global best solution
            this.getSolvers().get(id).setSolution( 
                    getBestSolution() );

//            //send the lambda of the current best solver
//            Integer idbest=getCurrBestSolver();
//
//            if( idbest!=null ){
//                double bestLambda=this.getSolvers().get(idbest).getLambda();
//                double lambda=this.getSolvers().get(id).getLambda();
//
//                if( lambda < bestLambda ){
//                    this.getSolvers().get(id).setLambda( Math.min( 1.0, lambda+DELTA_LAMBDA) );
//                }
//                else if( lambda > bestLambda ){
//                    this.getSolvers().get(id).setLambda( Math.max( 0.0, lambda-DELTA_LAMBDA) );
//                }
//            }

        }

        return true;
    }


    /**
     * run the solvers of the coordinator in threads and wait until finish all of
     * them
     *
     * @param coordinator
     * @param sleeptime : milliseconds to sleep between checks
     */
    public static void runSolvers( MetaCooperative coordinator, int sleeptime ){

        ExecutorService es =  Executors.newFixedThreadPool( getCores() );

        for( IMetaheuristic solver : coordinator.getSolvers() ){
            es.execute( new MetaheuristicRun(solver) );
        }

        es.shutdown();

        //wait all threads to finish
        while( !es.isTerminated() ){

            try {
                Thread.sleep(sleeptime);
            } catch (InterruptedException ex) {
                Logger.getLogger(MetaCooperative.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
}
