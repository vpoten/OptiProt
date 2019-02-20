/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.meme.impl;

import java.util.ArrayList;
import org.optiprot.metaheuristic.meme.IMemeOperators;
import org.optiprot.metaheuristic.meme.IMemeSolution;
import org.optiprot.metaheuristic.tabu.impl.DockConfSolution;
import org.optiprot.potential.docking.Docked;

/**
 * this class implements the operators for memetic algorithm applied
 * to docking
 * 
 * @author victor
 */
public class MemeOperatorDock implements IMemeOperators {
    
    private Docked docked;

    private int MIN_ITERATIONS = 10;
    private int INIT_ITERATIONS = 40;

    private int N_MUTANTS = 3;

    static private double RMSD_EQUALS = 0.75;
    static private int MAX_ITER_MUT = 10;

    private ArrayList<IMemeSolution> mutants=new ArrayList<IMemeSolution>();


    /**
     * do the crossover operations
     *
     * @param population
     * @param newPop
     */
    public void crossover(ArrayList<IMemeSolution> population, ArrayList<IMemeSolution> newPop) {
        return;//nothing to do
    }

    /**
     * create the population randomly and do a local search optimization
     *
     * @param population
     * @param nP : size of population
     */
    public void initPopulation(ArrayList<IMemeSolution> population, int nP) {
        
        for(int i=0;i<nP;i++){

            if( population.size()<nP ){
                //if is the first time then create a new population
                population.add(
                    new DockConfSolution( getDocked(), 
                    DockConfSolution.MODE_OPT_POSE, getInitIterations()) );
            }
            else{
                //reset and re-initialize
                ((DockConfSolution)population.get(i)).reset( getDocked(),
                        getInitIterations() );
            }
        }

        if( mutants.size()<getNMutants() ){
            //if is the first time then create mutants
            for(int i=0;i<getNMutants();i++){
                mutants.add( new DockConfSolution( getDocked(),
                    DockConfSolution.MODE_OPT_POSE, 0) );
            }
        }
        
    }

    /**
     * select some individuals and apply a mutation to get new individuals
     *
     * @param population
     * @param newPop
     */
    public void mutation(ArrayList<IMemeSolution> population, ArrayList<IMemeSolution> newPop) {

        IMemeSolution solution=null;
        int cont=0;
        int numiter=0;

        //generate new population
        while( newPop.size() < getNMutants() ){

            // select a solution randomly
            solution = population.get( (int)Math.floor(Math.random()*population.size()) );

            // copy the selected solution to the mutant
            ((DockConfSolution)mutants.get(cont)).copyParameters((DockConfSolution) solution);

            // generate a solution
            IMemeSolution newsol=mutants.get(cont).mutate();
            numiter++;

            // check if is in the population
            if( !isInPopulation(newsol,population) ){
                newPop.add( newsol );
                cont++;

                //optimize solution
                newsol.optimize(getInitIterations());
            }
            else if( numiter>MAX_ITER_MUT ){
                break;
            }
        }
    }

    /**
     * final operations
     *
     * @param population
     */
    public void postOperations(ArrayList<IMemeSolution> population) {
        //do local search for all individuals
        for( IMemeSolution sol : population )
            sol.optimize(MIN_ITERATIONS);
    }

    /**
     * operations previous to crossover and mutation
     *
     * @param population
     */
    public void preOperations(ArrayList<IMemeSolution> population) {
        //reset mutants population
        for( IMemeSolution sol : mutants )
            ((DockConfSolution)sol).reset( getDocked(), 0 );
    }

    /**
     * implements the selection mechanism:
     *
     * replace the worst solutions
     *
     * @param population : the current population sorted by fitness
     * @param newPop
     */
    public void selection(ArrayList<IMemeSolution> population, ArrayList<IMemeSolution> newPop) {

        int i=0;
        IMemeSolution temp=null;

        for( IMemeSolution newsol : newPop ){
            temp=population.get(i);
            population.set( i, newsol);
            mutants.set(i, temp);
            i++;
        }
    }

    /**
     * @return the docked
     */
    public Docked getDocked() {
        return docked;
    }

    /**
     * @param docked the docked to set
     */
    public void setDocked(Docked docked) {
        this.docked = docked;
    }

    /**
     * checks if the solution is already in the population
     * 
     * @param sol
     * @param population
     * @return
     */
    private boolean isInPopulation(IMemeSolution sol, ArrayList<IMemeSolution> population) {

        ///sol.optimize(5);

        for( IMemeSolution prevsol : population ){

            if( ((DockConfSolution)prevsol).equals(sol) )
                return true;

            ///if( ((DockConfSolution)prevsol).equalConformation( (DockConfSolution)sol,RMSD_EQUALS) )
            ///    return true;
        }

        return false;
    }

    /**
     * @return the INIT_ITERATIONS
     */
    public int getInitIterations() {
        return INIT_ITERATIONS;
    }

    /**
     * @param INIT_ITERATIONS the INIT_ITERATIONS to set
     */
    public void setInitIterations(int INIT_ITERATIONS) {
        this.INIT_ITERATIONS = INIT_ITERATIONS;
    }

    /**
     * @return the N_MUTANTS
     */
    public int getNMutants() {
        return N_MUTANTS;
    }

    /**
     * @param aN_MUTANTS the N_MUTANTS to set
     */
    public void setNMutants(int aN_MUTANTS) {
        N_MUTANTS = aN_MUTANTS;
    }


}
