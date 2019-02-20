/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.tabu.impl;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.metaheuristic.SolutionComparator;
import org.optiprot.metaheuristic.tabu.ITabuNeighManager;

/**
 * this class implements the neighbor manager for tabu search applied
 * to docking
 *
 * @author victor
 */
public class DockTabuNeighbor implements ITabuNeighManager {

    private int OPT_STEPS = 30;
    static private final int MAX_ITER = 12;
    private int numNeighbors=4;
    
    private ArrayList<DockConfSolution> neighbors=new ArrayList<DockConfSolution>();
    SolutionComparator solComp=new SolutionComparator();


    
    public IMetaheuristicSolution generate(IMetaheuristicSolution solution,
            Collection<IMetaheuristicSolution> tabuList) {

        if( !(solution instanceof DockConfSolution) ){
            throw new ClassCastException("bad DockConfSolution.");
        }

        DockConfSolution newsol=null;

        if( !neighbors.isEmpty() )
            neighbors.clear();

        int numiter=0;

        //generate neighbors
        while( neighbors.size() < getNumNeighbors() ){

            // generate a solution
            newsol=((DockConfSolution)solution).newNeighbor();
            numiter++;

            // check if is in tabu list
            if( !isInTabuList(newsol, tabuList) ){
                neighbors.add(newsol);
                
                //optimize solution
                newsol.optimize(getOptSteps());
            }
            else if( numiter>MAX_ITER ){
                break;
            }
        }

        if( neighbors.isEmpty() ){
            //if cannot find new neighbors return the current one
            ((DockConfSolution)solution).optimize(getOptSteps());
            return solution;
        }

        //get the best neighbor
        Collections.sort( neighbors, solComp);

        newsol=neighbors.get( neighbors.size()-1 );
        newsol.optimize(getOptSteps()*2);

        return newsol;
    }

    /**
     * @return the numNeighbors to explore
     */
    public int getNumNeighbors() {
        return numNeighbors;
    }

    /**
     * @param numNeighbors the numNeighbors to set
     */
    public void setNumNeighbors(int numNeighbors) {
        this.numNeighbors = numNeighbors;
    }

    /**
     * checks if the solution is in tabu list
     * 
     * @param sol
     * @param tabuList
     * @return
     */
    private boolean isInTabuList( DockConfSolution sol, Collection<IMetaheuristicSolution> tabuList ){

        for( IMetaheuristicSolution prevsol : tabuList ){

            if( ((DockConfSolution)prevsol).equals(sol) )
                return true;
        }

        return false;
    }

    /**
     * @return the OPT_STEPS
     */
    public int getOptSteps() {
        return OPT_STEPS;
    }

    /**
     * @param OPT_STEPS the OPT_STEPS to set
     */
    public void setOptSteps(int OPT_STEPS) {
        this.OPT_STEPS = OPT_STEPS;
    }

}
