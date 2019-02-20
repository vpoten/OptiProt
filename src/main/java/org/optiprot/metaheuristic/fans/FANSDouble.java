/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans;

import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.de.IDESolutionEval;

/**
 * FANS for DESolution (vector of doubles), optimize memory
 * reusing solutions
 *
 * @author victor
 */
public class FANSDouble extends FANS {


    public FANSDouble( int size, IDESolutionEval eval ) {
        this.bestSolution=new DESolution(size, eval);
        this.solution=new DESolution(size, eval);
    }


    /**
     * @param bestSolution the bestSolution to set
     */
    @Override
    protected void setBestSolution(IMetaheuristicSolution bestSolution) {
        ((DESolution)this.bestSolution).copy( (DESolution) bestSolution);
    }


    /**
     * @param solution the current solution to set
     */
    @Override
    public void setSolution(IMetaheuristicSolution solution) {

        ((DESolution)this.solution).copy( (DESolution) solution);

        if( this.getSolution().isBetter( this.getBestSolution() ) ){
            this.setBestSolution(this.getSolution());
        }
    }

    /**
     * reset the solutions
     */
    @Override
    public void reset() {
        this.getConditions().reset(this);
        ((DESolution)this.getSolution()).reset();
        ((DESolution)this.getBestSolution()).reset();
    }
}
