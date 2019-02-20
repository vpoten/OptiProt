/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans.impl;

import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.metaheuristic.fans.IFANSFuzzyEvaluator;
import org.optiprot.metaheuristic.fans.IFANSNeighManager;
import org.optiprot.metaheuristic.fans.IFANSOperator;

/**
 *
 * @author victor
 */
public class FANSNeighborBase implements IFANSFuzzyEvaluator, IFANSNeighManager {

    protected double lambda=0;
    protected double currentFitness=java.lang.Double.MAX_VALUE*0.5;
    protected double acceptRate=0;
    protected int maxTrials=0;

    /////////////////////

    /**
     *
     * @param p_acceptRate : the accept rate (between [0,1])
     * @param p_maxTrials : neighbor manager max trials to find a solution
     */
    public FANSNeighborBase( double p_acceptRate, int p_maxTrials, double initFitness) {
        acceptRate=p_acceptRate;
        maxTrials=p_maxTrials;
        currentFitness=initFitness;
    }

    public FANSNeighborBase( double p_acceptRate, int p_maxTrials) {
        acceptRate=p_acceptRate;
        maxTrials=p_maxTrials;
    }


    /**
     * adjust the fuzzy evaluator values to new solution
     *
     * @param newSolution
     */
    public void adjust(IMetaheuristicSolution newSolution) {
        currentFitness=newSolution.getFitness();
    }

    /**
     * returns true if the new solution is good
     *
     * @param newSolution
     * @return
     */
    public boolean isGood(IMetaheuristicSolution newSolution) {

        double newFitness = newSolution.getFitness();

        if( newFitness < currentFitness )
            return true;

        double beta = currentFitness*(1.0+acceptRate);

        if( (newFitness >= currentFitness) && (newFitness < beta) ){
            double val = (newFitness-beta)/(currentFitness-beta);

            if( val>=getLambda() )
                return true;
        }

        return false;
    }

    /**
     * set the fuzzy acceptability parameter
     *
     * @param p_lambda
     */
    public void setLambda(double p_lambda) {
        lambda=p_lambda;
    }

    /**
     * get the fuzzy acceptability parameter
     * @return
     */
    public double getLambda() {
        return lambda;
    }


    /**
     * generate a new solution (IFANSNeighManager)
     * neighbor manager type First: return the first S' / u(S,S')>=lambda
     *
     * @param operator
     * @param evaluator
     * @param solution
     * @return : null if doesnt find a good solution, a solution otherwise
     */
    public IMetaheuristicSolution generate( IFANSOperator operator,
            IFANSFuzzyEvaluator evaluator, IMetaheuristicSolution solution) {

        IMetaheuristicSolution newSol=null;

        for( int i=0; i<maxTrials; i++){
            newSol = operator.modify(solution);

            if( evaluator.isGood(newSol) )
                return newSol;

            newSol=null;
        }

        return newSol;
    }
    
}
