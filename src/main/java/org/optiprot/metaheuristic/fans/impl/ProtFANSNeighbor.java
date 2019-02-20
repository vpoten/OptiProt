/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans.impl;

import org.optiprot.maths.CalcClashes;
import org.optiprot.metaheuristic.fans.IFANSFuzzyEvaluator;
import org.optiprot.metaheuristic.fans.IFANSOperator;
import org.optiprot.metaheuristic.IMetaheuristicSolution;

/**
 * implements FuzzyEvaluator and NeighborManager for proteins
 *
 * @author victor
 */
public class ProtFANSNeighbor extends FANSNeighborBase {

    static private int maxFixTrials = 100;

    /////////////////////

    /**
     *
     * @param p_acceptRate : the accept rate (between [0,1])
     * @param p_maxTrials : neighbor manager max trials to find a solution
     */
    public ProtFANSNeighbor( double p_acceptRate, int p_maxTrials, double initFitness) {
        super(p_acceptRate, p_maxTrials, initFitness);
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
    @Override
    public IMetaheuristicSolution generate( IFANSOperator operator,
            IFANSFuzzyEvaluator evaluator, IMetaheuristicSolution solution) {

        ProtFANSSolution newSol=null;

        for( int i=0; i<maxTrials; i++){
            newSol=(ProtFANSSolution) operator.modify(solution);

            // fix the modified chain
            CalcClashes.fixClashes( 
                    Math.max( 0.25, 1.0-getLambda()),
                    newSol.getChain(),
                    newSol.getParameters().getRotLib(), maxFixTrials );

            if( evaluator.isGood(newSol) )
                return newSol;

            newSol=null;
        }

        return newSol;
    }


}
