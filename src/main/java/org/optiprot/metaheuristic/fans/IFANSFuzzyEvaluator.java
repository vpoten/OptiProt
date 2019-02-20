/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans;

import org.optiprot.metaheuristic.IMetaheuristicSolution;

/**
 *
 * @author victor
 */
public interface IFANSFuzzyEvaluator {

    public void adjust(IMetaheuristicSolution newSolution);

    public double getLambda();

    public boolean isGood(IMetaheuristicSolution newSolution);

    public void setLambda( double lambda );

}
