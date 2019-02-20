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
public interface IFANSNeighManager {

    public IMetaheuristicSolution generate( IFANSOperator operator,
            IFANSFuzzyEvaluator evaluator, IMetaheuristicSolution solution);

}
