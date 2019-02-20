/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans;

import org.optiprot.metaheuristic.IMetaheuristicSolution;

/**
 * Interface for FANS operators
 * @author victor
 */
public interface IFANSOperator {

    public IMetaheuristicSolution modify( IMetaheuristicSolution solution );
}
