/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.tabu;

import java.util.Collection;
import org.optiprot.metaheuristic.IMetaheuristicSolution;

/**
 *
 * @author victor
 */
public interface ITabuNeighManager {

    public IMetaheuristicSolution generate( IMetaheuristicSolution solution,
            Collection<IMetaheuristicSolution> tabuList );

}
