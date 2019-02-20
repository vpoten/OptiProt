/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.meme;

import org.optiprot.metaheuristic.IMetaheuristicSolution;

/**
 *
 * @author victor
 */
public interface IMemeSolution extends IMetaheuristicSolution {

    public void optimize(int evolutions);

    public IMemeSolution mutate();
}
