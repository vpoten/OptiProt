/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.meme;

import java.util.ArrayList;

/**
 *
 * @author victor
 */
public interface IMemeOperators {

    public void crossover(ArrayList<IMemeSolution> population, ArrayList<IMemeSolution> newPop);

    public void initPopulation(ArrayList<IMemeSolution> population, int nP);

    public void mutation(ArrayList<IMemeSolution> population, ArrayList<IMemeSolution> newPop);

    public void postOperations(ArrayList<IMemeSolution> population);

    public void preOperations(ArrayList<IMemeSolution> population);

    public void selection(ArrayList<IMemeSolution> population, ArrayList<IMemeSolution> newPop);

}
