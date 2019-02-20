/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic;

/**
 * class that runs a solver (IMetaheuristic)
 * 
 * @author victor
 */
public class MetaheuristicRun implements Runnable {

    private IMetaheuristic solver=null;

    public MetaheuristicRun(IMetaheuristic solver) {
        this.solver=solver;
    }

    public void run() {
        solver.doSearch();
    }

}
