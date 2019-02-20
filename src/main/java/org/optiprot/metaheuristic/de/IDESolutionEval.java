/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de;

/**
 *
 * @author victor
 */
public interface IDESolutionEval {
    
    public Double getFitness(DESolution sol);

    public boolean isBetter(double val1, double val2 );

    public int getNevaluations();
}
