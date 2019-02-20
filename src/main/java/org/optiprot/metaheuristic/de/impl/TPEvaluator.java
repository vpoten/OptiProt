/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de.impl;

import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.de.IDESolutionEval;

/**
 * Chebychev polynomial TN evaluator
 * 
 * @author victor
 */
public class TPEvaluator implements IDESolutionEval {

    int      evaluation_samples = 60;
    double   lower_limit = 5.9;
    int dim=5;
    double [] temp=null;
    private int nevaluations=0;

    public TPEvaluator(int p_dim){

        dim=p_dim;

        if(dim==9){
            lower_limit = 72.661;
        }

        temp=new double [dim];
    }

    /**
     * The actual objective function consists of the sum of squared
     * errors, where an error is the magnitude of deviation of the
     * polynomial at a specific argument value.		
     *
     * @param sol
     * @return
     */
    public Double getFitness(DESolution sol) {

        nevaluations++;

        for(int i=0;i<dim;i++)
            temp[i]=sol.getParameter(i);

        double y = 0.0;
        double x = -1.0;
        double z = 0.0, aux;

        double dx = 2 / ((double) evaluation_samples);
        for (int i = 0;  i <= evaluation_samples;  i++, x += dx)
        {
          if ((z = polynomial (temp, x, dim)) > 1.0)
          {
            aux = 1.0 - z;
            y += aux * aux;
          }
          else if (z < -1.0)
          {
            aux = z - 1.0;
            y += aux * aux;
          }
        }

        aux = lower_limit - z;
        aux *= aux;

        if (polynomial (temp, -1.2, dim) < lower_limit)
          y += aux;

        if (polynomial (temp, +1.2, dim) < lower_limit)
          y += aux;

        return y;
    }

    public int getNevaluations() {
        return nevaluations;
    }

    /**
     * 
     * @param val1
     * @param val2
     * @return
     */
    public boolean isBetter(double val1, double val2) {
         if(val1<val2)
            return true;

        return false;
    }

    /**
     * Evaluate the current polynomial.
     * 
     * @param temp
     * @param x
     * @param dim
     * @return
     */
    private double polynomial (double temp[], double x, int dim){
        double y = temp[0];

        for (int j = 1;  j < dim;  j++)
          y = x * y + temp[j];

        return y;
   }

}
