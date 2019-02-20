/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.fans.FANS;
import org.optiprot.metaheuristic.fans.IFANSConditions;
import org.optiprot.metaheuristic.fans.IFANSOperator;
import org.optiprot.metaheuristic.fans.IFANSOperatorManager;

/**
 *
 * @author victor
 */
public class SecStructFANSOper implements IFANSOperatorManager, IFANSOperator,
    IFANSConditions {

    private double K_STD_INI = 1.0;
    
    private double stdDev = K_STD_INI;///value of perturbation
    private int withoutImprove = 0;//iterations without improve

    private FANS fans = null;

    private int withoutImprLimit = 100;//limit of iterations without improve
    

    DESolution newSol = null;
    

    public SecStructFANSOper(FANS fans) {
        this.fans=fans;

        if( !(fans.getSolution() instanceof DESolution) ){
            throw new ClassCastException("bad SecStructFANSSolution.");
        }

        DESolution sol=(DESolution) fans.getSolution();
        newSol=new DESolution( sol.size() , sol.getEvaluator() );
    }

    public SecStructFANSOper(FANS fans, double stdev_ini) {
        this(fans);
        K_STD_INI=stdev_ini;
    }


    public void modify(IFANSOperator operator) {
        stdDev*=0.5;
    }

    /**
     * modify the solution (IFANSOperator)
     * 
     * @param solution : solution to modify
     * @return a new solution
     */
    public IMetaheuristicSolution modify(IMetaheuristicSolution solution) {

        if( !(solution instanceof DESolution) ){
            throw new ClassCastException("bad SecStructFANSSolution.");
        }

        DESolution currsol=(DESolution) solution;

        for(int i=0; i<currsol.size(); i++){
            newSol.setParameter( i, currsol.getParameter(i)+normal(0, stdDev) );
        }

        return newSol;
    }

    
    public boolean endCondition(FANS fans) {
        return fans.endCondition();
    }

    /**
     * restart the search (IFANSConditions)
     * 
     * @param fans
     */
    public void restart(FANS fans) {

        if( !(fans.getSolution() instanceof DESolution) ){
            throw new ClassCastException("bad SecStructFANSSolution.");
        }

        DESolution sol=(DESolution) fans.getSolution();

        //create a random solution
        ///DESolution newsol = new DESolution( sol.size(), sol.getEvaluator());
        stdDev=K_STD_INI*2;
        fans.setSolution( this.modify(sol) );

        //restart the parameter k and others
        reset(fans);
    }

    public void reset(FANS fans) {
        stdDev = K_STD_INI;
        withoutImprove = 0;
    }

    /**
     * checks if the search is blocked (IFANSConditions)
     *
     * @param fans
     * @return
     */
    public boolean searchBlocked(FANS fans) {

        if( stdDev<1e-6 )
            return true;

        if( withoutImprove>withoutImprLimit )
            return true;

        return false;
    }

    /**
     * update the conditions internals (IFANSConditions)
     *
     * @param fans
     */
    public void update(FANS fans) {
        if( !fans.getSolution().isBetter( fans.getBestSolution()) ){
            withoutImprove++;
        }
        else{
            withoutImprove=0;
        }
    }

    public int getIterations() {
        return fans.getIterations();
    }

    /**
	 * Generate a random number from a Gaussian (Normal) random variable.
	 *
	 * @param mu Mean of the random variable.
	 * @param sigma Standard deviation of the random variable.
	 * @return A double.
	 */
    protected double normal(double mu, double sigma) {
        double x = mu + sigma * Math.cos(2 * Math.PI * Math.random()) * Math.sqrt(-2 * Math.log(Math.random()));
		return x;
    }



}
