/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.meme;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import org.optiprot.metaheuristic.MetaheuristicBase;

/**
 *
 * @author victor
 */
public class Memetic extends MetaheuristicBase {

    // population
    private ArrayList<IMemeSolution> population=new ArrayList<IMemeSolution>();
    private int NP = 10;

    private IMemeOperators operators=null;

    
    public void doSearch() {

        if( getOperators()==null ){
            throw new RuntimeException("Memetic bad initialization.");
        }

        ArrayList<IMemeSolution> newPop=new ArrayList<IMemeSolution>();

        setStartTime( new Date() );

        getOperators().initPopulation( getPopulation(), getNP() );

        //sort population by fitness
        Collections.sort(getPopulation(), solComp);

        while( !endCondition() ){

            newPop.clear();

            getOperators().preOperations( getPopulation() );

            getOperators().crossover( getPopulation(), newPop );

            getOperators().mutation( getPopulation(), newPop );

            getOperators().selection( getPopulation(), newPop );

            getOperators().postOperations( getPopulation() );
            
            //sort population by fitness
            Collections.sort(getPopulation(), solComp);

            //get the current solution
            this.setSolution( getPopulation().get(getPopulation().size()-1) );

        }

        setEndTime( new Date() );
    }


    public double getLambda() {
        return 1.0;//nothing to do
    }

    public void setLambda(double val) {
        return;//nothing to do
    }

    /**
     * @return the NP
     */
    public int getNP() {
        return NP;
    }

    /**
     * @param NP the NP to set
     */
    public void setNP(int NP) {
        this.NP = NP;
    }

    /**
     * @return the population
     */
    public ArrayList<IMemeSolution> getPopulation() {
        return population;
    }

    /**
     * @param population the population to set
     */
    protected void setPopulation(ArrayList<IMemeSolution> population) {
        this.population = population;
    }

    /**
     * @return the operators
     */
    public IMemeOperators getOperators() {
        return operators;
    }

    /**
     * @param operators the operators to set
     */
    public void setOperators(IMemeOperators operators) {
        this.operators = operators;
    }

}
