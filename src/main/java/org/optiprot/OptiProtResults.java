/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot;

import java.util.Date;
import org.biojava.bio.structure.Chain;

/**
 *
 * @author victor
 */
public class OptiProtResults {
    private Chain bestSolution=null;
    private Chain initSolution = null;
    private double initFitness = 0;
    private double bestFitness = 0;
    private Date startTime = null;
    private Date endTime = null;
    private int iterations=0;

    /**
     * @return the bestSolution
     */
    public Chain getBestSolution() {
        return bestSolution;
    }

    /**
     * @param bestSolution the bestSolution to set
     */
    public void setBestSolution(Chain bestSolution) {
        this.bestSolution = bestSolution;
    }

    /**
     * @return the initSolution
     */
    public Chain getInitSolution() {
        return initSolution;
    }

    /**
     * @param initSolution the initSolution to set
     */
    public void setInitSolution(Chain initSolution) {
        this.initSolution = initSolution;
    }

    /**
     * @return the initFitness
     */
    public double getInitFitness() {
        return initFitness;
    }

    /**
     * @param initFitness the initFitness to set
     */
    public void setInitFitness(double initFitness) {
        this.initFitness = initFitness;
    }

    /**
     * @return the bestFitness
     */
    public double getBestFitness() {
        return bestFitness;
    }

    /**
     * @param bestFitness the bestFitness to set
     */
    public void setBestFitness(double bestFitness) {
        this.bestFitness = bestFitness;
    }

    /**
     * @return the startTime
     */
    public Date getStartTime() {
        return startTime;
    }

    /**
     * @param startTime the startTime to set
     */
    public void setStartTime(Date startTime) {
        this.startTime = startTime;
    }

    /**
     * @return the endTime
     */
    public Date getEndTime() {
        return endTime;
    }

    /**
     * @param endTime the endTime to set
     */
    public void setEndTime(Date endTime) {
        this.endTime = endTime;
    }

    /**
     * @return the iterations
     */
    public int getIterations() {
        return iterations;
    }

    /**
     * @param iterations the iterations to set
     */
    public void setIterations(int iterations) {
        this.iterations = iterations;
    }

}
