/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de;

/**
 * interface for crossover operator of DE
 *
 * @author victor
 */
public interface IDECrossover {

    public void cross(DESolution trial, DESolution target, DESolution mutant, double CR);

}
