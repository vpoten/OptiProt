/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import org.optiprot.neural.NNPattern;

/**
 * Interface for secondary structure NN Pattern generators
 *
 * @author victor
 */
public interface ISecStructPatternGen {

    public ISecStructPatternGen clone();

    /**
     * Creates a pattern that fits with the windowSize
     *
     * @param windowSize
     * @return
     */
    public NNPattern createPattern( int windowSize );

    /**
     * Generates a pattern using the window and the position of predicted
     * residue inside the window
     * 
     * @param window
     * @param resPredPos
     * @param pat : I/O generated pattern
     */
    public void genPattern( SecStructSubPattern []window, int resPredPos, NNPattern pat);

    public NNPattern genPattern( SecStructSubPattern []window, int resPredPos);
}
