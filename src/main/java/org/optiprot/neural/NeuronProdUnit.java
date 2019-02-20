/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.neural;

import org.optiprot.maths.CalcGeom;

/**
 * Product Unit Neuron
 * 
 * @author victor
 */
public class NeuronProdUnit extends Neuron {

    public NeuronProdUnit(int ninputs) {
        super(ninputs);
    }

    /**
     * calculates the output product: Prod(x[i]^w[i])
     *
     * @param inputs
     * @return
     */
    @Override
    public double output( double [] inputs ){

        double prod=1;

        for(int i=0; i<this.weights.length; i++){
            prod *= CalcGeom.pow( inputs[i], this.weights[i]);
        }
        
        return prod;
    }

    
}
