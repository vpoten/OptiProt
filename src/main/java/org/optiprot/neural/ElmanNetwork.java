/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.neural;

/**
 * Implements Elman recurrent neural network
 *
 * @author victor
 */
public class ElmanNetwork extends NeuralNetwork {

    double [] delayHiddOut = null;
    double [] inputs = null;

    public ElmanNetwork(int nin, int nout, int nhid) {
        super(nin+nhid, nout, nhid);

        inputs=new double [nin+nhid];
        delayHiddOut=new double [nhid];
    }

    /**
     * gets the real number of inputs (ninputs-noutputs)
     *
     * @return
     */
    @Override
    public int getRealNInputs() {
        return getNInputs()-getNHidden();
    }

    @Override
    public void evaluate( double [] inputs ){

        //the final inputs are the pattern and the previous output
        System.arraycopy(inputs, 0, this.inputs, 0, inputs.length);
        System.arraycopy(delayHiddOut, 0, this.inputs, inputs.length, getNHidden());

        super.evaluate(this.inputs);

        getHiddenOutput(delayHiddOut);
    }

    /**
     * reset delay output
     *
     */
    @Override
    public void init(){
        super.init();
        for(int i=0; i<delayHiddOut.length; i++){
            delayHiddOut[i]=0.0;
        }
    }
}
