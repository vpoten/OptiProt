/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.neural;

/**
 * Implements Jordan recurrent neural network
 * 
 * @author victor
 */
public class JordanNetwork extends NeuralNetwork {

    double [] delayOutput = null;
    double [] inputs = null;

    public JordanNetwork(int nin, int nout, int nhid) {
        super(nin+nout, nout, nhid);

        inputs=new double [nin+nout];
        delayOutput=new double [nout];
    }

    /**
     * gets the real number of inputs (ninputs-noutputs)
     *
     * @return
     */
    @Override
    public int getRealNInputs() {
        return getNInputs()-getNOutputs();
    }

    @Override
    public void evaluate( double [] inputs ){

        //the final inputs are the pattern and the previous output
        System.arraycopy(inputs, 0, this.inputs, 0, inputs.length);
        System.arraycopy(delayOutput, 0, this.inputs, inputs.length, getNOutputs());

        super.evaluate(this.inputs);
        
        getOutput(delayOutput);
    }

    /**
     * reset delay output
     *
     */
    @Override
    public void init(){
        super.init();
        for(int i=0; i<delayOutput.length; i++){
            delayOutput[i]=0.0;
        }
    }
}
