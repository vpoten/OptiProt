/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.neural;

/**
 * Neural Network with Product Unit hidden layer
 * 
 * @author victor
 */
public class ProdUnitNetwork extends NeuralNetwork {


    public ProdUnitNetwork(int nin, int nout, int nhid) {
        this.setNInputs(nin);
        this.setNOutputs(nout);
        this.setNHidden(nhid);

        for(int i=0; i<getNHidden(); i++){
            this.getHiddenLayer().add( new NeuronProdUnit(getNInputs()) );
        }

        outsHidden=new double [getNHidden()];

        for(int i=0; i<getNOutputs(); i++){
            this.getOutputLayer().add( new Neuron(this.getNHidden()) );
        }

        outputs=new double [getNOutputs()];
    }

    /**
     * evaluate the network over the inputs
     *
     * @param inputs
     */
    @Override
    public void evaluate( double [] inputs ){

        for(int i=0; i<getNHidden(); i++){
            outsHidden[i]=this.getHiddenLayer().get(i).output(inputs);
        }

        for(int i=0; i<getNOutputs(); i++){
            outputs[i]=this.getOutputLayer().get(i).output(outsHidden);
            outputs[i]=transferFunction(outputs[i]);
        }
    }

}
