/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.util.Iterator;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.de.IDESolutionEval;
import org.optiprot.neural.NNPattern;
import org.optiprot.neural.NeuralNetwork;

/**
 * DESolution Evaluator for secondary structure prediction neural networks
 *
 * @author victor
 */
public class SecStructDeEvaluator implements IDESolutionEval {

    private NeuralNetwork network = null;

    private int nevaluations=0;

    private SecStructPattManager pattManager=null;

    
    public SecStructDeEvaluator() {
    }

    public SecStructDeEvaluator( NeuralNetwork network,
            SecStructPattManager pattManager ) {
        this.setNetwork(network);
        this.setPattManager(pattManager);
    }
    

    public Double getFitness(DESolution sol) {
        
        double sse = 0;
        nevaluations++;

        getNetwork().setWeights( sol.getParameters() );

        int lim=getPattManager().getNumChains();

        for( int i=0; i<lim; i++){

            Iterator<NNPattern> it=getPattManager().getIterator4Chain(i);
            getNetwork().init();

            while(it.hasNext()){
                NNPattern pat=it.next();
                getNetwork().evaluate( pat.getInputs() );
                sse += getNetwork().calcSSE( pat.getOutputs() );
            }
        }

        return sse;

    }

    public boolean isBetter(double val1, double val2) {
        if(val1<val2)
            return true;

        return false;
    }

    public int getNevaluations() {
        return nevaluations;
    }

    /**
     * @return the network
     */
    protected NeuralNetwork getNetwork() {
        return network;
    }

    /**
     * @param network the network to set
     */
    public void setNetwork(NeuralNetwork network) {
        this.network = network;
    }

    /**
     * @return the pattManager
     */
    public SecStructPattManager getPattManager() {
        return pattManager;
    }

    /**
     * @param pattManager the pattManager to set
     */
    public void setPattManager(SecStructPattManager pattManager) {
        this.pattManager = pattManager;
    }

}
