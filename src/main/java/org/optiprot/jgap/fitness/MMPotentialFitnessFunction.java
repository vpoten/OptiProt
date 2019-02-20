/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.fitness;

import org.biojava.bio.structure.Chain;
import org.jgap.FitnessFunction;
import org.jgap.IChromosome;
import org.optiprot.OptiProtParameters;
import org.optiprot.jgap.chromosome.ProteinChromFactory;
import org.optiprot.potential.MolecularElementsImpl;

/**
 * Use the molecular mechanics function to evaluate the protein
 *
 * @author victor
 */
public class MMPotentialFitnessFunction extends FitnessFunction {

    private OptiProtParameters parameters=null;


    public MMPotentialFitnessFunction( OptiProtParameters par ) {

        this.setParameters(par);
    }


    @Override
    protected double evaluate(IChromosome chromosome) {
        
        Chain chain=ProteinChromFactory.toChain(chromosome, getParameters().isGeneratesH());

        return MolecularElementsImpl.calcEnergy(chain, getParameters().getGrid(),
                getParameters(),
                getParameters().isCalcMMechanics(),
                getParameters().isCalcGBorn(),
                getParameters().isCalcSASA() );
    }

    /**
     * @return the parameters
     */
    public OptiProtParameters getParameters() {
        return parameters;
    }

    /**
     * @param parameters the parameters to set
     */
    public void setParameters(OptiProtParameters parameters) {
        this.parameters = parameters;
    }
}
