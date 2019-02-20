/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap;

import org.optiprot.jgap.operator.OneCutCrossoverOperator;
import org.optiprot.jgap.operator.RotamerMutationOperator;
import org.optiprot.jgap.operator.AngleMutationOperator;
import org.optiprot.*;
import org.jgap.*;
import org.jgap.impl.*;

/**
 *
 * @author victor
 */
public class OptiProtConfiguration extends Configuration {

    private OptiProtParameters parameters=null;

    public OptiProtConfiguration( OptiProtParameters par )
        throws InvalidConfigurationException {
        super();

        reset();

        this.setParameters(par);

        setName( par.getName() );

        // TODO constructor
        setRandomGenerator( new StockRandomGenerator() );

        setPopulationSize( par.getPopulationSize() );

        addNaturalSelector( new WeightedRouletteSelector(this), false );
        //addNaturalSelector( new BestChromosomesSelector(this), false );
        setPreservFittestIndividual(true);
        setKeepPopulationSizeConstant(false);

        addGeneticOperator( new OneCutCrossoverOperator(this, par.getCrossOneCutRate()) );
        addGeneticOperator( new AngleMutationOperator(this, par.getMutRateAngle()) );
        addGeneticOperator( new RotamerMutationOperator(this, par.getMutRateRotamer()) );

        setFitnessEvaluator( new DeltaFitnessEvaluator() );
        //setFitnessFunction( new TinkerFitnessFunction(par) );
        
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
