/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans.impl;

import org.optiprot.aacids.AABasicFactory;
import org.biojava.bio.structure.Chain;
import org.optiprot.OptiProtParameters;
import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.aacids.IAABasic;
import org.optiprot.maths.CalcTransform;
import org.optiprot.potential.MolecularElementsImpl;

/**
 * implements FANS solution for proteins
 *
 * @author victor
 */
public class ProtFANSSolution implements IMetaheuristicSolution {

    private Double fitness=null;
    private IAABasic[] chain=null;
    private OptiProtParameters parameters=null;

    public ProtFANSSolution( IAABasic[] p_chain, OptiProtParameters p_parameters){
        chain=p_chain;
        parameters=p_parameters;
    }

    /**
     * center the chain
     */
    public void center() {
        CalcTransform.center( getChain() );
    }

    /**
     * 
     * @return
     */
    public Double getFitness() {

        if( fitness==null ){
            Chain l_chain=AABasicFactory.toChain( chain, getParameters().isGeneratesH());
            fitness =
                    MolecularElementsImpl.calcEnergy(l_chain, getParameters().getGrid(),
                        getParameters(),
                        getParameters().isCalcMMechanics(),
                        getParameters().isCalcGBorn(),
                        getParameters().isCalcSASA() );
            l_chain=null;
        }

        return fitness;
    }

    /**
     * @return the chain
     */
    public IAABasic[] getChain() {
        return chain;
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

    /**
     *
     * @param other : another solution
     * @return : true if this is better than other
     */
    public boolean isBetter(IMetaheuristicSolution other) {

        if( !(other instanceof ProtFANSSolution) ){
            throw new ClassCastException("bad ProtFANSSolution.");
        }

        if( this.getFitness()<other.getFitness() )
            return true;

        return false;
    }

    
}
