/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.fitness;

import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.optiprot.OptiProtParameters;
import org.optiprot.potential.TinkerPotentialEnergy;
import org.jgap.FitnessFunction;
import org.jgap.IChromosome;
import org.optiprot.jgap.chromosome.ProteinChromFactory;
import org.optiprot.io.XyzStructureWriter;

/**
 * Use tinker program (analize) to evaluate the energy
 *
 * @author victor
 */
public class TinkerFitnessFunction extends FitnessFunction {

    private OptiProtParameters parameters=null;
    
    
    public TinkerFitnessFunction( OptiProtParameters par ){

        this.setParameters(par);

        TinkerPotentialEnergy.setForceField( this.getParameters().getTinkerForceField() );
        TinkerPotentialEnergy.setTinkerPath( this.getParameters().getTinkerPath() );
        XyzStructureWriter.setForceField( this.getParameters().getTinkerForceField() );
        XyzStructureWriter.setTinkerPath( this.getParameters().getTinkerPath() );
    }

    @Override
    protected double evaluate(IChromosome chromosome) {

        String path=getParameters().getWorkDir() + File.separator +
                chromosome.getConfiguration().getName();

        String xyzFile=path+".xyz";
        String pdbFile=path+".pdb";

        try{
            XyzStructureWriter.writeStructurePDB(
                    ProteinChromFactory.toChain(chromosome, getParameters().isGeneratesH()), pdbFile);
            XyzStructureWriter.writeStructureXyz( pdbFile );
        }
        catch(Exception ex){
            Logger.getLogger(TinkerFitnessFunction.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        return TinkerPotentialEnergy.calcEnergy( xyzFile );

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
