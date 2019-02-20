/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.chromosome;


import java.util.logging.Level;
import java.util.logging.Logger;
import org.optiprot.aacids.IAABasic;
import org.jgap.*;
import org.optiprot.jgap.OptiProtConfiguration;
import org.optiprot.maths.CalcTransform;


/**
 *
 * @author victor
 */
public class ProteinChromosome extends Chromosome {

    public ProteinChromosome(Configuration theConfiguration)
            throws InvalidConfigurationException {

        super(theConfiguration);

        if( !(theConfiguration instanceof OptiProtConfiguration) ){

            throw new InvalidConfigurationException("Not OptiProtConfiguration");
        }
        
    }

   
    @Override
    public Object clone(){

        ProteinChromosome chr=null;
        try {
            chr = new ProteinChromosome(this.getConfiguration());
            chr.setGenes( this.getGenes().clone() );
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(ProteinChromosome.class.getName()).log(Level.SEVERE, null, ex);
            chr=null;
        }

        return chr;
    }

    /**
     * modify the #aaindex aminoacids's phi angle
     *
     * @param aaindex
     * @param percent
     */
    public void changePhi( int aaindex, double percent ){

        double angle=
                ((OptiProtConfiguration)this.getConfiguration()).getParameters().getMutAngleStepMax()*percent;

        CalcTransform.changePhi( this.toAAarray(), aaindex, angle);
    }

    /**
     * modify the #aaindex aminoacids's psi angle
     * 
     * @param aaindex
     * @param percent
     */
    public void changePsi( int aaindex, double percent ){

        double angle=
                ((OptiProtConfiguration)this.getConfiguration()).getParameters().getMutAngleStepMax()*percent;
       
        CalcTransform.changePsi( this.toAAarray(), aaindex, angle);
    }

    /**
     * cut the chromosomes and adds the appropiate mate's subchain
     *
     * @param cutPoint
     * @param chrMate
     */
    public void oneCutCross( int cutPoint,  ProteinChromosome chrMate )
        throws Exception {

        if( cutPoint==this.getGenes().length-1 ){
            return;
        }

        ProteinChromosome mateCopy=(ProteinChromosome)chrMate.clone();

        CalcTransform.placeSubChain( mateCopy.toAAarray(), this.toAAarray(),
                cutPoint, cutPoint, this.getGenes().length );

        for( int i=cutPoint+1; i<this.getGenes().length;i++){
            this.getGenes()[i]=mateCopy.getGenes()[i];
        }

    }

    /**
     * returns the alleles (getGenes()) in a IAABasic array
     * 
     * @return : the array of IAABasic
     */
    public IAABasic[] toAAarray(){

        Gene [] genes=this.getGenes();

        IAABasic [] array= new IAABasic[genes.length];

        for(int i=0;i<array.length;i++){
            array[i]=(IAABasic) genes[i].getAllele();
        }

        return array;
    }

}
