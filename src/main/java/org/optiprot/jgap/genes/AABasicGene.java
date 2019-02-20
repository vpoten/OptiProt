/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.genes;

import org.optiprot.aacids.*;
import org.jgap.*;

/**
 * Class for amino acid JGAP custom genes
 * @author victor
 */
public class AABasicGene extends BaseGene implements IAAGene, Gene, java.io.Serializable {

    private IAABasic AminoAcid;
    
    public AABasicGene(Configuration theConfiguration)
            throws InvalidConfigurationException {

        super(theConfiguration);

    }
    

    @Override
    protected Object getInternalValue() {

        return this.getAminoAcid();
    }

    @Override
    protected Gene newGeneInternal() {
        
        try {
            AABasicGene gene=new AABasicGene(getConfiguration());
            gene.setAminoAcid( (IAABasic)this.getAminoAcid().clone() );
            return gene;
        }
        catch (InvalidConfigurationException e) {
            throw new IllegalStateException(e);
        }
    }

    public void setAllele(Object o) {
        
        if (o instanceof IAABasic) {
            this.setAminoAcid((IAABasic) o);
        }
    }

    public String getPersistentRepresentation() throws UnsupportedOperationException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void setValueFromPersistentRepresentation(String arg0)
            throws UnsupportedOperationException, UnsupportedRepresentationException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void setToRandomValue(RandomGenerator arg0) {
        //without sense to randomize an amino acid with no other information
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void applyMutation(int i, double v) {

        this.getAminoAcid().applyMutation(i, v);
    }

    public int compareTo(Object o) {
        
        if (o instanceof IAABasic) {

            return 0;
        }
        else{
            throw new ClassCastException();
        }

    }

    /**
     * @return the AminoAcid
     */
    public IAABasic getAminoAcid() {
        return AminoAcid;
    }

    /**
     * @param AminoAcid the AminoAcid to set
     */
    public void setAminoAcid(IAABasic AminoAcid) {
        this.AminoAcid = AminoAcid;
    }

}
