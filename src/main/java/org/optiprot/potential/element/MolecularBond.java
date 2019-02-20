 /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import org.biojava.bio.structure.Atom;

/**
 *
 * @author victor
 */
public class MolecularBond extends MolecularPair {

    
    public enum BondType {SIMPLE,DOUBLE,TRIPLE,AROMATIC};

    private BondType bondType=BondType.SIMPLE;

    public MolecularBond() {
        super();
    }

    
    public MolecularBond(Atom at1, Atom at2) {
        super(at1,at2);
    }

    /**
     * @return the bondType
     */
    public BondType getBondType() {
        return bondType;
    }

    /**
     * @param bondType the bondType to set
     */
    public void setBondType(BondType type) {
        this.bondType = type;
    }

}
