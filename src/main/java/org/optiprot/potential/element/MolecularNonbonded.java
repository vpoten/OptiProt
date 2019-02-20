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
public class MolecularNonbonded extends MolecularPair {

    public MolecularNonbonded() {
        super();
    }

    public MolecularNonbonded(Atom at1, Atom at2) {
        super(at1,at2);
    }
    
}
