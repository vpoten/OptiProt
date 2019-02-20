/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;


import java.util.List;
import org.biojava.bio.structure.Chain;
import org.optiprot.potential.element.*;

/**
 * interface for molecular elements containers
 * @author victor
 */
public interface IMolecularElements {

    public boolean hasAngles();
    public boolean hasBonds();
    public boolean hasDihedrals();
    public boolean hasImpropers();
    public boolean hasNonbondeds();
    public boolean hasPairs();

    public List<MolecularAngle> getAngles();
    public List<MolecularBond> getBonds();
    public List<MolecularDihedral> getDihedrals();
    public List<MolecularImproper> getImpropers();
    public List<MolecularNonbonded> getNonbondeds();
    public List<MolecularPair> getPairs();

    public Chain getChain();
    
    
}
