/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import org.biojava.bio.structure.Atom;
import org.optiprot.potential.element.*;

/**
 * interface for force fields
 * @author victor
 */
public interface IForceField {

    public double[] getParBond(MolecularBond bond);
    public double[] getParAngle(MolecularAngle angle);
    public double[] getParDihedral(MolecularDihedral dihedral);
    public double[] getParImproper(MolecularImproper improper);
    public double[] getParVDW(MolecularNonbonded nonbond);

    public double getKBond(MolecularBond bond);
    public double getEqDistBond(MolecularBond bond);

    public double getKUreyBradley(MolecularAngle angle);
    public double getEqDistUreyBradley(MolecularAngle angle);

    public double getKBondAngle(MolecularAngle angle);
    public double getEqDistBondAngle(MolecularAngle angle);

    public double getKDihedral(MolecularDihedral dihedral);
    public double getNDihedral(MolecularDihedral dihedral);
    public double getDDihedral(MolecularDihedral dihedral);

    public double getKTorsionAngle(MolecularImproper improper);
    public double getEqDistTorsionAngle(MolecularImproper improper);

    public double getVDWRmin(MolecularNonbonded nonbond);
    public double getVDWEner(MolecularNonbonded nonbond);
    public double getInvVDWRadius(Atom atom);
    public double getVDWRadius(Atom a);

    public double getChargeA(MolecularNonbonded nonbond);
    public double getChargeB(MolecularNonbonded nonbond);
    ///public double getChargeAB(MolecularNonbonded nonbond);
    public double getChargeAB(MolecularPair pair);

    public double getKProtLigDielectric();
    public double getKWaterDielectric();
    public double getKSurfTensionWater();

}
