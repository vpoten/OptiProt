/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import org.biojava.bio.structure.Atom;
import org.optiprot.potential.element.MolecularAngle;
import org.optiprot.potential.element.MolecularBond;
import org.optiprot.potential.element.MolecularDihedral;
import org.optiprot.potential.element.MolecularImproper;
import org.optiprot.potential.element.MolecularNonbonded;
import org.optiprot.potential.element.MolecularPair;
import org.optiprot.rotamer.RotamerLibrary;

/**
 * Force field for some test cases
 *
 * @author victor
 */
public class TestForceField implements IForceField {

    public double getKBond(MolecularBond bond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getEqDistBond(MolecularBond bond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getKUreyBradley(MolecularAngle angle) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getEqDistUreyBradley(MolecularAngle angle) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getKBondAngle(MolecularAngle angle) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getEqDistBondAngle(MolecularAngle angle) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getKDihedral(MolecularDihedral dihedral) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getNDihedral(MolecularDihedral dihedral) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getDDihedral(MolecularDihedral dihedral) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getKTorsionAngle(MolecularImproper improper) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getEqDistTorsionAngle(MolecularImproper improper) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getVDWRmin(MolecularNonbonded nonbond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getVDWEner(MolecularNonbonded nonbond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getVDWRadius(Atom a) {

        char symbol='\0';

        if( RotamerLibrary.isHidrogen(a) ){
            symbol='H';
        }
        else{
            symbol=a.getName().charAt(0);
        }

        switch(symbol) {
           case 'C':
                return 1.70;
           case 'N':
                return 1.55;
           case 'O':
                return 1.52;
           case 'H':
                return 1.09;
           case 'S':
                return 1.80;
           default:
                return 2.00;
        }

    }

    public double getInvVDWRadius(Atom a) {

        char symbol='\0';

        if( RotamerLibrary.isHidrogen(a) ){
            symbol='H';
        }
        else{
            symbol=a.getName().charAt(0);
        }

        switch(symbol) {
           case 'C':
                return 0.58823;
           case 'N':
                return 0.64516;
           case 'O':
                return 0.65789;
           case 'H':
                return 0.91743;
           case 'S':
                return 0.55555;
           default:
                return 0.5;
        }

    }

    public double getChargeA(MolecularNonbonded nonbond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getChargeB(MolecularNonbonded nonbond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

//    public double getChargeAB(MolecularNonbonded nonbond) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }

    public double getChargeAB(MolecularPair pair) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double getKProtLigDielectric() {
        return 1.0;
    }

    public double getKWaterDielectric() {
        return 78.4;
    }

    public double getKSurfTensionWater() {
        return 0.0072;
    }

    public double[] getParBond(MolecularBond bond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double[] getParAngle(MolecularAngle angle) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double[] getParDihedral(MolecularDihedral dihedral) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double[] getParImproper(MolecularImproper improper) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double[] getParVDW(MolecularNonbonded nonbond) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
