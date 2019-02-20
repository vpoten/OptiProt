/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;

/**
 *
 * @author victor
 */
public class MolecularTorsional {

    private Atom [] atoms=new Atom [4];
    private int [] atomTypes=new int [4];
    private double [] atomCharges=new double [4];
    protected double m_angle=0;//angle in degrees

    
    
    public void setAtom(int idx, Atom at, double charge, int type){

        this.getAtoms()[idx]=at;
        this.atomCharges[idx]=charge;
        this.atomTypes[idx]=type;

    }

    public int getAtomType( int index ){

        return this.atomTypes[index];
    }

    public double getAtomCharge( int index ){

        return this.atomCharges[index];
    }

    /**
     *
     * @return : the dihedral angle in degrees between [-PI,PI]
     */
    public double getAngle() {
        return m_angle;
    }

    public void calcDihedralAng() throws StructureException {

        m_angle = Calc.torsionAngle(
                this.getAtoms()[0],
                this.getAtoms()[1],
                this.getAtoms()[2],
                this.getAtoms()[3] );

        if( Double.isNaN(m_angle) )
            m_angle=180;
    }

    /**
     * @return the atoms
     */
    protected Atom[] getAtoms() {
        return atoms;
    }

    public Atom getAtom( int idx ) {
        return this.atoms[idx];
    }

}
