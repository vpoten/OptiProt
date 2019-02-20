/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import org.biojava.bio.structure.Atom;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcIntegrals;
import org.optiprot.potential.IForceField;

/**
 * stores an atom and its residue index
 * @author victor
 */
public class AtomInfo{
    private Atom atom=null;
    private int resIdx = 0;
    private double bornRadii = 0;



    public AtomInfo(Atom at, int ridx) {
        atom=at;
        resIdx=ridx;
    }

    public void calcBornRadii(BSPTree btree, AtomGrid grid, IForceField ffield) {

        ////bornRadii=CalcIntegrals.calcBornRadiiVolInside(atom, grid, btree, ffield);
        setBornRadii(CalcIntegrals.calcBornRadiiVolOutside(getAtom(), grid, btree, true));
    }

    /**
     * @return the atom
     */
    public Atom getAtom() {
        return atom;
    }

    /**
     * @param atom the atom to set
     */
    public void setAtom(Atom atom) {
        this.atom = atom;
    }

    /**
     * @return the resIdx
     */
    public int getResIdx() {
        return resIdx;
    }

    /**
     * @param resIdx the resIdx to set
     */
    public void setResIdx(int resIdx) {
        this.resIdx = resIdx;
    }

    /**
     * @return the bornRadii
     */
    public double getBornRadii() {
        return bornRadii;
    }

    /**
     * @param bornRadii the bornRadii to set
     */
    public void setBornRadii(double bornRadii) {
        this.bornRadii = bornRadii;
    }

}
