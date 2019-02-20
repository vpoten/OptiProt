/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.Atom;

/**
 * class for aacid clashes
 * 
 * @author victor
 */
public class AAClash {

    public Integer idx1=null;
    public Integer idx2=null;
    public Atom atom1=null;
    public Atom atom2=null;
    public double sqDistLim=0;


    public AAClash(int id1, Atom at1, int id2, Atom at2, double sqlimit) {
        idx1=id1;
        idx2=id2;
        atom1=at1;
        atom2=at2;
        sqDistLim=sqlimit;
    }

    public boolean isClashed( double lambda ){

        return CalcClashes.isClashed(
            CalcGeom.squareDistance(atom1,atom2),
            sqDistLim, CalcClashes.acceptRate,
            lambda );
    }

    public int distance(){
        return idx2-idx1;
    }
}
