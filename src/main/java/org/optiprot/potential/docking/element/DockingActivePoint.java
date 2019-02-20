/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import org.biojava.bio.structure.Atom;
import org.optiprot.potential.docking.element.DockingAtomClassify.simpletype;

/**
 *
 * @author victor
 */
public class DockingActivePoint {
    private Atom point=null;
    private simpletype ligandType = null;

    public DockingActivePoint() {
    }

    public DockingActivePoint( Atom p, simpletype t) {
        setPoint(p);
        setLigandType(t);
    }


    /**
     * @return the point
     */
    public Atom getPoint() {
        return point;
    }

    /**
     * @param point the point to set
     */
    public void setPoint(Atom point) {
        this.point = point;
    }

    /**
     * @return the ligandType
     */
    public simpletype getLigandType() {
        return ligandType;
    }

    /**
     * @param ligandType the ligandType to set
     */
    public void setLigandType(simpletype ligandType) {
        this.ligandType = ligandType;
    }

    /**
     * checks if the types are compatible (equals)
     * 
     * @param t1
     * @param t2
     * @return
     */
    static public boolean isTypeComp( simpletype t1, simpletype t2 ){

        if( t1==t2 )
            return true;

        if( t1==simpletype.acceptor && t2==simpletype.donor_accept )
            return true;

        if( t2==simpletype.acceptor && t1==simpletype.donor_accept )
            return true;

        return false;
    }

}
