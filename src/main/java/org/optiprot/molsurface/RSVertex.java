/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.molsurface;

import java.util.ArrayList;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;

/**
 *
 * @author victor
 */
public class RSVertex extends AtomImpl {

    private double radii = 0;
    private ArrayList<RSEdge> edges=new ArrayList<RSEdge>();
    private boolean free=false;


    public RSVertex( Atom at, double rad) {
        setX( at.getX());
        setY( at.getY());
        setZ( at.getZ());

        setName( at.getName() );
        setFullName( at.getFullName() );
        
        setRadii(rad);
    }


    /**
     * @return the radii
     */
    public double getRadii() {
        return radii;
    }

    /**
     * @param radii the radii to set
     */
    public void setRadii(double radii) {
        this.radii = radii;
    }

    /**
     * @return the edges
     */
    public ArrayList<RSEdge> getEdges() {
        return edges;
    }

    public boolean isFree(){
        return free;
    }

    /**
     * @param free the free to set
     */
    public void setFree(boolean free) {
        this.free = free;
    }


}
