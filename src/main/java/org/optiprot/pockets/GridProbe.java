/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.pockets;

import org.biojava.bio.structure.Atom;
import org.optiprot.maths.AtomGrid;

/**
 *
 * @author victor
 */
public class GridProbe {
    private int index=-1;
    private Atom atom = null;
    private int buriedIdx = 0;

    public static final int NUM_NEIGHBORS = 26;

    public GridProbe( int id, Atom at, int bidx ) {
        index=id;
        atom=at;
        buriedIdx=bidx;
    }

    /**
     * @return the index
     */
    public int getIndex() {
        return index;
    }

    /**
     * @param index the index to set
     */
    public void setIndex(int index) {
        this.index = index;
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
     * @return the buriedIdx
     */
    public int getBuriedIdx() {
        return buriedIdx;
    }

    /**
     * @param buriedIdx the buriedIdx to set
     */
    public void setBuriedIdx(int buriedIdx) {
        this.buriedIdx = buriedIdx;
    }

    /**
     * increases the buried index
     * @param val
     */
    void incBuriedIdx(int val) {
        buriedIdx+=val;
    }

    /**
     * fills the array with the indices of the neighbors in the grid
     * 
     * @param grid
     * @param neighbors : array of 6 grid neighbors
     */
    public void getNeighbors( AtomGrid grid, int [] neighbors ){

       getNeighbors( this.getIndex(), grid, neighbors );

    }

    public static void getNeighbors( int idx, AtomGrid grid, int [] neighbors ){

        int x=grid.getX( idx );
        int y=grid.getY( idx );
        int z=grid.getZ( idx );

        
////        neighbors[0]=grid.getIndexSec( x-1, y, z);
////        neighbors[1]=grid.getIndexSec( x+1, y, z);
////
////        neighbors[2]=grid.getIndexSec( x, y-1, z);
////        neighbors[3]=grid.getIndexSec( x, y+1, z);
////
////        neighbors[4]=grid.getIndexSec( x, y, z-1);
////        neighbors[5]=grid.getIndexSec( x, y, z+1);

        
        neighbors[0]=grid.getIndexSec( x-1, y-1, z-1);
        neighbors[1]=grid.getIndexSec( x-1, y-1, z);
        neighbors[2]=grid.getIndexSec( x-1, y-1, z+1);

        neighbors[3]=grid.getIndexSec( x-1, y, z-1);
        neighbors[4]=grid.getIndexSec( x-1, y, z);
        neighbors[5]=grid.getIndexSec( x-1, y, z+1);

        neighbors[6]=grid.getIndexSec( x-1, y+1, z-1);
        neighbors[7]=grid.getIndexSec( x-1, y+1, z);
        neighbors[8]=grid.getIndexSec( x-1, y+1, z+1);

        neighbors[9]=grid.getIndexSec( x, y-1, z-1);
        neighbors[10]=grid.getIndexSec( x, y-1, z);
        neighbors[11]=grid.getIndexSec( x, y-1, z+1);

        neighbors[12]=grid.getIndexSec( x, y, z-1);
        neighbors[13]=grid.getIndexSec( x, y, z+1);

        neighbors[14]=grid.getIndexSec( x, y+1, z-1);
        neighbors[15]=grid.getIndexSec( x, y+1, z);
        neighbors[16]=grid.getIndexSec( x, y+1, z+1);

        neighbors[17]=grid.getIndexSec( x+1, y-1, z-1);
        neighbors[18]=grid.getIndexSec( x+1, y-1, z);
        neighbors[19]=grid.getIndexSec( x+1, y-1, z+1);

        neighbors[20]=grid.getIndexSec( x+1, y, z-1);
        neighbors[21]=grid.getIndexSec( x+1, y, z);
        neighbors[22]=grid.getIndexSec( x+1, y, z+1);

        neighbors[23]=grid.getIndexSec( x+1, y+1, z-1);
        neighbors[24]=grid.getIndexSec( x+1, y+1, z);
        neighbors[25]=grid.getIndexSec( x+1, y+1, z+1);
        

    }
    
}
