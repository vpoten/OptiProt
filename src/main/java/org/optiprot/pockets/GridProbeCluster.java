/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.pockets;

import java.util.ArrayList;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.CalcGeom;

/**
 *
 * @author victor
 */
public class GridProbeCluster {

    private int [] neighbors=new int [GridProbe.NUM_NEIGHBORS];
    
    private ArrayList<GridProbe> probes=new ArrayList<GridProbe>();

    /**
     * @return the probes
     */
    public ArrayList<GridProbe> getProbes() {
        return probes;
    }

    /**
     * 
     * @param probe
     * @param grid
     * @return : true if the probe is adjoin to the cluster
     */
     public boolean isAdjoin(GridProbe probe, AtomGrid grid) {

        probe.getNeighbors(grid, neighbors);

        for( GridProbe cprobe : getProbes() ){

            for(int i=0;i<neighbors.length;i++){
                if( neighbors[i]==cprobe.getIndex() )
                    return true;
            }
        }

        return false;
    }

     /**
      *
      * @param cluster
      * @param grid
      * @return : true if both clusters are adjoined
      */
    public boolean isAdjoin(GridProbeCluster cluster, AtomGrid grid) {

        for( GridProbe cprobe : cluster.getProbes() ){

            if( this.isAdjoin(cprobe, grid) )
                return true;
        }

        return false;
    }

    /**
     * join both clusters and clears the cluster passed by parameter
     *
     * @param cluster
     */
    public void join( GridProbeCluster cluster ){
        this.getProbes().addAll( cluster.getProbes() );
        cluster.getProbes().clear();
    }

    /**
     * fills the list with the points of the cluster
     * 
     * @param list
     */
    public void getAtoms( List<Atom> list ){

        for( GridProbe cprobe : getProbes() ){
            list.add( cprobe.getAtom() );
        }
    }

    /**
     * clacls the centroid of the cluster
     *
     * @return
     */
    public Atom getCentroid(){

        ArrayList<Atom> listAtom=new ArrayList<Atom>();
        getAtoms(listAtom);

        return CalcGeom.getCentroid( listAtom.toArray(new Atom [listAtom.size()]), 1 );

    }
    
}
