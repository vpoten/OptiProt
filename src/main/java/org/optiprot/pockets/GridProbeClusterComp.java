/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.pockets;

import java.util.Comparator;

/**
 *
 * @author victor
 */
public class GridProbeClusterComp implements Comparator<GridProbeCluster> {

    /**
     * Returns a negative integer, zero, or a positive integer as the first
     * argument is less than, equal to, or greater than the second
     * @param o1
     * @param o2
     * @return
     */
    public int compare(GridProbeCluster o1, GridProbeCluster o2) {

        if( o1.getProbes().size() < o2.getProbes().size() )
            return -1;
        
        if( o1.getProbes().size() > o2.getProbes().size() )
            return 1;

        return 0;
    }

}
