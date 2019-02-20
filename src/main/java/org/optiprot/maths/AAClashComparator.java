/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.Comparator;

/**
 *
 * @author victor
 */
public class AAClashComparator implements Comparator<AAClash> {

    /**
     * Returns a negative integer, zero, or a positive integer as the first
     * argument is less than, equal to, or greater than the second
     * @param o1
     * @param o2
     * @return
     */
    public int compare(AAClash o1, AAClash o2) {

        int dist1=o1.distance();
        int dist2=o2.distance();

        //sort in reverse order
        if( dist1>dist2 )
            return -1;

        if( dist1<dist2 )
            return 1;

        return 0;
    }

}
