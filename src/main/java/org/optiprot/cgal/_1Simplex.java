/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.cgal;

import org.biojava.bio.structure.*;

/**
 *
 * @author victor
 */
public class _1Simplex extends BaseSimplex  {

     
    _1Simplex(double x, double y, double z,
            double x2, double y2, double z2) {

        super();

        vertices[0]=new AtomImpl();
        vertices[0].setCoords( new double [] {x,y,z} );
        vertices[1]=new AtomImpl();
        vertices[1].setCoords( new double [] {x2,y2,z2} );
    }


    public int getNumVertices() {
        return 2;
    }

}
