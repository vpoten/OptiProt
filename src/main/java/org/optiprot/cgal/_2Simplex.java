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
public class _2Simplex extends BaseSimplex {

  
    _2Simplex(double x, double y, double z,
            double x2, double y2, double z2,
            double x3, double y3, double z3) {

        super();

        vertices[0]=new AtomImpl();
        vertices[0].setCoords( new double [] {x,y,z} );
        vertices[1]=new AtomImpl();
        vertices[1].setCoords( new double [] {x2,y2,z2} );
        vertices[2]=new AtomImpl();
        vertices[2].setCoords( new double [] {x3,y3,z3} );
    }


    public int getNumVertices() {
        return 3;
    }

   

}
