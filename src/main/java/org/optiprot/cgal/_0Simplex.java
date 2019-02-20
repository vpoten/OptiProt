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
public class _0Simplex extends BaseSimplex {


    public _0Simplex(double x, double y, double z) {

        super();
        
        vertices[0]=new AtomImpl();
        vertices[0].setCoords( new double [] {x,y,z} );
    }


    public int getNumVertices() {
        return 1;
    }

    
}
