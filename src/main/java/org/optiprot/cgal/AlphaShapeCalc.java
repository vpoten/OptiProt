/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.cgal;

import org.biojava.bio.structure.Atom;

/**
 *
 * @author victor
 */

class AlphaShapeCalc {

    private long ptr2cgalAlphaShape=0;

    private native long doTriangulation( Atom [] atoms, double [] radii, double addRadii );
    public native BaseSimplex [] getBoundary();
    public native double calcArea(int depth);
    public native double [] calcAreaVol(int depth, int iterMont );


    static{      
        String libname=System.mapLibraryName("cgal_native");
        System.load("/usr/local/lib/"+libname);
    }


    public AlphaShapeCalc() {
    }

    public void triangulate( Atom [] atoms, double [] radii, double addRadii ){

        setPtr2cgalAlphaShape(doTriangulation(atoms, radii, addRadii));
    }

    /**
     * @return the ptr2cgalAlphaShape
     */
    private long getPtr2cgalAlphaShape() {
        return ptr2cgalAlphaShape;
    }

    /**
     * @param ptr2cgalAlphaShape the ptr2cgalAlphaShape to set
     */
    private void setPtr2cgalAlphaShape(long ptr2cgalAlphaShape) {
        this.ptr2cgalAlphaShape = ptr2cgalAlphaShape;
    }

    
    static public BaseSimplex createSimplex(double x,double y,double z){

        return new _0Simplex(x,y,z);
    }

    static public BaseSimplex createSimplex(double x,double y,double z,
            double x2,double y2,double z2){

        return new _1Simplex(x,y,z,x2,y2,z2);
    }

    static public BaseSimplex createSimplex(double x,double y,double z,
            double x2,double y2,double z2,
            double x3,double y3,double z3 ){

        return new _2Simplex(x,y,z,x2,y2,z2,x3,y3,z3);
    }
}
