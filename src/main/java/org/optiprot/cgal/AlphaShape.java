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
public class AlphaShape {

    private BaseSimplex [] simplices=null;
    private AlphaShapeCalc as=null;
    private double area=0;
    private double volume=0;

    public AlphaShape( Atom [] molecule, double [] radii, double radiiAdd ) {

        as=new AlphaShapeCalc();

        //calculates the alpha-shape (with alpha=0)
        as.triangulate(molecule, radii, radiiAdd);
    }

    /**
     * @return the simplices
     */
    public BaseSimplex[] getSimplices() {

        if( simplices==null ){
            //get the boundary of alpha-complex
            this.setSimplices(as.getBoundary());
        }

        return simplices;
    }

    /**
     * @param simplices the simplices to set
     */
    public void setSimplices(BaseSimplex[] simplices) {
        this.simplices = simplices;
    }

    /**
     * calc the area
     * @param depth : max depth of recursivity (4 for accuracy)
     * @return
     */
    public double getArea(int depth){

        if( area==0 ){
            area=as.calcArea(depth);
        }

        return area;
    }

    /**
     * calc the volume
     * @param iterMont : iteration of montecarlo method (250000 for accuracy)
     * @return
     */
    public double getVolume(int iterMont){

        if( volume==0 ){
            double [] res = as.calcAreaVol(0, iterMont);
            volume=res[1];
        }

        return volume;
    }

}