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
public abstract class BaseSimplex {

    protected Atom [] vertices=null;
    private double [] weights=null;
    private double outAngle=-1.0;
    private Atom normal=null;

    abstract public int getNumVertices();

    public BaseSimplex() {
        weights=new double [getNumVertices()];
        vertices=new Atom [getNumVertices()];
    }


    public Atom[] getVertices() {
        return vertices;
    }

    /**
     * @return the weights
     */
    public double[] getWeights() {
        return weights;
    }

    /**
     * @return the outAngle
     */
    public double getOutAngle() {
        return outAngle;
    }

    /**
     * @param outAngle the outAngle to set
     */
    public void setOutAngle(double outAngle) {
        this.outAngle = outAngle;
    }

    /**
     * @param weights the weights to set
     */
    public void setWeights(double[] weights) {
        this.weights=weights;
    }

    /**
     * @return the normal
     */
    public Atom getNormal() {
        return normal;
    }

    /**
     * @param normal the normal to set
     */
    public void setNormal(Atom normal) {
        this.normal = normal;
    }

    public void setNormal(double x,double y, double z) {
        this.normal = new AtomImpl();
        this.normal.setCoords(new double [] {x,y,z});
    }

    public void printSimplex(){
        if( this.getNumVertices()==1 ){
            System.out.println("0Simplex");
        }
        else if( this.getNumVertices()==2 ){
            System.out.println("1Simplex");
        }
        else if( this.getNumVertices()==3 ){
            System.out.println("2Simplex");
        }

        System.out.println("Out Angle: "+this.getOutAngle());

        for( Atom a : this.getVertices() ){
            System.out.println(a.getX()+","+a.getY()+","+a.getZ());
        }
    }
}
