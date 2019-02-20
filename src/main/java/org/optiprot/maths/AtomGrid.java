/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;

/**
 *
 * @author victor
 */
public class AtomGrid {

    private Atom [] grid=null;
    private boolean [] inVolume=null;
    private double [] values=null;
    private int npointx=0;
    private int npointy=0;
    private int npointz=0;
    private int size=0;
    private double resol=0;

    public AtomGrid( int nX, int nY, int nZ ) {
        grid=new Atom [nX*nY*nZ];
        inVolume=new boolean [nX*nY*nZ];
        values=new double [nX*nY*nZ];

        npointx=nX;
        npointy=nY;
        npointz=nZ;
        size=nX*nY*nZ;

        for(int i=0;i<size();i++){
            grid[i]=new AtomImpl();
            inVolume[i]=false;
            values[i]=0.0;
        }
    }

    public void setDimension( int nX, int nY, int nZ ) {

        if( nX*nY*nZ>size() )
            throw new RuntimeException("Grid dimension exceeded");

        npointx=nX;
        npointy=nY;
        npointz=nZ;
    }

    public int getIndex(int i, int j, int k){
        return i*getNpointy()*getNpointz()+j*getNpointz()+k;
    }

    /**
     *
     * @param i
     * @param j
     * @param k
     * @return : the index or -1 if i,j,k are out of bounds
     */
    public int getIndexSec(int i, int j, int k){

        if( i<0 || j<0 || k<0 ){
            return -1;
        }
        if( i>=getNpointx() || j>=getNpointy() || k>=getNpointz() ){
            return -1;
        }

        return i*getNpointy()*getNpointz()+j*getNpointz()+k;
    }

    public int getX( int index ){
        return index/(getNpointy()*getNpointz());
    }

    public int getY( int index ){
        return (index%(getNpointy()*getNpointz()))/getNpointz();
    }

    public int getZ( int index ){
        return (index%(getNpointy()*getNpointz()))%getNpointz();
    }

    public void setPoint( int idx, Atom at, boolean inVol ){
        grid[idx].setX( at.getX() );
        grid[idx].setY( at.getY() );
        grid[idx].setZ( at.getZ() );
        inVolume[idx]=inVol;
    }

    public Atom getPoint( int idx ){
        return grid[idx];
    }

    public boolean isInVolume( int idx ){
        return inVolume[idx];
    }

    public void setInVolume( int idx, boolean val ){
        inVolume[idx]=val;
    }

    public double getValue( int idx ){
        return values[idx];
    }

    public void setValue( int idx, double val ){
        values[idx]=val;
    }

    public int size(){
        return size;
    }

    /**
     * @return the npointx
     */
    public int getNpointx() {
        return npointx;
    }

    /**
     * @return the npointy
     */
    public int getNpointy() {
        return npointy;
    }

    /**
     * @return the npointz
     */
    public int getNpointz() {
        return npointz;
    }

    /**
     * clear the elements of the grid
     */
    public void clear(){

        for(int i=0;i<size();i++){
            //grid[i].setX(0);
            //grid[i].setY(0);
            //grid[i].setZ(0);
            inVolume[i]=false;
            values[i]=0.0;
        }
    }

    /**
     * @return the resol
     */
    public double getResol() {
        return resol;
    }

    /**
     * @param resol the resol to set
     */
    public void setResol(double resol) {
        this.resol = resol;
    }

}
