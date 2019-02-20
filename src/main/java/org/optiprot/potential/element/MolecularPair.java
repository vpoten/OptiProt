/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import org.biojava.bio.structure.*;
import org.optiprot.maths.CalcGeom;

/**
 *
 * @author victor
 */
public class MolecularPair {

    private Atom m_atomA=null;
    private Atom m_atomB=null;
//    private String m_atomTypeA="";
//    private String m_atomTypeB = "";
    private int m_atomTypeA =-1;
    private int m_atomTypeB =-1;
    private double m_chargeA=0;
    private double m_chargeB=0;
    private double m_sqdistance=0;
    private double m_distance=0;
    private double m_bornradiiAB=0;


    public MolecularPair() {
        m_atomA=new AtomImpl();
        m_atomB=new AtomImpl();
    }

    public MolecularPair( Atom atA, Atom atB ) {

        m_atomA=atA;
        m_atomB=atB;
        
        calcDistance();
    }

    public double getBornRadiusAB() {
        return m_bornradiiAB;
    }

    public void setBornRadiusAB( double rad ) {
        m_bornradiiAB=rad;
    }

    public double getDistance() {
        return m_distance;
    }

    public double getSqrDistance() {
        return m_sqdistance;
    }

    public void setAtoms( double [] coordsA,  double [] coordsB ){
        m_atomA.setX(coordsA[0]);
        m_atomA.setY(coordsA[1]);
        m_atomA.setZ(coordsA[2]);

        m_atomB.setX(coordsB[0]);
        m_atomB.setY(coordsB[1]);
        m_atomB.setZ(coordsB[2]);

        calcDistance();
    }

    public void setAtomA(double [] coords){
        m_atomA.setX(coords[0]);
        m_atomA.setY(coords[1]);
        m_atomA.setZ(coords[2]);

        calcDistance();
    }

    public void setAtomB(double [] coords){
        m_atomB.setX(coords[0]);
        m_atomB.setY(coords[1]);
        m_atomB.setZ(coords[2]);

        calcDistance();
    }

    public void setAtomA(String name){
        m_atomA.setName(name);
    }

    public void setAtomB(String name){
        m_atomB.setName(name);
    }

    public String getNameAtomA(){
        return m_atomA.getName();
    }

    public String getNameAtomB(){
        return m_atomB.getName();
    }

//    /**
//     * @return the m_atomTypeA
//     */
//    public String getAtomTypeA() {
//        return m_atomTypeA;
//    }
//
//    /**
//     * @param m_atomTypeA the m_atomTypeA to set
//     */
//    public void setAtomTypeA(String m_atomTypeA) {
//        this.m_atomTypeA = m_atomTypeA;
//    }
//
//    /**
//     * @return the m_atomTypeB
//     */
//    public String getAtomTypeB() {
//        return m_atomTypeB;
//    }
//
//    /**
//     * @param m_atomTypeB the m_atomTypeB to set
//     */
//    public void setAtomTypeB(String m_atomTypeB) {
//        this.m_atomTypeB = m_atomTypeB;
//    }

    /**
     * @return the m_atomTypeA
     */
    public int getAtomTypeA() {
        return m_atomTypeA;
    }

    /**
     * @param m_atomTypeA the m_atomTypeA to set
     */
    public void setAtomTypeA(int m_atomTypeA) {
        this.m_atomTypeA = m_atomTypeA;
    }

    /**
     * @return the m_atomTypeB
     */
    public int getAtomTypeB() {
        return m_atomTypeB;
    }

    /**
     * @param m_atomTypeB the m_atomTypeB to set
     */
    public void setAtomTypeB(int m_atomTypeB) {
        this.m_atomTypeB = m_atomTypeB;
    }

    /**
     * @return the m_chargeA
     */
    public double getChargeA() {
        return m_chargeA;
    }

    /**
     * @param m_chargeA the m_chargeA to set
     */
    public void setChargeA(double m_chargeA) {
        this.m_chargeA = m_chargeA;
    }

    /**
     * @return the m_chargeB
     */
    public double getChargeB() {
        return m_chargeB;
    }

    /**
     * @param m_chargeB the m_chargeB to set
     */
    public void setChargeB(double m_chargeB) {
        this.m_chargeB = m_chargeB;
    }

    private void calcDistance(){
        m_sqdistance = CalcGeom.squareDistance(m_atomA, m_atomB);
        m_distance=Math.sqrt(m_sqdistance);
    }
}
