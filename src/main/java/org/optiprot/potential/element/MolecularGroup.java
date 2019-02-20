/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import java.util.ArrayList;
import java.util.HashMap;
import org.biojava.bio.structure.*;

/**
 *
 * @author victor
 */
public class MolecularGroup {

    private ArrayList<Atom> m_elements=new ArrayList<Atom>();
    private HashMap<String,Double> m_charges=new HashMap<String,Double>();
    //private HashMap<String,String> m_atomTypes=new HashMap<String,String>();
    private HashMap<String,Integer> m_atomTypes=new HashMap<String,Integer>();
    

    /**
     * @return the m_elements
     */
    public ArrayList<Atom> getElements() {
        return m_elements;
    }

    public int size() {
        return getElements().size();
    }

    /**
     * @param m_elements the m_elements to set
     */
    private void setElements(ArrayList<Atom> m_elements) {
        this.m_elements = m_elements;
    }

    /**
     * @return the m_charges
     */
    protected HashMap<String, Double> getCharges() {
        return m_charges;
    }

    /**
     * @param m_charges the m_charges to set
     */
    private void setCharges(HashMap<String, Double> m_charges) {
        this.m_charges = m_charges;
    }

    public double getCharge( Atom atom ){

        return this.getCharges().get(atom.getName());
    }

//    /**
//     *
//     * @param at
//     * @param charge
//     * @param type
//     */
//    public void addAtom( Atom at, double charge, String type ){
//        this.getElements().add(at);
//        this.getCharges().put(at.getName(), charge);
//        this.getAtomTypes().put(at.getName(), type);
//    }

    /**
     *
     * @param at
     * @param charge
     * @param type
     */
    public void addAtom( Atom at, double charge, int type ){
        this.getElements().add(at);
        this.getCharges().put(at.getName(), charge);
        this.getAtomTypes().put(at.getName(), type);
    }

//    /**
//     * @return the m_atomTypes
//     */
//    protected HashMap<String, String> getAtomTypes() {
//        return m_atomTypes;
//    }
//
//    /**
//     * @param m_atomTypes the m_atomTypes to set
//     */
//    private void setAtomTypes(HashMap<String, String> m_atomTypes) {
//        this.m_atomTypes = m_atomTypes;
//    }

    /**
     * @return the m_atomTypes
     */
    protected HashMap<String, Integer> getAtomTypes() {
        return m_atomTypes;
    }

    /**
     * @param m_atomTypes the m_atomTypes to set
     */
    private void setAtomTypes(HashMap<String, Integer> m_atomTypes) {
        this.m_atomTypes = m_atomTypes;
    }

//    public String getAtomType( Atom atom ){
//
//        return this.getAtomTypes().get(atom.getName());
//    }

    public Integer getAtomType( Atom atom ){

        return this.getAtomTypes().get(atom.getName());
    }
    
}
