/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import java.util.ArrayList;
import java.util.HashMap;
import org.biojava.bio.structure.Atom;

/**
 * contains a residue of CHARMM topology file
 *
 * @author victor
 */
public class CharmmResidue {

    private String name="";
    private String description="";
    private double charge=0;

    private ArrayList<MolecularBond> m_bonds=new ArrayList<MolecularBond>();
    private ArrayList<MolecularGroup> m_groups=new ArrayList<MolecularGroup>();
    private ArrayList<MolecularAngle> m_angles=new ArrayList<MolecularAngle>();
    private ArrayList<MolecularImproper> m_impropers=new ArrayList<MolecularImproper>();
    private ArrayList<MolecularDihedral> m_dihedrals=new ArrayList<MolecularDihedral>();

    //atomtypes table
    private HashMap<String,Integer> m_atomTypes=new HashMap<String,Integer>();

    //charges table
    private HashMap<String,Double> m_atomCharges=new HashMap<String,Double>();

    public CharmmResidue() {
        
    }

//    /**
//     * search in the groups list for the atom type
//     * @param name
//     * @return
//     */
//    public String getAtomType(String name) {
//
//        String type=null;
//
//        if( name.startsWith("+") || name.startsWith("-") )
//            name=name.substring(1);
//
//        for( MolecularGroup group : this.getGroups() ){
//
//            type=group.getAtomTypes().get(name);
//
//            if( type!=null )
//                break;
//        }
//
//        return type;
//    }

    /**
     * search in the atom types table for the atom type
     * @param name
     * @return
     */
    public Integer getAtomType(String name) {

        if( name.startsWith("+") || name.startsWith("-") )
            name=name.substring(1);

        return this.getAtomTypes().get(name);
    }

    /**
     * search in the charges table for the atom charge
     * @param name
     * @return
     */
    public Double getAtomCharge(String name) {

        if( name.startsWith("+") || name.startsWith("-") )
            name=name.substring(1);

        return this.getAtomCharges().get(name);
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the charge
     */
    public double getCharge() {
        return charge;
    }

    /**
     * @param charge the charge to set
     */
    public void setCharge(double charge) {
        this.charge = charge;
    }

    /**
     * @return the m_bonds
     */
    public ArrayList<MolecularBond> getBonds() {
        return m_bonds;
    }

    /**
     * @param m_bonds the m_bonds to set
     */
    private void setBonds(ArrayList<MolecularBond> bonds) {
        this.m_bonds = bonds;
    }

    /**
     * @return the m_groups
     */
    private ArrayList<MolecularGroup> getGroups() {
        return m_groups;
    }

    /**
     * adds group to residue and update the tables of charges and atom types
     * 
     * @param group
     */
    public void addGroup(MolecularGroup group) {

        this.getGroups().add(group);

        for( Atom atom : group.getElements() ){
            this.getAtomTypes().put(atom.getName(), group.getAtomType(atom));
            this.getAtomCharges().put(atom.getName(), group.getCharge(atom));
        }
    }

    /**
     * @param m_groups the m_groups to set
     */
    private void setGroups(ArrayList<MolecularGroup> groups) {
        this.m_groups = groups;
    }

    /**
     * @return the m_angles
     */
    public ArrayList<MolecularAngle> getAngles() {
        return m_angles;
    }

    /**
     * @param m_angles the m_angles to set
     */
    private void setAngles(ArrayList<MolecularAngle> m_angles) {
        this.m_angles = m_angles;
    }

    /**
     * @return the m_impropers
     */
    public ArrayList<MolecularImproper> getImpropers() {
        return m_impropers;
    }

    /**
     * @param m_impropers the m_impropers to set
     */
    private void setImpropers(ArrayList<MolecularImproper> m_impropers) {
        this.m_impropers = m_impropers;
    }

    /**
     * @return the m_dihedrals
     */
    public ArrayList<MolecularDihedral> getDihedrals() {
        return m_dihedrals;
    }

    /**
     * @param m_dihedrals the m_dihedrals to set
     */
    private void setDihedrals(ArrayList<MolecularDihedral> m_dihedrals) {
        this.m_dihedrals = m_dihedrals;
    }

    /**
     * @return the description
     */
    public String getDescription() {
        return description;
    }

    /**
     * @param description the description to set
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * @return the m_atomTypes
     */
    private HashMap<String, Integer> getAtomTypes() {
        return m_atomTypes;
    }

    /**
     * @return the m_atomCharges
     */
    private HashMap<String, Double> getAtomCharges() {
        return m_atomCharges;
    }

}
