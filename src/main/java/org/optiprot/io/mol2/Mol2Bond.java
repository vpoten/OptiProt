/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io.mol2;

/**
 *
 * @author victor
 */
public class Mol2Bond {
    private Integer idOrigin=null;
    private Integer idTarget = null;
    private String type = "";
    private Mol2Atom atOrigin = null;
    private Mol2Atom atTarget = null;

    /**
     * @return the idOrigin
     */
    public Integer getIdOrigin() {
        return idOrigin;
    }

    /**
     * @param idOrigin the idOrigin to set
     */
    public void setIdOrigin(Integer idOrigin) {
        this.idOrigin = idOrigin;
    }

    /**
     * @return the idTarget
     */
    public Integer getIdTarget() {
        return idTarget;
    }

    /**
     * @param idTarget the idTarget to set
     */
    public void setIdTarget(Integer idTarget) {
        this.idTarget = idTarget;
    }

    /**
     * @return the type
     */
    public String getType() {
        return type;
    }

    /**
     * @param type the type to set
     */
    public void setType(String type) {
        this.type = type;
    }

    /**
     * @return the atOrigin
     */
    public Mol2Atom getAtOrigin() {
        return atOrigin;
    }

    /**
     * @param atOrigin the atOrigin to set
     */
    public void setAtOrigin(Mol2Atom atOrigin) {
        this.atOrigin = atOrigin;
    }

    /**
     * @return the atTarget
     */
    public Mol2Atom getAtTarget() {
        return atTarget;
    }

    /**
     * @param atTarget the atTarget to set
     */
    public void setAtTarget(Mol2Atom atTarget) {
        this.atTarget = atTarget;
    }

    /**
     * get the neighbor of atom in this bond
     * @param at
     * @return
     */
    public Mol2Atom getNeighbor( Mol2Atom at ) {

        if( getAtOrigin()==at )
            return getAtTarget();

        if( getAtTarget()==at )
            return getAtOrigin();

        return null;
    }

    public boolean isSingle(){
        return getType().equals("1");
    }

    public boolean isDouble(){
        return getType().equals("2");
    }

    public boolean isTriple(){
        return getType().equals("3");
    }

    public boolean isAromatic(){
        return getType().equals("ar");
    }

    public boolean isAmide(){
        return getType().equals("am");
    }

}
