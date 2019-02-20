/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io.mol2;

/**
 *
 * @author victor
 */
public class Mol2SubStructure {
    private String name="";
    private Integer rootAtom = null;
    private String type = "";
    private String chain = "";
    private String subType = "";
    private Integer intBonds = null;

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
     * @return the rootAtom
     */
    public Integer getRootAtom() {
        return rootAtom;
    }

    /**
     * @param rootAtom the rootAtom to set
     */
    public void setRootAtom(Integer rootAtom) {
        this.rootAtom = rootAtom;
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
     * @return the chain
     */
    public String getChain() {
        return chain;
    }

    /**
     * @param chain the chain to set
     */
    public void setChain(String chain) {
        this.chain = chain;
    }

    /**
     * @return the subType
     */
    public String getSubType() {
        return subType;
    }

    /**
     * @param subType the subType to set
     */
    public void setSubType(String subType) {
        this.subType = subType;
    }

    /**
     * @return the intBonds
     */
    public Integer getIntBonds() {
        return intBonds;
    }

    /**
     * @param intBonds the intBonds to set
     */
    public void setIntBonds(Integer intBonds) {
        this.intBonds = intBonds;
    }



}
