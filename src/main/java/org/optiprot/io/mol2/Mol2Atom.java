/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io.mol2;

import java.util.ArrayList;
import org.biojava.bio.structure.AtomImpl;

/**
 *
 * @author victor
 */
public class Mol2Atom extends AtomImpl {
    private String type="";
    private int atomicNum=0;
    private Integer substr = null;
    private Float charge = null;
    private ArrayList<Mol2Bond> bonds=new ArrayList<Mol2Bond>();

    
    public int getAtomicNum(){
        return atomicNum;
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

        setAtomicNum(0);

        if( type.length()<=2 ){
            if( Elements.getElement(type)!=null )
                setAtomicNum( Elements.getElement(type).getNum() );
        }
        else if( type.indexOf('.')>=0 ){

            ElementEntry ent =
                    Elements.getElement(type.substring(0, type.indexOf('.')));

            if(ent!=null)
                setAtomicNum( ent.getNum() );
        }
        
    }

    /**
     * @return the substr
     */
    public Integer getSubstr() {
        return substr;
    }

    /**
     * @param substr the substr to set
     */
    public void setSubstr(Integer substr) {
        this.substr = substr;
    }

    /**
     * @return the charge
     */
    public Float getCharge() {
        return charge;
    }

    /**
     * @param charge the charge to set
     */
    public void setCharge(Float charge) {
        this.charge = charge;
    }

    /**
     * @param atomicNum the atomicNum to set
     */
    protected void setAtomicNum(int atomicNum) {
        this.atomicNum = atomicNum;
    }

    /**
     * @return the bonds
     */
    public ArrayList<Mol2Bond> getBonds() {
        return bonds;
    }

    /**
     * @param bonds the bonds to set
     */
    public void setBonds(ArrayList<Mol2Bond> bonds) {
        this.bonds = bonds;
    }

    public boolean isCarbon() {
        return (getAtomicNum()==6);
    }

    public boolean isOxygen() {
        return (getAtomicNum()==8);
    }

    public boolean isNitrogen() {
        return (getAtomicNum()==7);
    }

    public boolean isSulfur() {
        return (getAtomicNum()==16);
    }

    public boolean isFluorine() {
        return (getAtomicNum()==9);
    }

    public boolean isHydrogen() {
        return (getAtomicNum()==1);
    }

    public boolean isHbondAcceptor() {

        if( isOxygen() || isFluorine() )
            return true;//O, F

        if( isNitrogen() ){
            // N+ ions and sp2 hybrid N with 3 valences should not be Hbond acceptors
            if( !((getValence() == 4 && getHyb() == 3) ||
                    (getValence() == 3 && getHyb() == 2)) )
              return true;
        }

        return false;
    }

    public boolean isHbondDonor() {

        if( !(isOxygen() || isFluorine() || isNitrogen()) )
            return false;//O, F, N

        //explore neighbors atoms
        for( Mol2Bond bond : getBonds() ){
            if( bond.getNeighbor(this).isHydrogen() )
                return true;
        }

        return false;
    }

    public boolean isHbondDonorH() {

        if( !isHydrogen() )
            return false;

        //explore neighbors atoms
        for( Mol2Bond bond : getBonds() ){
            if( bond.getNeighbor(this).isHbondDonor() )
                return true;
        }

        return false;
    }


    /**
     * The current number of explicit connections
     *
     * @return
     */
    public int getValence() {
        return getBonds().size();
    }

    /**
     * The hybridization of this atom (i.e. 1 for sp, 2 for sp2, 3 for sp3)
     * used only for C, N & O
     *
     * @return
     */
    public int getHyb() {
        
        if( isCarbon() ){
            if( getType().equals("C.3"))
                return 3;
            if( getType().equals("C.2"))
                return 2;
            if( getType().equals("C.1"))
                return 1;
            if( getType().equals("C.ar"))
                return 2;
            return 3;
        }
        else if( isNitrogen() ){
            if( getType().equals("N.3"))
                return 3;
            if( getType().equals("N.2"))
                return 2;
            if( getType().equals("N.1"))
                return 1;
            if( getType().equals("N.ar"))
                return 2;
            if( getType().equals("N.am"))
                return 3;
            if( getType().equals("N.pl3"))
                return 3;
            if( getType().equals("N.4"))
                return 3;
            return 3;
        }
        else if( isOxygen() ){
            if( getType().equals("O.3"))
                return 3;
            if( getType().equals("O.2"))
                return 2;
            if( getType().equals("O.co2"))
                return 2;
            return 3;
        }

        return 0;
    }

    /**
     * Is this atom an element in the 15th or 16th main groups
     * (i.e., N, O, P, S ...) ?
     *
     * @return
     */
    boolean isHeteroatom(){

        switch( getAtomicNum() ){
          case 7: return true;
          case 8: return true;
          case 15: return true;
          case 16: return true;
          case 33: return true;
          case 34: return true;
          case 51: return true;
          case 52: return true;
          case 83: return true;
          case 84: return true;
        }

       return false;
    }

}
