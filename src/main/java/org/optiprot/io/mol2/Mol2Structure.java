/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io.mol2;

import java.util.ArrayList;
import java.util.HashMap;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.StandardAminoAcid;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.io.PDBParseException;

/**
 *
 * @author victor
 */
public class Mol2Structure extends StructureImpl {
    private String comment="";
    private String chargeType = "";
    private String molType="";

    private ArrayList<Mol2Atom> atoms=new ArrayList<Mol2Atom>();
    private ArrayList<Mol2Bond> bonds=new ArrayList<Mol2Bond>();
    private ArrayList<Mol2SubStructure> subStructs=new ArrayList<Mol2SubStructure>();


    static private final String AMINOACID = "RESIDUE";
    static private final String HETEROATOM = "GROUP";

    /**
     * @return the comment
     */
    public String getComment() {
        return comment;
    }

    /**
     * @param comment the comment to set
     */
    public void setComment(String comment) {
        this.comment = comment;
    }

    /**
     * @return the chargeType
     */
    public String getChargeType() {
        return chargeType;
    }

    /**
     * @param chargeType the chargeType to set
     */
    public void setChargeType(String chargeType) {
        this.chargeType = chargeType;
    }

    /**
     * @return the molType
     */
    public String getMolType() {
        return molType;
    }

    /**
     * @param molType the molType to set
     */
    public void setMolType(String molType) {
        this.molType = molType;
    }

    /**
     * @return the atoms
     */
    public ArrayList<Mol2Atom> getAtoms() {
        return atoms;
    }

    /**
     * @return the bonds
     */
    public ArrayList<Mol2Bond> getBonds() {
        return bonds;
    }

    /**
     * @return the subStructs
     */
    public ArrayList<Mol2SubStructure> getSubStructs() {
        return subStructs;
    }

    
    public void arrange(){

        HashMap<String,Chain> chainTable=new HashMap<String,Chain>();
        ArrayList<Group> groups=new ArrayList<Group>();

        //create chains
        for( Mol2SubStructure substr : getSubStructs() ){

            if( !chainTable.containsKey( substr.getChain() ) ){
                Chain ch=new ChainImpl();
                ch.setName( substr.getChain() );
                chainTable.put( substr.getChain(), ch );

                this.addChain(ch);
            }

            Group grp=null;

            //add structure to chain
            if( substr.getType().equals(AMINOACID) ){
               
                Group grp2=StandardAminoAcid.getAminoAcid( substr.getSubType() );

                if( grp2!=null ){
                    grp=(Group) grp2.clone();
                    grp.clearAtoms();
                }
                else{
                    grp=new HetatomImpl();
                    grp.setPDBCode("");
                    
                    try {
                        grp.setPDBName( substr.getName() );
                    } catch (PDBParseException ex) {}

                }
            }
            else{
                grp=new HetatomImpl();
                grp.setPDBCode("");

                try {
                    grp.setPDBName( substr.getName() );
                } catch (PDBParseException ex) {}
            }

            chainTable.get(substr.getChain() ).addGroup( grp );
            groups.add(grp);

        }

        if( getSubStructs().isEmpty() ){
            //if theres not substructures create one and add atoms

            Chain ch=new ChainImpl();
            ch.setName("A");
            this.addChain(ch);

            Group grp=new HetatomImpl();
            grp.setPDBCode("");
            ch.addGroup(grp);
            groups.add(grp);

            for( Mol2Atom atom: getAtoms() ){
                grp.addAtom(atom);
            }
        }
        else{
            //fill groups
            for( Mol2Atom atom: getAtoms() ){

                if( atom.getSubstr()<=0 )
                    continue;//if the atom dont has substructure

                Group grp=groups.get( atom.getSubstr()-1 );
                grp.addAtom(atom);
            }
        }

        //fill bonds
        for( Mol2Bond bond: getBonds() ){
            bond.setAtOrigin( getAtoms().get(bond.getIdOrigin()-1) );
            bond.setAtTarget( getAtoms().get(bond.getIdTarget()-1) );
            bond.getAtOrigin().getBonds().add(bond);
            bond.getAtTarget().getBonds().add(bond);
        }
    }
    
    
}
