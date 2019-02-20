/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import java.util.ArrayList;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;

/**
 *
 * @author victor
 */
public class DockingLigand extends DockingMolecule {

    private ArrayList<Atom> atoms=new ArrayList<Atom>();

    public DockingLigand( Chain chain ) {
        super(chain);
    }
    
    /**
     *
     * @return
     */
    public List<Atom> getAtoms(){

        if( getListAtoms().isEmpty() ){

            for(Group group : this.getChain().getAtomGroups() ){
                for( Atom at : group.getAtoms() ){
                    getListAtoms().add(at);
                }
            }
        }

        return getListAtoms();
    }


    /**
     * @return the listAtoms
     */
    private ArrayList<Atom> getListAtoms() {
        return atoms;
    }

}
