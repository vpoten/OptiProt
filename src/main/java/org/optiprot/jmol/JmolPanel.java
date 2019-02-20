/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jmol;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.gui.BiojavaJmol;

/**
 *
 * @author victor
 */
public class JmolPanel {

    private Structure struct=null;
    BiojavaJmol jmolPanel=null;

    public JmolPanel() {
        jmolPanel = new BiojavaJmol();
        struct=new StructureImpl();
    }

  
    /**
     * @param struc the struc to set
     */
    public void setStruct() {
        jmolPanel.setStructure(struct);
    }

    public void evalCommand( String command ){
        jmolPanel.evalString(command);
    }

    public void addChain( Chain chain ){
        struct.addChain(chain);
    }

    public void addHetatm( int chainidx, HetatomImpl het ){
        struct.getChain(chainidx).addGroup(het);
    }

}
