/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import java.util.ArrayList;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.PDBParseException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.jmol.JmolPanel;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class DockingProteinTest {

    public DockingProteinTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of calcActivePoints method, of class DockingProtein.
     */
    @Test
    public void testCalcActivePoints() throws StructureException, PDBParseException {
        System.out.println("calcActivePoints");

        String pdbcode="1acm";

        DockingProtein instance=ConfigOptiProtTest.readAstexProtein( pdbcode, new int [] {0,1}, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom center =ligand.getCoM();
        
        instance.calcActivePoints(center, ligand.calcRadii(), 
                ConfigOptiProtTest.getParameters().getGrid(),true);

        

        DockingProtein protein2=ConfigOptiProtTest.readAstexProtein( pdbcode, new int [] {0,1}, true);

        ArrayList<Atom> list=(ArrayList<Atom>) instance.getActSite();
        Atom centroid=Calc.getCentroid( list.toArray(new Atom[list.size()]));

        drawProtLigand( protein2.getChain() , list );

        assertTrue( Calc.getDistance(center, centroid)<4.0 );
    }

    
    static public void drawProtLigand( Chain chain , List<Atom> ligand )
            throws PDBParseException {

        //draw the chain and the ligand (jmol)
        JmolPanel panel=new JmolPanel();
        panel.addChain(chain);

        HetatomImpl het=new HetatomImpl();

        for( Atom at : ligand ){
            het.addAtom(at);
        }
        ///het.addAtom( center );
        het.setPDBCode("");
        het.setPDBName("");
        panel.addHetatm(0, het);

        panel.setStruct();

        panel.evalCommand("select * ; colour chain;");
        panel.evalCommand("select hetero ; colour red;");
        panel.evalCommand("select amino ; isosurface ignore(hetero) molecular");
        ///panel.evalCommand("select *; spacefill off; wireframe off; backbone 0.4;");

    }
    

}