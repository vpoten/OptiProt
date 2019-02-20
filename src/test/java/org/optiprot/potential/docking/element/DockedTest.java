/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import org.optiprot.potential.docking.Docked;
import java.util.Collection;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.CalcGeom;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class DockedTest {

    public DockedTest() {
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
     * Test of getProteinAtoms method, of class Docked.
     */
    @Test
    public void testGetProteinAtoms() {
        System.out.println("getProteinAtoms");
        
        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( "1a07", new int [] {0,1}, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( "1a07", false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid());


        Collection<Atom> result = instance.getProteinAtoms();
        int size=result.size();
        
        assertTrue( !result.isEmpty() );
       
    }

   

    /**
     * Test of setTransform method, of class Docked.
     */
    @Test
    public void testSetTransform() {
        System.out.println("setTransform");

        Atom axis = new AtomImpl();
        axis.setCoords(new double [] {1,0,0});
        double angle = 0.0;

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( "1a07", new int [] {0,1}, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( "1a07", false );

        Atom actSiteCoM=ligand.getCoM();

        //move ligand to origin
        Atom shift=CalcGeom.product(actSiteCoM, -1);

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid());
        instance.setTransform(axis, angle, shift);

        Collection<Atom> result = instance.getLigandAtoms();

        Atom com=Calc.getCentroid( result.toArray(new Atom [result.size()]) );

        assertTrue( Calc.amount(com)<1e-9 );

    }


}