/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking;

import java.util.ArrayList;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.jama.Matrix;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;

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
     * Test of getLigandAtoms method, of class Docked.
     */
    @Test
    public void testGetLigandAtoms() {
        System.out.println("getLigandAtoms");

        String pdbcode="1a0q";
        int numAtoms=43;

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        assertTrue( instance.getLigandAtoms().size()==numAtoms );
    }

    /**
     * Test of getLigandAtomsCopy1 method, of class Docked.
     */
    @Test
    public void testGetLigandAtomsCopy1() {
        System.out.println("getLigandAtomsCopy1");

        String pdbcode="1a1b";
        int numAtoms=78;

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        assertTrue( instance.getLigandAtomsCopy1().size()==numAtoms );
    }


    /**
     * Test of setLigand method, of class Docked.
     */
    @Test
    public void testSetLigand() {
        System.out.println("setLigand");

        String pdbcode="1a1e";
        int numActAtoms=14;

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        assertTrue( instance.getNumActAtoms()==numActAtoms );

        ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, true );

        instance.setLigand(ligand);

        assertTrue( instance.getNumActAtoms()==numActAtoms );
    }
   

    /**
     * Test of setTransform method, of class Docked.
     */
    @Test
    public void testSetTransform_Matrix() {
        System.out.println("setTransform");
        Matrix mat = Matrix.identity(4, 4);

        mat.set(0, 0, 4);
        mat.set(1, 0, 3);
        mat.set(2, 0, 2);

        String pdbcode="1a1e";
        
        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        instance.setTransform(mat);

        assertTrue( instance.getMatTransf().get(0,0)==mat.get(0, 0) );
        assertTrue( instance.getMatTransf().get(1,0)==mat.get(1, 0) );
        assertTrue( instance.getMatTransf().get(2,0)==mat.get(2, 0) );
    }


    /**
     * Test of getNumActAtoms method, of class Docked.
     */
    @Test
    public void testGetNumActAtoms() {
        System.out.println("getNumActAtoms");

        String pdbcode="1a1e";
        int numActAtoms=14;

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        assertTrue( instance.getNumActAtoms()==numActAtoms );
    }

    /**
     * Test of getCompActAtom method, of class Docked.
     */
    @Test
    public void testGetCompActAtom() {
        System.out.println("getCompActAtom");

        String pdbcode="1acm";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        int id_point = (int) Math.floor( Math.random()*instance.getNumActPoints() );

        Integer id_atom = instance.getCompActAtom(id_point);

        assertTrue( id_atom!=null );
    }

    /**
     * Test of getCompActPoint method, of class Docked.
     */
    @Test
    public void testGetCompActPoint() {
        System.out.println("getCompActPoint");

        String pdbcode="1acm";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        int id_atom = (int) Math.floor( Math.random()*instance.getNumActAtoms() );

        Integer id_point = instance.getCompActPoint(id_atom);

        assertTrue( id_point!=null );
    }

    
    /**
     * Test of calcOverlapFactor method, of class Docked.
     */
    @Test
    public void testCalcOverlapFactor() {
        System.out.println("calcOverlapFactor");

        String pdbcode="1aco";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        double expResult = 1.0;
        double result = instance.calcOverlapFactor();

        assertTrue( Math.abs(result-expResult)<1e-6 );
    }

   

    /**
     * Test of isModePolar method, of class Docked.
     */
    @Test
    public void testIsModePolar() {
        System.out.println("isModePolar");

        String pdbcode="1fen";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        assertEquals( instance.isModePolar(), false );

        pdbcode="1acm";

        protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        actSiteCoM=ligand.getCoM();

        instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        assertEquals( instance.isModePolar(), true );
    }


    /**
     * Test of getLigandAtomsInActivePoints method, of class Docked.
     */
    @Test
    public void testGetLigandAtomsInActivePoints() {
        System.out.println("getLigandAtomsInActivePoints");

        String pdbcode="1lmo";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        List<Atom> list = new ArrayList<Atom>();
        double dist = 0.7;

        Docked instance = new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        instance.getLigandAtomsInActivePoints(list, dist);

        assertTrue( list.size()>0 );
    }

}