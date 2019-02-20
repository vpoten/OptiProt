/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class DockingMoleculeTest {

    public DockingMoleculeTest() {
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
     * Test of getInternalEnergy method, of class DockingMolecule.
     */
    @Test
    public void testGetInternalEnergy() {
        System.out.println("getInternalEnergy");

        DockingLigand instance=ConfigOptiProtTest.readAstexLigand( "1a07", false );
        
        Double result = instance.getInternalEnergy();
        
        assertTrue( result!=null );

    }

    /**
     * Test of readCDKMol2 method, of class DockingMolecule.
     */
    @Test
    public void testReadCDKMol2() {
        System.out.println("readCDKMol2");

        DockingLigand instance=ConfigOptiProtTest.readAstexLigand( "1a07", false );

        instance.readCDKMol2( ConfigOptiProtTest.getAstexLigPath( "1a07", false ));

        assertTrue( instance.getMolecule()!=null );
    }

    @Test
    public void testCalcObProperties() {
        System.out.println("calcObProperties");

        DockingLigand instance=ConfigOptiProtTest.readAstexLigand( "1a07", false );

        assertTrue( instance.getNumRotableBonds()>0 );

        assertTrue( instance.getInternalEnergy()!=null );
    }

    
}