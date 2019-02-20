/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.*;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.potential.IForceField;
import org.optiprot.potential.TestForceField;

/**
 *
 * @author victor
 */
public class CalcSurfVolTest {

    public CalcSurfVolTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {

        ConfigOptiProtTest.setName("CalcSurfVolTest");
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of SAsurface method, of class CalcSurfVol.
     */
    @Test
    public void testSAsurface() throws StructureException {
        System.out.println("SAsurface");

        Chain chain1=ConfigOptiProtTest.getChain1();
        Atom [] molecule=CalcTransform.toAtomArray(chain1);
        IForceField ff = new TestForceField();

        double result=0.0;
        double expResult=0.0;
        double error=0.0;

        expResult = 9561.13346695531;//for eif6_swiss
        result = CalcSurfVol.SAsurface(molecule, ff, 4);
        error=expResult*0.15;
        assertTrue( Math.abs(result-expResult)<error);

        expResult = 21594.5643783;//for eif6_swiss
        result = CalcSurfVol.VDWsurface(molecule, ff, 4);
        //error is a percentage
        error=expResult*0.1;
        assertTrue( Math.abs(result-expResult)<error);

        int [] array_index=new int [] {0};

        chain1=ConfigOptiProtTest.readChain("1ECA.pdb", array_index);
        molecule=CalcTransform.toAtomArray(chain1);
        expResult = 7100;//1ECA
        result = CalcSurfVol.SAsurface(molecule, ff, 4);
        error=expResult*0.15;
        assertTrue( Math.abs(result-expResult)<error);

        expResult = 13852;//1ECA
        result = CalcSurfVol.VDWsurface(molecule, ff, 4);
        error=expResult*0.1;
        assertTrue( Math.abs(result-expResult)<error);

        chain1=ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1});
        molecule=CalcTransform.toAtomArray(chain1);
        expResult = 37453.0;//1CAU
        result = CalcSurfVol.VDWsurface(molecule, ff, 4);
        error=expResult*0.1;
        assertTrue( Math.abs(result-expResult)<error);

        chain1=ConfigOptiProtTest.readChain("1NXB.pdb", array_index);
        molecule=CalcTransform.toAtomArray(chain1);
        expResult = 5673.7;//1NXB
        result = CalcSurfVol.VDWsurface(molecule, ff, 4);
        error=expResult*0.1;
        assertTrue( Math.abs(result-expResult)<error);
    }

}