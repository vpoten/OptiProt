/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcSurfVol;
import org.optiprot.maths.CalcTransform;
import static org.junit.Assert.*;
import org.optiprot.rotamer.RotamerLibrary;

/**
 *
 * @author victor
 */
public class NVSASACalcTest {

    static BSPTree bsptree=null;

    public NVSASACalcTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.setName("NVSASACalcTest");
        ConfigOptiProtTest.loadRotamerLibrary();

        bsptree=new BSPTree(
                ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1}) );
    }

    @Before
    public void setUp() {
        
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of calc method, of class NVSASACalc.
     */
    @Test
    public void testCalc() throws StructureException {
        System.out.println("calc");

        Chain chain=ConfigOptiProtTest.getChain1();
        RotamerLibrary rotLib = ConfigOptiProtTest.getParameters().getRotLib();

        double expResult = 0.0;
        double result = 0;
        double error = 0;
        Atom [] molecule=null;
        IForceField ff = new TestForceField();

//        result = NVSASACalc.calc(chain, rotLib);
//        molecule=CalcTransform.toAtomArray(chain);
//        expResult = CalcSurfVol.SAsurface(molecule, ff, 3);
//        error=expResult*0.25;
//        assertTrue( Math.abs(result-expResult)<error);

        int [] array_index=new int [] {0};

//        chain=ConfigOptiProtTest.readChain("1ECA.pdb", array_index);
//        result = NVSASACalc.calc(chain, rotLib);
//        molecule=CalcTransform.toAtomArray(chain);
//        expResult = CalcSurfVol.SAsurface(molecule, ff, 3);
//        error=expResult*0.25;
//        assertTrue( Math.abs(result-expResult)<error);

        chain=ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1});
        result = NVSASACalc.calc(chain, rotLib);
        molecule=CalcTransform.toAtomArray(chain);
        expResult = CalcSurfVol.SAsurface(molecule, ff, 3);
        error=expResult*0.25;
        assertTrue( Math.abs(result-expResult)<error);
    }

    @Test
    public void testCalc2() throws StructureException {
        System.out.println("calc2");

        Chain chain=ConfigOptiProtTest.getChain1();
        RotamerLibrary rotLib = ConfigOptiProtTest.getParameters().getRotLib();

        double expResult = 0.0;
        double result = 0;
        double error = 0;

        chain=ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1});
        result = NVSASACalc.calc(chain, rotLib);
        expResult = result;
        error=expResult*1e-12;
        assertTrue( Math.abs(result-expResult)<error);
    }

    @Test
    public void testCalcBSP() throws StructureException {
        System.out.println("calcBSP");

        Chain chain=ConfigOptiProtTest.getChain1();
        RotamerLibrary rotLib = ConfigOptiProtTest.getParameters().getRotLib();

        double expResult = 0.0;
        double result = 0;
        double error = 0;

        chain=ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1});
        
        result = NVSASACalc.calc(chain, bsptree, rotLib);
        expResult = result;
        error=expResult*1e-12;
        assertTrue( Math.abs(result-expResult)<error);
    }

}