/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
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
public class CalcAlignmentTest {

    public CalcAlignmentTest() {
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
     * Test of align method, of class CalcAlignment.
     */
    @Test
    public void testAlign() throws Exception {
        System.out.println("align");
        
        Chain chain1 = ConfigOptiProtTest.readChainPDB("1bid", new int [] {0});
        Chain chain2 = ConfigOptiProtTest.readChainPDB("3tms", new int [] {0});
  

        double rmsd=CalcRmsd.bruteRmsd(chain1, chain2);

        CalcAlignment.align(chain1, chain2);
        
        double rmsd2=CalcRmsd.bruteRmsd(chain1, chain2);

        assertTrue( rmsd2<=rmsd );
    }

}