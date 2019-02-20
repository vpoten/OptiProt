/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.aacids;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.optiprot.OptiProtParameters;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.CalcTransform;

/**
 *
 * @author victor
 */
public class AABasicFactoryTest {

    public AABasicFactoryTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.setName("AABasicFactoryTest");
        ConfigOptiProtTest.loadRotamerLibrary();
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
     * Test of create method, of class AABasicFactory.
     */
    @Test
    public void testCreate_String_OptiProtParameters() throws Exception {
        System.out.println("create");
        
        String sequence =
                "MPVAVGPYGQSQPSCFDRVKMGFVMGCAVGMAAGALFGTFSCLRIGMRGRELMGGIGKTMMQSGGTFGTF" +
                          "MAIGMGIRC";

//        String sequence = "MPVAVGPYGQSQPSCFDRVK";

        OptiProtParameters parameters = ConfigOptiProtTest.getParameters();
        IAABasic[] result = AABasicFactory.create(sequence, parameters);

        assertTrue( CalcTransform.checksBackbone(result, 0, result.length) );

//        ConfigOptiProtTest.writeStructure(
//                AABasicFactory.toChain(result, false), "AABasicFactoryTest");

    }

}