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

/**
 *
 * @author victor
 */
public class CalcRmsdTest {

    

    public CalcRmsdTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        
        
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        ConfigOptiProtTest.setName("CalcRmsdTest");

        //rotate chain1 randomly

        double angle=Math.random()*Math.PI;
        Atom vaxis=new AtomImpl();
        vaxis.setX(0.5);
        vaxis.setY(Math.random()-0.5 );
        vaxis.setZ(Math.random()-0.5 );

        ConfigOptiProtTest.setChain2( (Chain) ConfigOptiProtTest.getChain1().clone());

        Atom center= CalcGeom.getCentroid(ConfigOptiProtTest.getChain1());

        CalcTransform.rotateSubChain( center, vaxis, angle,
            ConfigOptiProtTest.getChain1(), 0,
            ConfigOptiProtTest.getChain1().getAtomGroups().size() );

    }

    @After
    public void tearDown() {
    }

    /**
     * Test of rmsd method, of class CalcRmsd.
     */
    @Test
    public void testRmsd() throws Exception {
        System.out.println("Test Calc.rmsd");
        
        
        double result = CalcRmsd.rmsd(ConfigOptiProtTest.getChain1(),
                ConfigOptiProtTest.getChain2());

        assertTrue(result<0.01);
    }

    /**
     * Test of rmsdCA method, of class CalcRmsd.
     */
    @Test
    public void testRmsdCA() throws Exception {
        System.out.println("Test Calc.rmsdCA");
        
        
        double result = CalcRmsd.rmsdCA(ConfigOptiProtTest.getChain1(),
                ConfigOptiProtTest.getChain2());

        assertTrue(result<0.01);
    }

    
   
}