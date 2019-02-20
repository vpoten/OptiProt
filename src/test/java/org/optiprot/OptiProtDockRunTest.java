/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot;

import java.util.ArrayList;
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
public class OptiProtDockRunTest {

    public OptiProtDockRunTest() {
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
     * Test of run method, of class OptiProtDockRun.
     */
    @Test
    public void testRun() {
        System.out.println("run");
        
        ArrayList<String> pdbcodes=new ArrayList<String>();
         pdbcodes.add("1mup");
        //pdbcodes.add("1lkk");
        pdbcodes.add("1blh");
        pdbcodes.add("1fen");
        pdbcodes.add("1a1b");
        pdbcodes.add("1lpm");
        //pdbcodes.add("1a0q");
        //pdbcodes.add("1rt2");
        //pdbcodes.add("1sln");

        //OptiProtDockRun instance = new OptiProtDockRun( 2, 5, 30 );
        OptiProtDockRun instance = new OptiProtDockRun( 50, 0.9, 0.8 );
        instance.setPdbcodes(pdbcodes);
        instance.setAstexDir( ConfigOptiProtTest.astex_path );
        instance.setTimeLimit(30);

        instance.run();

        assertTrue( instance.getHits()>0 );
        
    }

    

}