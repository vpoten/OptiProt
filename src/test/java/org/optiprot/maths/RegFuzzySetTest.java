/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class RegFuzzySetTest {

    public RegFuzzySetTest() {
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
     * Test of menbership method, of class RegFuzzySet.
     */
    @Test
    public void testMenbership() {
        System.out.println("menbership");
        
        double[] outputs = new double [5];

        double in = 0.0;
        RegFuzzySet.menbership(in, outputs);
        assertTrue( outputs[0]==1.0 );
        assertTrue( outputs[1]==0.0 );

        in = 0.25;
        RegFuzzySet.menbership(in, outputs);
        assertTrue( outputs[0]==0.0 );
        assertTrue( outputs[1]==1.0 );
        assertTrue( outputs[2]==0.0 );

        in = 0.5;
        RegFuzzySet.menbership(in, outputs);
        assertTrue( outputs[0]==0.0 );
        assertTrue( outputs[1]==0.0 );
        assertTrue( outputs[2]==1.0 );

        in = 0.75;
        RegFuzzySet.menbership(in, outputs);
        assertTrue( outputs[2]==0.0 );
        assertTrue( outputs[3]==1.0 );
        assertTrue( outputs[4]==0.0 );

        in = 1;
        RegFuzzySet.menbership(in, outputs);
        assertTrue( outputs[3]==0.0 );
        assertTrue( outputs[4]==1.0 );

        in = 0.125;
        RegFuzzySet.menbership(in, outputs);
        assertTrue( outputs[0]==0.5 );
        assertTrue( outputs[1]==0.5 );

        in = 0.25+0.0625;
        RegFuzzySet.menbership(in, outputs);
        assertTrue( outputs[1]==0.75 );
        assertTrue( outputs[2]==0.25 );
    }

}