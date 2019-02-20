/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.optiprot.neural.NNPattern;

/**
 *
 * @author victor
 */
public class SecStructPattGen1Test {

    public SecStructPattGen1Test() {
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
     * Test of createPattern method, of class SecStructPattGen1.
     */
    @Test
    public void testCreatePattern() {
        System.out.println("createPattern");
        int windowSize = 5;
        SecStructPattGen1 instance = new SecStructPattGen1();

        NNPattern result = instance.createPattern(windowSize);

        assertTrue( result.getNInputs()==5*4 );
        assertTrue( result.getNOutputs()==3 );
    }

    /**
     * Test of genPattern method, of class SecStructPattGen1.
     */
    @Test
    public void testGenPattern() {
        System.out.println("genPattern");
        SecStructSubPattern[] window = new SecStructSubPattern[3];

        window[0]=new SecStructSubPattern( 0, SecStructPattManager.TYPE_COIL, -1, 0, 0 );
        window[1]=new SecStructSubPattern( 5, SecStructPattManager.TYPE_HELIX, -1, 0, 0 );
        window[2]=new SecStructSubPattern( 10, SecStructPattManager.TYPE_STRAND, -1, 0, 0 );
        
        int resPredPos = 1;

        SecStructPattGen1 instance = new SecStructPattGen1();
        NNPattern pat = instance.createPattern(window.length);

        
        instance.genPattern(window, resPredPos, pat);

        assertTrue( pat.getOutputs()[SecStructPattManager.TYPE_HELIX]==0.9 );
        assertTrue( pat.getOutputs()[SecStructPattManager.TYPE_STRAND]==0.1 );
        assertTrue( pat.getOutputs()[SecStructPattManager.TYPE_COIL]==0.1 );

        assertTrue( pat.getInputs()[0]==1 );
        assertTrue( pat.getInputs()[1]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( pat.getInputs()[2]==0.5 );
        assertTrue( pat.getInputs()[4]==1 );
        assertTrue( pat.getInputs()[5]==1 );
        assertTrue( pat.getInputs()[6]==0.5 );
        assertTrue( pat.getInputs()[9]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( pat.getInputs()[10]==0.5 );

    }

   

}