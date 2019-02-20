/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.util.ArrayList;
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
public class SecStructPattIteratorTest {

    public SecStructPattIteratorTest() {
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
     * Test of hasNext method, of class SecStructPattIterator.
     */
    @Test
    public void testHasNext() {
        System.out.println("hasNext");

        
        ArrayList<SecStructSubPattern> list=new ArrayList<SecStructSubPattern>();

        SecStructPattIterator instance = new SecStructPattIterator(list, 5, 0, new SecStructPattGen1(), true);

        boolean result = instance.hasNext();
        assertEquals( false, result);

        list.add( new SecStructSubPattern(0, SecStructPattManager.TYPE_COIL,-1,0,0) );
        instance = new SecStructPattIterator(list, 5, 0, new SecStructPattGen1(), false);

        result = instance.hasNext();
        assertEquals( true, result);

        instance.next();
        result = instance.hasNext();
        assertEquals( false, result);

    }

    /**
     * Test of next method, of class SecStructPattIterator.
     */
    @Test
    public void testNext() {
        System.out.println("next");
        
        ArrayList<SecStructSubPattern> list=new ArrayList<SecStructSubPattern>();
        
        list.add( new SecStructSubPattern(0, SecStructPattManager.TYPE_COIL,-1,0,0) );
        list.add( new SecStructSubPattern(0, SecStructPattManager.TYPE_HELIX,-1,0,0) );

        SecStructPattIterator instance = new SecStructPattIterator(list, 5, 2, new SecStructPattGen1(), false);
        
        NNPattern result = instance.next();

        assertTrue( result.getInputs()[4]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( result.getInputs()[16]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( result.getOutputs()[SecStructPattManager.TYPE_COIL]==0.9 );

        result = instance.next();

        assertTrue( result.getInputs()[5]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( result.getInputs()[13]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( result.getInputs()[15]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( result.getOutputs()[SecStructPattManager.TYPE_HELIX]==0.9 );

    }

    
}