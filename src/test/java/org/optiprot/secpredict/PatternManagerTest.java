/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.neural.NNPattern;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class PatternManagerTest {

    public PatternManagerTest() {
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

    @Test
    public void testConstructor() 
            throws FileNotFoundException, IOException {

        System.out.println("constructor");

        SecStructPattManager instance=new SecStructPattManager( "/home/victor/recentpdb_codes100.txt",
                ConfigOptiProtTest.getParameters().getPDBDir());

        assertTrue( instance.getNumChains()>100*0.9 );

    }

    @Test
    public void testGetIterator4Chain()
            throws FileNotFoundException, IOException {

        System.out.println("getIterator4Chain");

        SecStructPattManager instance=new SecStructPattManager( "/home/victor/recentpdb_codes100.txt",
                ConfigOptiProtTest.getParameters().getPDBDir());

        int windowSize=7;
        int resPredPos=0;

        instance.genPatterns(windowSize, resPredPos, new SecStructPattGen1() );

        Iterator<NNPattern> it=instance.getIterator4Chain(0);

        assertTrue(it!=null);

        NNPattern pat=null;

        while(it.hasNext()){
            pat=it.next();
        }

        assertTrue( pat.getInputs()[4]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( pat.getInputs()[5]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( pat.getInputs()[13]==SecStructPattGen1.ZERO_VALUE );
        assertTrue( pat.getInputs()[14]==SecStructPattGen1.ZERO_VALUE );
    }


}