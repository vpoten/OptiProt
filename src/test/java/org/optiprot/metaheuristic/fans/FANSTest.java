/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans;

import java.util.Date;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.aacids.AABasicFactory;
import org.optiprot.metaheuristic.fans.impl.ProtFANSNeighbor;
import org.optiprot.metaheuristic.fans.impl.ProtFANSOperator;
import org.optiprot.metaheuristic.fans.impl.ProtFANSSolution;

/**
 *
 * @author victor
 */
public class FANSTest {

    public FANSTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.setName("FANSTest");
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
     * Test of doSearch method, of class FANS.
     */
    @Test
    public void testDoSearch() throws StructureException {
        System.out.println("doSearch");

        Chain chain = ConfigOptiProtTest.readChain("romo1_itasser.pdb", new int [] {0});
        
        ProtFANSSolution initSol = new ProtFANSSolution(
                AABasicFactory.create(chain, ConfigOptiProtTest.getParameters()),
                ConfigOptiProtTest.getParameters() );

        FANS instance = new FANS();

        instance.setSolution(initSol);

        int maxTrials=5;
        ProtFANSNeighbor neigFuzzMang = 
                new ProtFANSNeighbor( 0.2, maxTrials, initSol.getFitness() );
        neigFuzzMang.setLambda( 0.9 );

        instance.setEvaluator(neigFuzzMang);
        instance.setNeighManager(neigFuzzMang);

        int length=chain.getAtomGroups().size();
        int iterLimit=10;
        int secondsLimit=0;

        ProtFANSOperator operMang =
                new ProtFANSOperator(instance, length , iterLimit, secondsLimit);

        instance.setConditions(operMang);
        instance.setOperator(operMang);
        instance.setOperManager(operMang);

        instance.doSearch();

        assertTrue( instance.getBestSolution().isBetter(initSol) );
    }

    @Test
    public void testDoSearchTime() throws StructureException {
        System.out.println("doSearchTime");

        Chain chain = ConfigOptiProtTest.readChain("romo1_itasser.pdb", new int [] {0});

        ProtFANSSolution initSol = new ProtFANSSolution(
                AABasicFactory.create(chain, ConfigOptiProtTest.getParameters()),
                ConfigOptiProtTest.getParameters() );

        FANS instance = new FANS();

        instance.setSolution(initSol);

        int maxTrials=5;
        ProtFANSNeighbor neigFuzzMang =
                new ProtFANSNeighbor( 0.2, maxTrials, initSol.getFitness() );
        neigFuzzMang.setLambda( 0.9 );

        instance.setEvaluator(neigFuzzMang);
        instance.setNeighManager(neigFuzzMang);

        int length=chain.getAtomGroups().size();
        int iterLimit=0;
        int secondsLimit=10;

        ProtFANSOperator operMang =
                new ProtFANSOperator(instance, length , iterLimit, secondsLimit);

        instance.setConditions(operMang);
        instance.setOperator(operMang);
        instance.setOperManager(operMang);

        Date start=new Date();
        instance.doSearch();
        Date end=new Date();

        double time=(end.getTime()-start.getTime())/1000.0;

        assertTrue( Math.abs( time-secondsLimit ) < secondsLimit*0.20 );
    }
    

}