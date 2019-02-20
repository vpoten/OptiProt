/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.io.FileNotFoundException;
import java.io.IOException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.fans.FANS;
import org.optiprot.metaheuristic.fans.FANSDouble;
import org.optiprot.metaheuristic.fans.impl.FANSNeighborBase;
import org.optiprot.neural.JordanNetwork;
import org.optiprot.neural.NNPattern;
import org.optiprot.neural.NeuralNetwork;
import org.optiprot.neural.ProdUnitNetwork;

/**
 *
 * @author victor
 */
public class SecStructFANSOperTest {

    public SecStructFANSOperTest() {
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

//    @Test
//    public void testFANSTrainNetwork()
//            throws FileNotFoundException, IOException {
//
//        System.out.println("testFANSTrainNetwork");
//        
//
//        SecStructPattManager patman=new SecStructPattManager( "/home/victor/recentpdb_codes25.txt",
//                ConfigOptiProtTest.getParameters().getPDBDir());
//
//        int windowSize=11;
//        int resPredPos=5;
//        int nhid=75;
//
//        SecStructPattGen1 patgen=new SecStructPattGen1();
//        patgen.setMode( SecStructPattGen1.MODE_COD3 );
//
//        patman.genPatterns(windowSize, resPredPos, patgen );
//
//        NNPattern pattern=patman.getPatGen().createPattern(windowSize);
//
//        int nin=pattern.getNInputs();
//        int nout=pattern.getNOutputs();
//
//        //NeuralNetwork network = new JordanNetwork(nin, nout, nhid);
//        ////NeuralNetwork network = new NeuralNetwork(nin, nout, nhid);
//        NeuralNetwork network = new ProdUnitNetwork(nin, nout, nhid);
//
//        SecStructDeEvaluator evaluator= new SecStructDeEvaluator(network, patman);
//
//        int dim = network.getNInputs()*network.getNHidden() +
//                network.getNOutputs()*network.getNHidden();
//
//        DESolution initSol=new DESolution( dim, evaluator, -1, 1);
//
//        //prepare FANS
//
//        FANS fans = new FANSDouble(dim, evaluator);
//
//        fans.setSolution(initSol);
//
//        int maxTrials=5;
//        FANSNeighborBase neigFuzzMang =
//                new FANSNeighborBase( 0.12, maxTrials, initSol.getFitness() );
//        neigFuzzMang.setLambda( 0.9 );
//
//        fans.setEvaluator(neigFuzzMang);
//        fans.setNeighManager(neigFuzzMang);
//
//
//        SecStructFANSOper operMang = new SecStructFANSOper(fans, 0.05);
//
//        fans.setConditions(operMang);
//        fans.setOperator(operMang);
//        fans.setOperManager(operMang);
//
//
//        double fitness=java.lang.Double.MAX_VALUE;
//        double min_mse=0.1;
//
//        while( fitness > min_mse  ){
//            fans.setIterLimit( fans.getIterLimit()+100 );
//            fans.doSearch();
//            fitness=fans.getBestSolution().getFitness();
//            fitness /= (double)patman.getNumResidues();
//        }
//
//        int neval=evaluator.getNevaluations();
//
//        assertTrue( fitness<min_mse );
//
//     }

}