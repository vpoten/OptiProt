/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.metaheuristic.de.DESolution;
import static org.junit.Assert.*;
import org.optiprot.metaheuristic.de.DifferentialEvolution;
import org.optiprot.metaheuristic.de.IDESolutionEval;
import org.optiprot.metaheuristic.de.impl.BinCrossover;
import org.optiprot.metaheuristic.de.impl.RandMutation;
import org.optiprot.metaheuristic.fans.FANSDouble;
import org.optiprot.metaheuristic.fans.impl.FANSNeighborBase;
import org.optiprot.neural.ElmanNetwork;
import org.optiprot.neural.JordanNetwork;
import org.optiprot.neural.NNPattern;
import org.optiprot.neural.NeuralNetwork;
import org.optiprot.neural.ProdUnitNetwork;

/**
 *
 * @author victor
 */
public class SecStructDeEvaluatorTest {

    public SecStructDeEvaluatorTest() {
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
     * Test of isBetter method, of class SecStructDeEvaluator.
     */
    @Test
    public void testIsBetter() {
        System.out.println("isBetter");
        double val1 = 0.0;
        double val2 = 1.0;
        SecStructDeEvaluator instance = new SecStructDeEvaluator();
        
        boolean result = instance.isBetter(val1, val2);
        assertEquals(true, result);
    }

    /**
     * Test of getNevaluations method, of class SecStructDeEvaluator.
     */
    @Test
    public void testGetNevaluations() {
        System.out.println("getNevaluations");
        SecStructDeEvaluator instance = new SecStructDeEvaluator();
        int expResult = 0;
        int result = instance.getNevaluations();
        assertEquals(expResult, result);
    }


    /**
     * Test of setNetwork method, of class SecStructDeEvaluator.
     */
    @Test
    public void testSetNetwork() {
        System.out.println("setNetwork");

        NeuralNetwork network = new NeuralNetwork(5, 3, 5);
        SecStructDeEvaluator instance = new SecStructDeEvaluator();
        instance.setNetwork(network);

        assertTrue( instance.getNetwork()==network );
    }

//    @Test
//    public void testDETrainNetwork()
//            throws FileNotFoundException, IOException {
//
//        System.out.println("testDETrainNetwork()");
//
//        SecStructPattManager patman=new SecStructPattManager( "/home/victor/recentpdb_codes25.txt",
//                ConfigOptiProtTest.getParameters().getPDBDir());
//
//        int windowSize=11;
//        int resPredPos=5;
//
//        SecStructPattGen1 patgen=new SecStructPattGen1();
//        patgen.setMode( SecStructPattGen1.MODE_COD4 );
//
//        patman.genPatterns(windowSize, resPredPos, patgen );
//
//        NNPattern pattern=patman.getPatGen().createPattern(windowSize);
//
//        int nin=pattern.getNInputs();
//        int nout=pattern.getNOutputs();
//        int nhid=50;
//
//        ///NeuralNetwork network = new JordanNetwork(nin, nout, nhid);
//        ////NeuralNetwork network = new NeuralNetwork(nin, nout, nhid);
//        NeuralNetwork network = new ProdUnitNetwork(nin, nout, nhid);
//
//
//        //prepare DE
//        SecStructDeEvaluator evaluator= new SecStructDeEvaluator(network, patman);
//
//        double CR=0.9;
//        double F=0.8;
//        int dim = network.getNInputs()*network.getNHidden() +
//                network.getNOutputs()*network.getNHidden();
//        int NP = 100;
//
//        //BufferedReader r = new BufferedReader(
//         //       new FileReader("/home/victor/work_bio/optiprot/depopNN_100.txt") );
//
//        //DifferentialEvolution devol = new DifferentialEvolution( r, evaluator);
//
//        DifferentialEvolution devol = new DifferentialEvolution();
//        devol.createRandomPopulation( dim, evaluator, -1, 1 );
//
//
//        devol.setCR( CR );
//        devol.setF( F );
//        devol.setNP( NP );
//        devol.setCrossover( new BinCrossover() );
//        devol.setMutation( new RandMutation() );
//
//        //devol.setModifyF(true);
//        devol.setHybridFCM(true);
//
//        //add FANS local search to DE
//        devol.setFansFCM(true, null, 0.05);
//        devol.setFcmPeriod(15);
//        devol.setFansIter(25);
//        //
//
//
//        ArrayList<IDESolutionEval> listEvals=new ArrayList<IDESolutionEval>();
//
//        int nworkers=3;
//        listEvals.add( evaluator );
//
//        for(int i=1;i<nworkers;i++){
//            ///listEvals.add( new SecStructDeEvaluator( new JordanNetwork(nin, nout, nhid), patman.clone()) );
//            listEvals.add( new SecStructDeEvaluator( new ProdUnitNetwork(nin, nout, nhid), patman.clone()) );
//        }
//
//        devol.setWorkers( listEvals, 250 );
//
//
//        double fitness=java.lang.Double.MAX_VALUE;
//        double min_mse=0.1;
//
//        while( fitness > min_mse  ){
//            devol.setIterLimit( devol.getIterLimit()+10 );
//            devol.doSearch();
//            fitness=devol.getBestSolution().getFitness();
//            fitness /= (double)patman.getNumResidues();
//        }
//
//        int neval=evaluator.getNevaluations();
//
//        assertTrue( fitness<min_mse );
//
//    }

//    @Test
//    public void testClassify()
//            throws FileNotFoundException, IOException {
//
//        System.out.println("testClassify");
//
//        String file_pdbcodes="/home/victor/recentpdb_codes25.txt";
//        String file_pdbcodes2="/home/victor/recentpdb_codes100.txt";
//
//        SecStructPattManager patman=new SecStructPattManager( file_pdbcodes,
//                ConfigOptiProtTest.getParameters().getPDBDir());
//
//        int windowSize=11;
//        int resPredPos=5;
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
//        int nhid=75;
//
//        ////NeuralNetwork network = new JordanNetwork(nin, nout, nhid);
//        NeuralNetwork network = new NeuralNetwork(nin, nout, nhid);
//        ////NeuralNetwork network = new ProdUnitNetwork(nin, nout, nhid);
//        ////NeuralNetwork network = new ElmanNetwork(nin, nout, nhid);
//
//        //prepare DE
//        SecStructDeEvaluator evaluator= new SecStructDeEvaluator(network, patman);
//
//        double CR=0.9;
//        double F=0.8;
//        int dim = network.getNInputs()*network.getNHidden() +
//                network.getNOutputs()*network.getNHidden();
//        int NP = 100;
//
//        BufferedReader r = new BufferedReader(
//                new FileReader("/home/victor/work_bio/optiprot/depopMLPa_100_75.txt") );
//
//        //DifferentialEvolution devol = new DifferentialEvolution();
//        //devol.setNP(NP);
//        //devol.createRandomPopulation( dim, evaluator, -1, 1 );
//
//        DifferentialEvolution devol = new DifferentialEvolution( r, evaluator);
//
//        devol.setCrossover( new BinCrossover() );
//        devol.setMutation( new RandMutation() );
//
//        ////devol.setIterLimit( 1 );
//        ////devol.doSearch();
//
//        network.setWeights( ((DESolution)devol.getBestSolution()).getParameters() );
//
//        double [][]mat=patman.classify( network);
//        //double [][]mat=patman.classify( file_pdbcodes2,
//         //       ConfigOptiProtTest.getParameters().getPDBDir(), network);
//
//        int numres=patman.getNumResidues();
//        int total=0;
//        int correct=0;
//
//        for(int i=0;i<3;i++){
//            for(int j=0;j<3;j++){
//                total+=mat[i][j];
//
//                if(i==j)
//                    correct+=mat[i][j];
//            }
//        }
//
//        double rate=correct/(double)total;
//
//        assertTrue( mat!=null );
//        //assertTrue( total==numres );
//    }
}