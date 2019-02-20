/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.neural;

import java.util.ArrayList;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.de.DifferentialEvolution;
import org.optiprot.metaheuristic.de.impl.BinCrossover;
import org.optiprot.metaheuristic.de.impl.NeuralEvaluator;
import org.optiprot.metaheuristic.de.impl.RandMutation;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class JordanNetworkTest {

    static ArrayList<NNPattern> patterns=ConfigOptiProtTest.readSNNSFile( "eight_160.pat" );

    public JordanNetworkTest() {
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
     * Test of getNInputs method, of class JordanNetwork.
     */
    @Test
    public void testGetNInputs() {
        System.out.println("getNInputs");

        int nin=patterns.get(0).getNInputs();
        int nout=patterns.get(0).getNOutputs();
        int nhid=14;

        JordanNetwork instance = new JordanNetwork(nin, nout, nhid);

        int expResult = nin+nout;
        int result = instance.getNInputs();
        assertEquals(expResult, result);

        expResult = nin;
        result = instance.getRealNInputs();
        assertEquals(expResult, result);
    }

    /**
     * Test of evaluate method, of class JordanNetwork.
     */
    @Test
    public void testEvaluate() {
        System.out.println("evaluate");

        int nin=patterns.get(0).getNInputs();
        int nout=patterns.get(0).getNOutputs();
        int nhid=14;

        NNPattern pat=patterns.get( (int)(Math.random()*patterns.size()) );

        
        JordanNetwork instance = new JordanNetwork(nin, nout, nhid);
        instance.evaluate( pat.getInputs() );
        
        assertTrue( instance.calcSSE(pat.getOutputs())>=0.0 );

    }

    @Test
    public void testDETraining() {
        System.out.println("DETraining");

        int nin=patterns.get(0).getNInputs();
        int nout=patterns.get(0).getNOutputs();
        int nhid=8;

        JordanNetwork instance = new JordanNetwork(nin, nout, nhid);

        //prepare DE
        double CR=0.9;
        double F=0.8;
        int dim = nhid*instance.getRealNInputs() + 2*nhid*nout;
        int NP = dim*5;
        NP=100;

        DifferentialEvolution devol = new DifferentialEvolution();

        devol.setCR( CR );
        devol.setF( F );
        devol.setNP( NP );
        devol.setCrossover( new BinCrossover() );
        devol.setMutation( new RandMutation() );

        devol.setModifyF(true);
        ///devol.setHybridFCM(true);

        NeuralEvaluator evaluator=new  NeuralEvaluator();
        evaluator.setNetwork(instance);
        evaluator.setPatterns(patterns);

        devol.createRandomPopulation( dim, evaluator, -0.05, 0.05);

        double fitness=java.lang.Double.MAX_VALUE;
        double min_mse=0.06;

        while( fitness > min_mse  ){
            devol.setIterLimit( devol.getIterLimit()+15 );
            devol.doSearch();
            fitness=devol.getBestSolution().getFitness();
            fitness /= (double)patterns.size();
        }

        int neval=evaluator.getNevaluations();

        assertTrue( fitness<min_mse );

        //prediction of the test set
        double [] vector = new double [dim];
        ((DESolution)devol.getBestSolution()).getParameters(vector);

        instance.setWeights(vector);
        instance.init();

        ArrayList<NNPattern> test=ConfigOptiProtTest.readSNNSFile( "eight_016.pat" );

        double sse=0;

        for( NNPattern pat : test ){
            instance.evaluate( pat.getInputs() );
            sse += instance.calcSSE( pat.getOutputs() );
        }

        sse /= (double)test.size();

        assertTrue( sse < min_mse*1.25 );
    }

}