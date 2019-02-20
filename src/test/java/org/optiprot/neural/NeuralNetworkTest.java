/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.neural;

import java.io.IOException;
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
public class NeuralNetworkTest {

    static ArrayList<NNPattern> patterns=ConfigOptiProtTest.readSNNSFile( "iris_tr3.pat" );

    public NeuralNetworkTest() {
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
    public void testDETraining() {
        System.out.println("DETraining");

        int nin=patterns.get(0).getNInputs();
        int nout=patterns.get(0).getNOutputs();
        int nhid=4;

        NeuralNetwork instance = new NeuralNetwork(nin, nout, nhid);

        //prepare DE
        double CR=0.9;
        double F=0.8;
        int dim = nhid*instance.getNInputs() + nhid*nout;
        int NP = dim*5;
        NP=100;

        DifferentialEvolution devol = new DifferentialEvolution();

        devol.setCR( CR );
        devol.setF( F );
        devol.setNP( NP );
        devol.setCrossover( new BinCrossover() );
        devol.setMutation( new RandMutation() );

        devol.setModifyF(true);

        NeuralEvaluator evaluator=new  NeuralEvaluator();
        evaluator.setNetwork(instance);
        evaluator.setPatterns(patterns);

        devol.createRandomPopulation( dim, evaluator, -0.05, 0.05);

        double fitness=java.lang.Double.MAX_VALUE;
        double min_mse=0.05;

        while( fitness > min_mse  ){
            devol.setIterLimit( devol.getIterLimit()+15 );
            devol.doSearch();
            fitness=devol.getBestSolution().getFitness();
            fitness /= (double)patterns.size();
        }

        int neval=evaluator.getNevaluations();

        ///assertTrue( fitness<min_mse );

        //prediction of the test set
        double [] vector = new double [dim];
        ((DESolution)devol.getBestSolution()).getParameters(vector);

        instance.setWeights(vector);
        instance.init();

        ArrayList<NNPattern> test=ConfigOptiProtTest.readSNNSFile( "iris_v3.pat" );

        double sse=0;

        for( NNPattern pat : test ){
            instance.evaluate( pat.getInputs() );
            sse += instance.calcSSE( pat.getOutputs() );
        }

        sse /= (double)test.size();

        assertTrue( sse < min_mse*1.45 );
    }

    @Test
    public void testToSaveString() throws IOException {

        System.out.println("ToSaveString");

        int nin=5;
        int nout=3;
        int nhid=6;

        NeuralNetwork instance = new NeuralNetwork(nin, nout, nhid);

        String cad=instance.toSaveString();

        int lines=0;

        //count result lines
        for(int i=0;i<cad.length();i++){
            if( cad.charAt(i)=='\n')
                lines++;
        }

        assertTrue( lines==(3+nin*nhid+nout*nhid) );

        instance = new NeuralNetwork(cad);

        assertTrue( instance.getNInputs()==nin );
        assertTrue( instance.getNOutputs()==nout );

        String cad2=instance.toSaveString();

        assertTrue( cad2.equals(cad) );
    }
}