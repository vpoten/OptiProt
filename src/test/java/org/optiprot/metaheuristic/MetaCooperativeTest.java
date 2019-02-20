/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic;

import java.util.ArrayList;
import java.util.Date;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.OptiProtParameters;
import org.optiprot.aacids.AABasicFactory;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.metaheuristic.fans.FANS;
import org.optiprot.metaheuristic.fans.impl.ProtFANSNeighbor;
import org.optiprot.metaheuristic.fans.impl.ProtFANSOperator;
import org.optiprot.metaheuristic.fans.impl.ProtFANSSolution;

/**
 *
 * @author victor
 */
public class MetaCooperativeTest {

    public MetaCooperativeTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.setName("MetaCooperativeTest");
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
     * Test of getSolMemory method, of class MetaCooperative.
     */
    @Test
    public void testGetSolMemory() {
        System.out.println("getSolMemory");
        MetaCooperative instance = new MetaCooperative();

        ArrayList<IMetaheuristicSolution> result = instance.getSolMemory();

        assertTrue( result.size()==MetaCooperative.MEMORY_LENGTH );
        assertTrue( result.get(MetaCooperative.MEMORY_LENGTH-1) == null );
    }

    /**
     * Test of getImproveMemory method, of class MetaCooperative.
     */
    @Test
    public void testGetImproveMemory() {
        System.out.println("getImproveMemory");
        MetaCooperative instance = new MetaCooperative();
        
        ArrayList<Double> result = instance.getImproveMemory();

        assertTrue( result.size()==MetaCooperative.MEMORY_LENGTH );
        assertTrue( result.get(MetaCooperative.MEMORY_LENGTH-1) ==
                -java.lang.Double.MAX_VALUE );
    }


    /**
     * Test of addSolver method, of class MetaCooperative.
     */
    @Test
    public void testAddSolver() {
        System.out.println("addSolver");

        IMetaheuristic solver = new FANS();
        MetaCooperative instance = new MetaCooperative();

        instance.addSolver(solver);

        assertTrue( instance.getSolvers().size()==1 );
        assertTrue( instance.getSolversInfo().size()==1 );

        solver = new FANS();
        instance.addSolver(solver);

        assertTrue( solver.getId()==1 );
        assertTrue( instance.getSolversInfo().size()==2 );
    }

    /**
     * Test of storeInMemory method, of class MetaCooperative.
     */
    @Test
    public void testStoreInMemory() {
        System.out.println("storeInMemory");
        IMetaheuristicSolution sol = null;

        Double improve = 1000.0;
        MetaCooperative instance = new MetaCooperative();
        instance.storeInMemory(sol, improve);

        assertTrue(
                instance.getImproveMemory().get(MetaCooperative.MEMORY_LENGTH-1) == improve );

        improve = 1300.0;
        instance.storeInMemory(sol, improve);

        assertTrue(
                instance.getImproveMemory().get(MetaCooperative.MEMORY_LENGTH-1) == improve );

    }

    /**
     * Test of getPercentile method, of class MetaCooperative.
     */
    @Test
    public void testGetPercentile_double() {
        System.out.println("getPercentile");
        Double improve = 0.0;
        MetaCooperative instance = new MetaCooperative();
        
        improve=1200.0;
        instance.storeInMemory(null, improve);
        improve=1100.0;
        instance.storeInMemory(null, improve);

        double expResult = 100.0*(MetaCooperative.MEMORY_LENGTH-1) /
                (double)MetaCooperative.MEMORY_LENGTH;
        improve=1150.0;
        double result = instance.getPercentile(improve);

        assertTrue( expResult==result );

        expResult = 100.0*(MetaCooperative.MEMORY_LENGTH-2) /
                (double)MetaCooperative.MEMORY_LENGTH;
        improve=1000.0;
        result = instance.getPercentile(improve);

        assertTrue( expResult==result );
    }

    /**
     * Test of getPercentile method, of class MetaCooperative.
     */
    @Test
    public void testGetPercentile_IMetaheuristicSolution() {
        System.out.println("getPercentile");
        IMetaheuristicSolution sol = null;
        MetaCooperative instance = new MetaCooperative();
        
        double expResult = 100.0;
        double result = instance.getPercentile(sol);

        assertTrue( expResult==result );
    }

    /**
     * Test of isBad method, of class MetaCooperative.
     */
    @Test
    public void testIsBad() {
        System.out.println("isBad");
        double value = MetaCooperative.BETA+0.01;
        MetaCooperative instance = new MetaCooperative();
        double expResult = 0.0;
        double result = instance.isBad(value);

        assertTrue(expResult==result);

        value = MetaCooperative.ALPHA-0.01;
        expResult = 1.0;
        result = instance.isBad(value);

        assertTrue(expResult==result);

        value = (MetaCooperative.ALPHA+MetaCooperative.BETA)*0.5;
        expResult = 0.5;
        result = instance.isBad(value);

        assertTrue(expResult==result);
        
    }

    /**
     * Test of evalFuzzyAntecedent method, of class MetaCooperative.
     */
    @Test
    public void testEvalFuzzyAntecedent() {
        System.out.println("evalFuzzyAntecedent");
        IMetaheuristicSolution sol = null;
        double improve = 10.0;
        double lambda = MetaCooperative.LAMBDA;
        MetaCooperative instance = new MetaCooperative();
        boolean expResult = false;
        boolean result = instance.evalFuzzyAntecedent(sol, improve, lambda);

        assertEquals(expResult, result);
    }

    /**
     * Test of getCurrBestSolver method, of class MetaCooperative.
     */
    @Test
    public void testGetCurrBestSolver() {
        System.out.println("getCurrBestSolver");
        MetaCooperative instance = new MetaCooperative();
        Integer expResult = null;
        Integer result = instance.getCurrBestSolver();

        assertEquals(expResult, result);
    }

    /**
     * Test of sendReport method, of class MetaCooperative.
     */
    @Test
    public void testSendReport() throws StructureException {
        System.out.println("sendReport");

        ProtFANSSolution initSol = new ProtFANSSolution(
                AABasicFactory.create(
                ConfigOptiProtTest.getChain1(),
                ConfigOptiProtTest.getParameters()),
                ConfigOptiProtTest.getParameters() );

        IMetaheuristic solver = new FANS();
        solver.setSolution(initSol);
        MetaCooperative instance = new MetaCooperative();

        instance.addSolver(solver);
        boolean command = 
                instance.sendReport( solver.getId(), solver.getSolution(), new Date());

        assertTrue( !command );

        
    }

    /**
     * Test of runSolvers method, of class MetaCooperative.
     */
    @Test
    public void testRunSolvers() throws StructureException, Exception {
        System.out.println("runSolvers");
        
        Chain chain = ConfigOptiProtTest.readChain("romo1_itasser.pdb", new int [] {0});

        MetaCooperative coordinator =  new MetaCooperative();

        //parameters for solvers
        int maxTrials=4;
        int length=chain.getAtomGroups().size();
        int iterLimit=50;
        int secondsLimit=0;

        //create the solvers

        //solver 1

        //create a object parameter for each solver
        OptiProtParameters par = OptiProtParameters.createParameters(
                ConfigOptiProtTest.getParameters().getWorkDir(),
                "/home/victor/work_bio/rotPDBs" );

        par.setCalcGBorn(false);
        par.setCalcMMechanics(false);

        ProtFANSSolution initSol = new ProtFANSSolution(
                AABasicFactory.create(chain, par), par );

        FANS instance = new FANS();
        ProtFANSNeighbor neigFuzzMang =
                new ProtFANSNeighbor( 0.2, maxTrials, initSol.getFitness() );
        neigFuzzMang.setLambda( 0.9 );

        instance.setEvaluator(neigFuzzMang);
        instance.setNeighManager(neigFuzzMang);
        ProtFANSOperator operMang =
                new ProtFANSOperator(instance, length , iterLimit, secondsLimit);

        instance.setConditions(operMang);
        instance.setOperator(operMang);
        instance.setOperManager(operMang);
        
        instance.setSolution( initSol );

        coordinator.addSolver(instance);

        //solver 2
        instance = new FANS();
        neigFuzzMang =
                new ProtFANSNeighbor( 0.2, maxTrials, initSol.getFitness() );
        neigFuzzMang.setLambda( 0.7 );

        instance.setEvaluator(neigFuzzMang);
        instance.setNeighManager(neigFuzzMang);
        operMang =
                new ProtFANSOperator(instance, length , iterLimit, secondsLimit);

        instance.setConditions(operMang);
        instance.setOperator(operMang);
        instance.setOperManager(operMang);

        //create a object parameter for each solver
        par = OptiProtParameters.createParameters(
                ConfigOptiProtTest.getParameters().getWorkDir(),
                "/home/victor/work_bio/rotPDBs" );

        par.setCalcGBorn(false);
        par.setCalcMMechanics(false);

        instance.setSolution( new ProtFANSSolution(
                AABasicFactory.create(chain, par), par )
                );

        coordinator.addSolver(instance);
        

        //run the solvers in cooperative mode
        try{
            MetaCooperative.runSolvers(coordinator, 6000);
        }
        catch(Exception e){
            fail("Test failed.");
        }

        assertTrue( coordinator.getBestSolution().isBetter(initSol) );
    }

}