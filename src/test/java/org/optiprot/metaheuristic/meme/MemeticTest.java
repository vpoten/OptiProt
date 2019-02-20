/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.meme;

import org.biojava.bio.structure.Atom;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcRmsd;
import org.optiprot.metaheuristic.de.impl.DockScoreEval;
import org.optiprot.metaheuristic.meme.impl.MemeOperatorDock;
import org.optiprot.metaheuristic.tabu.impl.DockConfSolution;
import org.optiprot.potential.docking.Docked;
import org.optiprot.potential.docking.SimpleDockScore;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class MemeticTest {

    public MemeticTest() {
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
     * Test of doSearch method, of class Memetic.
     */
    @Test
    public void testDoSearch() {
        System.out.println("doSearch");
        Memetic instance = new Memetic();

        String pdbcode="1acm";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        //read reference ligand to compare
        DockingLigand ligand_ref=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        //translates ligand to origin
        ligand.shift( CalcGeom.product(actSiteCoM, -1));

        Docked docked=new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        //constructs the evaluator and set it to docked
        DockScoreEval evaluator = new DockScoreEval();
        evaluator.setScoreFunc( new SimpleDockScore() );
        evaluator.setDocked( docked );
        docked.setEvaluator( evaluator );


        MemeOperatorDock operator=new MemeOperatorDock();
        operator.setInitIterations(40);
        operator.setNMutants(2);
        operator.setDocked(docked);


        instance.setOperators(operator);

        instance.setNP(5);
        instance.setTimeLimit(50);

        instance.doSearch();
        
        DockConfSolution bestSol=(DockConfSolution) instance.getBestSolution();
        bestSol.optimize(40);

        docked.assignAtomToActivePoint( bestSol.getLigActiveAtom(), bestSol.getActivePoint());
        docked.setTransform( bestSol.getMTrans() );

        double rmsd=CalcRmsd.bruteRmsd(ligand_ref.getAtoms(), docked.getLigandAtoms());
        double overlap=docked.calcOverlapFactor();
        double fitness=bestSol.recalcFitness();

        assertTrue( rmsd<2.0 );
    }

    

}