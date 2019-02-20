/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.tabu;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.io.PDBParseException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcRmsd;
import org.optiprot.metaheuristic.de.impl.DockScoreEval;
import static org.junit.Assert.*;
import org.optiprot.metaheuristic.tabu.impl.DockConfSolution;
import org.optiprot.metaheuristic.tabu.impl.DockTabuNeighbor;
import org.optiprot.potential.docking.Docked;
import org.optiprot.potential.docking.SimpleDockScore;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;
import org.optiprot.potential.docking.element.DockingProteinTest;

/**
 *
 * @author victor
 */
public class TabuSearchTest {

    public TabuSearchTest() {
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
     * Test of doSearch method, of class TabuSearch.
     */
    @Test
    public void testDoSearch() throws PDBParseException {
        
        System.out.println("doSearch");
        TabuSearch instance = new TabuSearch();

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

        //create initial solution
        DockConfSolution initSol = new DockConfSolution( docked,
                DockConfSolution.MODE_OPT_POSE, 50 );

        instance.setSolution(initSol);
        instance.setNeighManager( new DockTabuNeighbor() );
        instance.setTenure(10);

        instance.setTimeLimit(50);
        //instance.setIterLimit(10);

        instance.doSearch();
        
        
        DockConfSolution bestSol=(DockConfSolution) instance.getBestSolution();
        bestSol.optimize(10);

        docked.assignAtomToActivePoint( bestSol.getLigActiveAtom(), bestSol.getActivePoint());
        docked.setTransform( bestSol.getMTrans() );

        double rmsd=CalcRmsd.bruteRmsd(ligand_ref.getAtoms(), docked.getLigandAtoms());
        double overlap=docked.calcOverlapFactor();
        double fitness=bestSol.recalcFitness();

        DockingProtein protein2=ConfigOptiProtTest.readAstexProtein( pdbcode, null, true);
        //DockingProteinTest.drawProtLigand( protein2.getChain() , docked.getLigandAtoms() );
        ////DockingProteinTest.drawProtLigand( protein2.getChain() , ligand_ref.getAtoms() );

        assertTrue( rmsd<2.0 );

    }

    

}