/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de.impl;

import org.biojava.bio.structure.Atom;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcRmsd;
import static org.junit.Assert.*;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.potential.docking.SimpleDockScore;
import org.optiprot.potential.docking.Docked;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;

/**
 *
 * @author victor
 */
public class DockScoreEvalTest {

    public DockScoreEvalTest() {
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
     * Test of isBetter method, of class DockScoreEval.
     */
    @Test
    public void testIsBetter() {
        System.out.println("isBetter");

        String pdbcode="1a07";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        Atom actSiteCoM=ligand.getCoM();

        //translates ligand to origin
        ligand.shift( CalcGeom.product(actSiteCoM, -1));

        double boxLength=8.0;

        double [] lowers=new double [] { -1, -1, -1, -1,
            actSiteCoM.getX()-boxLength, actSiteCoM.getY()-boxLength, actSiteCoM.getZ()-boxLength };

        double [] uppers=new double [] { 1, 1, 1, 1,
            actSiteCoM.getX()+boxLength, actSiteCoM.getY()+boxLength, actSiteCoM.getZ()+boxLength};

        Docked docked=new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );
        
        DockScoreEval instance = new DockScoreEval();
        instance.setDocked(docked);
        instance.setScoreFunc( new SimpleDockScore() );

        //optimal solution
        DESolution solOpt=new DESolution( DockScoreEval.NUM_PARAMETERS, instance);
        solOpt.setParameter(0, 1);//cos(0)
        solOpt.setParameter(1, 1);
        solOpt.setParameter(2, 0);
        solOpt.setParameter(3, 0);
        solOpt.setParameter(4, actSiteCoM.getX());
        solOpt.setParameter(5, actSiteCoM.getY());
        solOpt.setParameter(6, actSiteCoM.getZ());
        
        //random solution
        DESolution solRand=new DESolution( DockScoreEval.NUM_PARAMETERS, instance,
                lowers, uppers);

        double val1=instance.getFitness(solOpt);
        double val2=instance.getFitness(solRand);

        //read reference ligand to compare
        DockingLigand ligand_ref=ConfigOptiProtTest.readAstexLigand( pdbcode, false );
        docked.setTransform( instance.getMTrans(solOpt) );
        double rmsd=CalcRmsd.bruteRmsd(ligand_ref.getAtoms(), docked.getLigandAtoms());
        double overlap=docked.calcOverlapFactor();
        
        assertTrue( rmsd<0.001 );

        assertTrue( instance.isBetter(val1, val2) );
        
    }

   
   
}