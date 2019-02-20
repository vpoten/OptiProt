/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import org.biojava.bio.structure.Atom;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcRmsd;
import org.optiprot.metaheuristic.de.impl.BinCrossover;
import org.optiprot.metaheuristic.de.impl.DockScoreEval;
import org.optiprot.metaheuristic.de.impl.RandMutation;
import org.optiprot.metaheuristic.de.impl.TPEvaluator;
import org.optiprot.potential.docking.SimpleDockScore;
import static org.junit.Assert.*;
import org.optiprot.potential.docking.Docked;
import org.optiprot.potential.docking.element.DockingLigand;
import org.optiprot.potential.docking.element.DockingProtein;

/**
 *
 * @author victor
 */
public class DifferentialEvolutionTest {

    public DifferentialEvolutionTest() {
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
    public void testDoSearchTP() {

        System.out.println("DE doSearch Chebychev polynomial");
        DifferentialEvolution instance = new DifferentialEvolution();

        double CR=0.9;
        double F=0.8;
        int NP=50;
        int dim=9;

        double [] lowers=null;
        double [] uppers=null;

        if( dim==5){
            lowers=new double [] { -100, -100, -100, -100, -100 };
            uppers=new double [] { 100, 100, 100, 100, 100 };
        }
        else{
            lowers=new double [] { -1000, -1000, -1000, -1000, -1000,
                    -1000, -1000, -1000, -1000};
            uppers=new double [] { 1000, 1000, 1000, 1000, 1000,
                    1000, 1000, 1000, 1000};
        }


        instance.setCR( CR );
        instance.setF( F );
        instance.setNP( NP );
        instance.setCrossover( new BinCrossover() );
        instance.setMutation( new RandMutation() );
        instance.setBounds(lowers, uppers);

        TPEvaluator evaluator=new TPEvaluator(dim);

        instance.createRandomPopulation( dim, evaluator);

        //instance.setTimeLimit(5);
        instance.setIterLimit(1500);
        instance.doSearch();

        int neval=evaluator.getNevaluations();
        double fitness=instance.getBestSolution().getFitness();

        assertTrue( fitness<0.1 );
    }


    @Test
    public void testToSaveString() throws IOException {

        System.out.println("toSaveString");
        DifferentialEvolution instance = new DifferentialEvolution();

        double CR=0.9;
        double F=0.8;
        int NP=50;
        int dim=9;

        instance.setCR( CR );
        instance.setF( F );
        instance.setNP( NP );

        instance.createRandomPopulation( dim, null, -10, 10);

        String cad=instance.toSaveString();

        int lines=0;

        //count result lines
        for(int i=0;i<cad.length();i++){
            if( cad.charAt(i)=='\n')
                lines++;
        }

        assertTrue( lines==(4+NP*dim) );

        BufferedReader r = new BufferedReader( new StringReader(cad) );
        instance = new DifferentialEvolution( r, null);

        String cad2=instance.toSaveString();

        assertTrue( cad2.equals(cad) );
        assertTrue( instance.getCR()==CR );
        assertTrue( instance.getNP()==NP );

    }

    
    /**
     * Test of doSearch method, of class DifferentialEvolution.
     */
    @Test
    public void testDoSearch() {
        System.out.println("DE doSearch");
        DifferentialEvolution instance = new DifferentialEvolution();

        String pdbcode="1lyb";
        //String pdbcode="1lmo";

        DockingProtein protein=ConfigOptiProtTest.readAstexProtein( pdbcode, null, false);
        DockingLigand ligand=ConfigOptiProtTest.readAstexLigand( pdbcode, true );

        Atom actSiteCoM=ligand.getCoM();

        //read reference ligand to compare
        DockingLigand ligand_ref=ConfigOptiProtTest.readAstexLigand( pdbcode, false );

        //translates ligand to origin
        ligand.shift( CalcGeom.product(actSiteCoM, -1));
        

        double CR=0.9;
        double F=0.8;
        int NP=50;

        double boxLength=ligand.calcRadii()*0.5;
        //double boxLength=0.0;

        double [] lowers=new double [] { -1, -1, -1, -1,
            actSiteCoM.getX()-boxLength, actSiteCoM.getY()-boxLength, actSiteCoM.getZ()-boxLength };
        
        double [] uppers=new double [] { 1, 1, 1, 1,
            actSiteCoM.getX()+boxLength, actSiteCoM.getY()+boxLength, actSiteCoM.getZ()+boxLength};

        instance.setCR( CR );
        instance.setF( F );
        instance.setNP( NP );
        instance.setCrossover( new BinCrossover() );
        instance.setMutation( new RandMutation() );
        instance.setBounds(lowers, uppers);

        Docked docked=new Docked( protein, ligand, actSiteCoM,
                ConfigOptiProtTest.getParameters().getGrid() );

        DockScoreEval evaluator=new DockScoreEval();
        evaluator.setScoreFunc( new SimpleDockScore() );
        evaluator.setDocked(docked);


        instance.createRandomPopulation( DockScoreEval.NUM_PARAMETERS, evaluator);

        ///instance.setHybridFCM(true);
        ///instance.setModifyF(true);

        instance.setTimeLimit(20);
        instance.doSearch();

        //instance.setIterLimit( instance.getIterLimit()+10 );
        //instance.doSearch();

        int neval=evaluator.getScoreFunc().getNevaluations();
        DESolution bestsol=(DESolution) instance.getBestSolution();

        docked.setTransform( evaluator.getMTrans(bestsol) );

        double rmsd=CalcRmsd.bruteRmsd(ligand_ref.getAtoms(), docked.getLigandAtoms());
        double overlap=docked.calcOverlapFactor();
        double fitness=bestsol.recalcFitness();

        assertTrue( rmsd<2.0 );
    }

    
    
}