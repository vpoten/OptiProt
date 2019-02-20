/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.operator;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.*;
import org.jgap.*;
import org.jgap.impl.*;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.jgap.chromosome.ProteinChromFactory;
import org.optiprot.maths.CalcRmsd;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class OneCutCrossoverOperatorTest {


    public OneCutCrossoverOperatorTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.loadRotamerLibrary();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        ConfigOptiProtTest.setName("OneCutCrossoverOperatorTest");
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of doCrossover method, of class OneCutCrossoverOperator.
     */
    @Test
    public void testDoCrossover() {
        System.out.println("Test OneCutCrossoverOperator.doCrossover");
        
        IChromosome firstMate=null;
        IChromosome secondMate=null;

        try {
            firstMate = ProteinChromFactory.create(ConfigOptiProtTest.getConfiguration(),
                    ConfigOptiProtTest.getChain1());
            secondMate = ProteinChromFactory.create(ConfigOptiProtTest.getConfiguration(),
                    ConfigOptiProtTest.getChain1());
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(OneCutCrossoverOperatorTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testDoCrossover -  ");
        } catch (StructureException ex) {
            Logger.getLogger(OneCutCrossoverOperatorTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testDoCrossover -  ");
        }
        
        
        List a_candidateChromosomes = new ArrayList<IChromosome>();
        RandomGenerator generator = new StockRandomGenerator();
        OneCutCrossoverOperator instance=null;
        
        try {
            instance = new OneCutCrossoverOperator(ConfigOptiProtTest.getConfiguration(),
                    ConfigOptiProtTest.getParameters().getCrossOneCutRate());
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(OneCutCrossoverOperatorTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testDoCrossover - OneCutCrossoverOperator() ");
        }

        instance.doCrossover(firstMate, secondMate, a_candidateChromosomes, generator);

        IChromosome offspr1=(IChromosome) a_candidateChromosomes.get(0);
        IChromosome offspr2=(IChromosome) a_candidateChromosomes.get(1);
        double result=0;

        try {
            result = CalcRmsd.rmsd(ProteinChromFactory.toChain(offspr1, false),
                    ProteinChromFactory.toChain(offspr2, false));
        } catch (StructureException ex) {
            Logger.getLogger(OneCutCrossoverOperatorTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testDoCrossover -  CalcRmsd");
        }

        assertTrue(result<0.01);
    }

}