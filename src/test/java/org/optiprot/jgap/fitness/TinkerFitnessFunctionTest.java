/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.fitness;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.*;
import org.jgap.IChromosome;
import org.jgap.InvalidConfigurationException;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;
import org.optiprot.jgap.chromosome.ProteinChromFactory;
import org.optiprot.configtests.ConfigOptiProtTest;

/**
 *
 * @author victor
 */
public class TinkerFitnessFunctionTest {


    public TinkerFitnessFunctionTest() {

    }

    @Before
    public void setUp() {
        ConfigOptiProtTest.setName("TinkerFitnessFunctionTest");
        ConfigOptiProtTest.loadRotamerLibrary();
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of evaluate method, of class TinkerFitnessFunction.
     */
    @Test
    public void testEvaluate() {
        
        System.out.println("Test TinkerFitnessFunction.evaluate");
        
        IChromosome chromosome=null;

        try {
            chromosome = ProteinChromFactory.create(ConfigOptiProtTest.getConfiguration(),
                    ConfigOptiProtTest.getChain1());
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(TinkerFitnessFunctionTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testEvaluate ");
        } catch (StructureException ex) {
            Logger.getLogger(TinkerFitnessFunctionTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testEvaluate ");
        }

        TinkerFitnessFunction instance = new TinkerFitnessFunction( ConfigOptiProtTest.getParameters() );


        Double result = instance.evaluate(chromosome);

        assertNotNull(result);
        
    }


}