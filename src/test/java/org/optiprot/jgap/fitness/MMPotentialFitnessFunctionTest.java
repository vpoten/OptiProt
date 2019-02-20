/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.fitness;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.jgap.IChromosome;
import org.jgap.InvalidConfigurationException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.jgap.chromosome.ProteinChromFactory;
import org.optiprot.potential.MMPotentialEnergy;
import org.optiprot.potential.MolecularElementsImpl;

/**
 *
 * @author victor
 */
public class MMPotentialFitnessFunctionTest {

    public MMPotentialFitnessFunctionTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.setName("TinkerFitnessFunctionTest");
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
     * Test of evaluate method, of class MMPotentialFitnessFunction.
     */
    @Test
    public void testEvaluate() throws InvalidConfigurationException {
        System.out.println("evaluate");

        IChromosome chromosome=null;
        Chain chain=ConfigOptiProtTest.readChain("1NXB.pdb", new int [] {0});

        try {
            chromosome = ProteinChromFactory.create(ConfigOptiProtTest.getConfiguration(),
                    chain );
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(MMPotentialFitnessFunctionTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testEvaluate ");
        } catch (StructureException ex) {
            Logger.getLogger(MMPotentialFitnessFunctionTest.class.getName()).log(Level.SEVERE, null, ex);
            fail("testEvaluate ");
        }

        MMPotentialFitnessFunction instance = new MMPotentialFitnessFunction( ConfigOptiProtTest.getParameters() );
        ConfigOptiProtTest.getParameters().setGeneratesH(true);


        Double result = instance.evaluate(chromosome);

        chain=ProteinChromFactory.toChain(chromosome, ConfigOptiProtTest.getParameters().isGeneratesH());

        MolecularElementsImpl molEle=new MolecularElementsImpl(chain, ConfigOptiProtTest.getParameters().getGrid(),
                ConfigOptiProtTest.getParameters() );
        double expResult=MMPotentialEnergy.calcEnergy(molEle, ConfigOptiProtTest.getParameters());

        assertTrue( Math.abs(result-expResult)<1e-9);
    }

   
}