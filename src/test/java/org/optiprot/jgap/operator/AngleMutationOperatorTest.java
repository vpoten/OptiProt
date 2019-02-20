/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.operator;

import java.util.ArrayList;
import java.util.List;
import org.biojava.bio.structure.StructureException;
import org.jgap.IChromosome;
import org.jgap.InvalidConfigurationException;
import org.jgap.Population;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.jgap.OptiProtConfiguration;
import org.optiprot.jgap.chromosome.ProteinChromFactory;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class AngleMutationOperatorTest {

    public AngleMutationOperatorTest() {
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
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of operate method, of class AngleMutationOperator.
     */
    @Test
    public void testOperate() throws InvalidConfigurationException, StructureException {
        System.out.println("AngleMutationOperator.operate");


        List a_candidateChromosomes = new ArrayList<IChromosome>();

        OptiProtConfiguration conf=new OptiProtConfiguration( ConfigOptiProtTest.getParameters() );

        //new operator with mutation rate 1/n (n=1)
        AngleMutationOperator instance = new AngleMutationOperator(conf, 1);

        //list of chromosome for population
        List<IChromosome> listChrom=new ArrayList<IChromosome>();
        listChrom.add( ProteinChromFactory.create(conf, ConfigOptiProtTest.getChain1()) );

        Population a_population = new Population(conf, 
                (IChromosome [])listChrom.toArray(new IChromosome [1]) );
        
        instance.operate(a_population, a_candidateChromosomes);

        assertTrue( a_candidateChromosomes.size()==1 );
    }

}