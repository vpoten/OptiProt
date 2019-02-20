/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.List;
import org.biojava.bio.structure.StructureException;
import org.jgap.InvalidConfigurationException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.jgap.OptiProtConfiguration;
import org.optiprot.jgap.chromosome.ProteinChromFactory;
import org.optiprot.jgap.chromosome.ProteinChromosome;
import org.optiprot.metaheuristic.fans.FANS;
import org.optiprot.metaheuristic.fans.impl.ProtFANSOperator;
import org.optiprot.metaheuristic.fans.impl.ProtFANSSolution;
import static org.junit.Assert.*;
import org.optiprot.rotamer.RotamerLibrary;

/**
 *
 * @author victor
 */
public class CalcClashesTest {

    public CalcClashesTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.setName("CalcClashesTest");
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
     * Test of isSelfClash method, of class CalcClashes.
     */
    @Test
    public void testIsSelfClash() throws InvalidConfigurationException, StructureException {
        System.out.println("isSelfClash");

        OptiProtConfiguration conf=new OptiProtConfiguration( ConfigOptiProtTest.getParameters() );
        RotamerLibrary rot = ConfigOptiProtTest.getParameters().getRotLib();

        //eif6
        ProteinChromosome chrom=(ProteinChromosome) ProteinChromFactory.create( conf, ConfigOptiProtTest.getChain1());

        int start = chrom.getGenes().length/2;
        int end = chrom.getGenes().length;
        
        boolean expResult = false;
        boolean result = CalcClashes.isSelfClash(chrom.toAAarray(), start, end, rot);

        assertEquals(expResult, result);

        //1CAU
        chrom=(ProteinChromosome) ProteinChromFactory.create( conf,
                ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1}) );

        start = chrom.getGenes().length/2;
        end = chrom.getGenes().length;

        expResult = false;
        result = CalcClashes.isSelfClash(chrom.toAAarray(), start, end, rot);

        assertEquals(expResult, result);

        //1NXB
        chrom=(ProteinChromosome) ProteinChromFactory.create( conf,
                ConfigOptiProtTest.readChain("1NXB.pdb", new int [] {0}) );

        start = chrom.getGenes().length/2;
        end = chrom.getGenes().length;

        expResult = false;
        result = CalcClashes.isSelfClash(chrom.toAAarray(), start, end, rot);

        assertEquals(expResult, result);

        //3TRX
        chrom=(ProteinChromosome) ProteinChromFactory.create( conf,
                ConfigOptiProtTest.readChain("3TRX.pdb", new int [] {0}) );

        start = chrom.getGenes().length/2;
        end = chrom.getGenes().length;

        expResult = false;
        result = CalcClashes.isSelfClash(chrom.toAAarray(), start, end, rot);

        assertEquals(expResult, result);
        
    }

    @Test
    public void testClashDetection() throws InvalidConfigurationException, StructureException{
        System.out.println("clashDetection");

        double lambda=0.1;

        OptiProtConfiguration conf=new OptiProtConfiguration( ConfigOptiProtTest.getParameters() );
        RotamerLibrary rot = ConfigOptiProtTest.getParameters().getRotLib();

        //eif6
        ProteinChromosome chrom=(ProteinChromosome) ProteinChromFactory.create( conf, ConfigOptiProtTest.getChain1());

        List<AAClash> clashes = CalcClashes.clashDetection( lambda, chrom.toAAarray(), rot);

        assertTrue( clashes.size()==0 );

        //1CAU
        chrom=(ProteinChromosome) ProteinChromFactory.create( conf,
                ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1}) );

        clashes = CalcClashes.clashDetection( lambda, chrom.toAAarray(), rot);

        assertTrue( clashes.size()==0 );
    }

    @Test
    public void testFixClashes() throws InvalidConfigurationException, StructureException{
        System.out.println("fixClashes");

        double lambda=0.1;

        OptiProtConfiguration conf=new OptiProtConfiguration( ConfigOptiProtTest.getParameters() );
        RotamerLibrary rot = ConfigOptiProtTest.getParameters().getRotLib();
        int maxTrials=1;

        //eif6
        ProteinChromosome chrom=(ProteinChromosome) ProteinChromFactory.create( conf, ConfigOptiProtTest.getChain1());

        CalcClashes.fixClashes( lambda, chrom.toAAarray(), rot, maxTrials );

        List<AAClash> clashes = CalcClashes.clashDetection( lambda, chrom.toAAarray(), rot);

        assertTrue( clashes.size()==0 );

        //1CAU
        chrom=(ProteinChromosome) ProteinChromFactory.create( conf,
                ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1}) );

        CalcClashes.fixClashes( lambda, chrom.toAAarray(), rot, maxTrials );

        clashes = CalcClashes.clashDetection( lambda, chrom.toAAarray(), rot);

        assertTrue( clashes.size()==0 );

        //modify 1CAU
        lambda=0.6;
        
        ProtFANSOperator op=new ProtFANSOperator( new FANS(), chrom.toAAarray().length, 0, 0);

        ProtFANSSolution initSol = new ProtFANSSolution( chrom.toAAarray(),
                ConfigOptiProtTest.getParameters() );

        ProtFANSSolution newSol=(ProtFANSSolution) op.modify(initSol);

        maxTrials=100;
        CalcClashes.fixClashes( lambda, newSol.getChain(), rot, maxTrials );
        clashes = CalcClashes.clashDetection( lambda, newSol.getChain(), rot);

        assertTrue( clashes.size()<15 );

    }

}