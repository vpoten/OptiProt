/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io;

import java.io.File;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.io.mol2.Mol2Atom;
import static org.junit.Assert.*;
import org.optiprot.io.mol2.Mol2Structure;
import org.optiprot.potential.docking.element.DockingLigand;

/**
 *
 * @author victor
 */
public class Mol2ReaderTest {

    public Mol2ReaderTest() {
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
     * Test of parseMol2File method, of class Mol2Reader.
     */
    @Test
    public void testParseMol2File() throws Exception {
        System.out.println("parseMol2File");
        String path = 
                ConfigOptiProtTest.getParameters().getWorkDir()+File.separator+"ligand_1a07.mol2";

        Mol2Structure expResult = null;
        Mol2Structure result = Mol2Reader.parseMol2File(path);
        
        assertTrue( result!=null );
        assertTrue( result.getAtoms().size()==74 );
        assertTrue( result.getChain(0).getAtomLength()==4 );
        assertTrue( result.getAtoms().get(4).getAtomicNum()==6 );
        assertTrue( result.getAtoms().get(68).isHydrogen() );
        assertTrue( result.getBonds().get(5).isDouble() );
        
    }

    @Test
    public void testGenerObConformer() throws Exception {
        System.out.println("generObConformer");
        String path =
                ConfigOptiProtTest.getParameters().getWorkDir()+File.separator+"ligand_1a07.mol2";

        DockingLigand lig=Mol2Reader.generObConformer(path, 5, 15);

        assertTrue( lig!=null );
        assertTrue( lig.getNumRotableBonds()==20 );
        List<Atom> list=lig.getAtoms();

        assertTrue( list.size()==74 );
        assertTrue( ((Mol2Atom)list.get(4)).getAtomicNum()==6 );
        assertTrue( ((Mol2Atom)list.get(68)).isHydrogen() );

    }

}