/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.potential.CharmmForceField;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class MolecularGroupTest {

    static CharmmForceField ffield =
            (CharmmForceField) ConfigOptiProtTest.getParameters().getForceField();

    public MolecularGroupTest() {
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
     * Test of size method, of class MolecularGroup.
     */
    @Test
    public void testSize() {
        System.out.println("size");
        MolecularGroup instance = new MolecularGroup();
        int expResult = 0;
        int result = instance.size();
        assertEquals(expResult, result);
    }

    

    /**
     * Test of getCharge method, of class MolecularGroup.
     */
    @Test
    public void testGetCharge() {
        System.out.println("getCharge");
        Atom atom = new AtomImpl();
        MolecularGroup instance = new MolecularGroup();

        Atom at1=new AtomImpl();
        at1.setCoords(new double [] {0,0,0});
        at1.setName("CA1");
        instance.addAtom(at1, 0.0, ffield.getAtomType("CA"));

        Atom at2=new AtomImpl();
        at2.setCoords(new double [] {1,0,1});
        at2.setName("CA2");
        instance.addAtom(at2, 0.0, ffield.getAtomType("CA"));

        Atom at3=new AtomImpl();
        at3.setCoords(new double [] {2,0,0});
        at3.setName("CA3");
        instance.addAtom(at3, 0.0, ffield.getAtomType("CA"));

        Atom at4=new AtomImpl();
        at4.setCoords(new double [] {3,0,1});
        at4.setName("CA4");
        instance.addAtom(at4, -0.15, ffield.getAtomType("CA"));

        double expResult = 0.0;
        atom.setName("CA1");
        double result = instance.getCharge(atom);
        assertTrue( Math.abs(result-expResult)<1e-15 );

        expResult = -0.15;
        atom.setName("CA4");
        result = instance.getCharge(atom);
        assertTrue( Math.abs(result-expResult)<1e-15 );
        
    }

    

    /**
     * Test of getAtomType method, of class MolecularGroup.
     */
    @Test
    public void testGetAtomType() {
        System.out.println("getAtomType");
        Atom atom = new AtomImpl();
        MolecularGroup instance = new MolecularGroup();

        Atom at1=new AtomImpl();
        at1.setCoords(new double [] {0,0,0});
        at1.setName("CA1");
        instance.addAtom(at1, 0.0, ffield.getAtomType("CA"));

        Atom at2=new AtomImpl();
        at2.setCoords(new double [] {1,0,1});
        at2.setName("N2");
        instance.addAtom(at2, 0.0, ffield.getAtomType("N"));

        Atom at3=new AtomImpl();
        at3.setCoords(new double [] {2,0,0});
        at3.setName("CA3");
        instance.addAtom(at3, 0.0, ffield.getAtomType("CA"));

        Atom at4=new AtomImpl();
        at4.setCoords(new double [] {3,0,1});
        at4.setName("CA4");
        instance.addAtom(at4, -0.15, ffield.getAtomType("CA"));

        Integer expResult = ffield.getAtomType("CA");
        atom.setName("CA1");
        Integer result = instance.getAtomType(atom);
        assertEquals(expResult, result);

        expResult = ffield.getAtomType("N");
        atom.setName("N2");
        result = instance.getAtomType(atom);
        assertEquals(expResult, result);
        
    }

}