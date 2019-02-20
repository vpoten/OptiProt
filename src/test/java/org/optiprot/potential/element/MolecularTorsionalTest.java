/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.StructureException;
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
public class MolecularTorsionalTest {

    static CharmmForceField ffield =
            (CharmmForceField) ConfigOptiProtTest.getParameters().getForceField();


    public MolecularTorsionalTest() {
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
     * Test of getAngle method, of class MolecularGroup.
     */
    @Test
    public void testGetAngle() throws StructureException {
        System.out.println("getAngle");
        MolecularTorsional instance = new MolecularTorsional();

        Atom at1=new AtomImpl();
        at1.setCoords(new double [] {0,0,0});
        at1.setName("CA1");
        instance.setAtom(0, at1, 0, ffield.getAtomType("CA"));

        Atom at2=new AtomImpl();
        at2.setCoords(new double [] {1,0,1});
        at2.setName("CA2");
        instance.setAtom(1, at2, 0, ffield.getAtomType("CA"));

        Atom at3=new AtomImpl();
        at3.setCoords(new double [] {2,0,0});
        at3.setName("CA3");
        instance.setAtom(2, at3, 0, ffield.getAtomType("CA"));

        Atom at4=new AtomImpl();
        at4.setCoords(new double [] {3,0,1});
        at4.setName("CA4");
        instance.setAtom(3, at4, 0, ffield.getAtomType("CA"));

        instance.calcDihedralAng();

        double expResult = 180.0;
        double result = instance.getAngle();

        assertTrue( Math.abs(result-expResult)<0.01 );
    }


}