/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.potential.CharmmForceField;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class CharmmResidueTest {

    static CharmmForceField ffield =
            (CharmmForceField) ConfigOptiProtTest.getParameters().getForceField();

    public CharmmResidueTest() {
    }

    
    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of getAtomType method, of class CharmmResidue.
     */
    @Test
    public void testGetAtomType() {
        System.out.println("getAtomType");

        CharmmResidue instance = new CharmmResidue();
        MolecularGroup group = new MolecularGroup();
        

        Atom atom=new AtomImpl();
        atom.setName( "CA1" );
        group.addAtom(atom, -0.2, ffield.getAtomType("CA"));

        atom=new AtomImpl();
        atom.setName( "CA2" );
        group.addAtom(atom, 0.25, ffield.getAtomType("CB"));

        atom=new AtomImpl();
        atom.setName( "CA3" );
        group.addAtom(atom, 0.5, ffield.getAtomType("CB"));

        atom=new AtomImpl();
        atom.setName( "CA4" );
        group.addAtom(atom, 0.25, ffield.getAtomType("CZ"));
        
        instance.addGroup(group);

        
        Integer expResult = ffield.getAtomType("CA");
        Integer result = instance.getAtomType("CA1");
        assertEquals(expResult, result);

        expResult = ffield.getAtomType("CZ");
        result = instance.getAtomType("CA4");
        assertEquals(expResult, result);

        expResult = ffield.getAtomType("CB");
        result = instance.getAtomType("CA2");
        assertEquals(expResult, result);
        
    }

    /**
     * Test of getAtomCharge method, of class CharmmResidue.
     */
    @Test
    public void testGetAtomCharge() {
        System.out.println("getAtomCharge");
        String name = "";
        
        CharmmResidue instance = new CharmmResidue();
        MolecularGroup group = new MolecularGroup();
        
        
        Atom atom=new AtomImpl();
        atom.setName( "CA1" );
        group.addAtom(atom, -0.2, ffield.getAtomType("CA"));
        
        atom=new AtomImpl();
        atom.setName( "CA2" );
        group.addAtom(atom, 0.25, ffield.getAtomType("CB"));
        
        atom=new AtomImpl();
        atom.setName( "CA3" );
        group.addAtom(atom, 0.5, ffield.getAtomType("CB"));
        
        atom=new AtomImpl();
        atom.setName( "CA4" );
        group.addAtom(atom, 0.25, ffield.getAtomType("CZ"));

        instance.addGroup(group);

        double expResult = 0.5;
        double result = instance.getAtomCharge("CA3");
        assertTrue( Math.abs(result-expResult)<1e-15 );

        expResult = -0.2;
        result = instance.getAtomCharge("CA1");
        assertTrue( Math.abs(result-expResult)<1e-15 );

        expResult = 0.25;
        result = instance.getAtomCharge("CA4");
        assertTrue( Math.abs(result-expResult)<1e-15 );
        
    }

   

}