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
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class MolecularPairTest {

    public MolecularPairTest() {
    }


    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

   

    /**
     * Test of getDistance method, of class MolecularPair.
     */
    @Test
    public void testGetDistance() {
        System.out.println("getDistance");

        Atom at1=new AtomImpl();
        at1.setCoords(new double [] {0,0,0});
        at1.setName("CA1");

        Atom at2=new AtomImpl();
        at2.setCoords(new double [] {1,0,1});
        at2.setName("CA2");

        MolecularPair instance = new MolecularPair(at1,at2);
        double expResult = Math.sqrt(2);
        double result = instance.getDistance();
        assertTrue( Math.abs(result-expResult)<1e-9 );

        instance.setAtomA(new double [] {-1,0,-1} );
        instance.setAtomB(new double [] {1,0,1} );
        expResult = 2*Math.sqrt(2);
        result = instance.getDistance();
        assertTrue( Math.abs(result-expResult)<1e-9 );

    }

    /**
     * Test of getSqrDistance method, of class MolecularPair.
     */
    @Test
    public void testGetSqrDistance() {
        System.out.println("getSqrDistance");

        Atom at1=new AtomImpl();
        at1.setCoords(new double [] {0,0,0});
        at1.setName("CA1");

        Atom at2=new AtomImpl();
        at2.setCoords(new double [] {1,0,1});
        at2.setName("CA2");

        MolecularPair instance = new MolecularPair(at1,at2);
        double expResult = 2;
        double result = instance.getSqrDistance();
        assertTrue( Math.abs(result-expResult)<1e-9 );

        instance.setAtomA(new double [] {-1,0,-1} );
        instance.setAtomB(new double [] {1,0,1} );
        expResult = 8;
        result = instance.getSqrDistance();
        assertTrue( Math.abs(result-expResult)<1e-9 );
    }

   

}