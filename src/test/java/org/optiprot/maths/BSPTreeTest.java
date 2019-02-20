/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.ArrayList;
import java.util.List;
import org.biojava.bio.structure.*;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.potential.IForceField;
import org.optiprot.potential.TestForceField;

/**
 *
 * @author victor
 */
public class BSPTreeTest {

    static BSPTree instance = new BSPTree( ConfigOptiProtTest.getChain1() );

    public BSPTreeTest() {
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
     * Test of boundingBox method, of class BSPTree.
     */
    @Test
    public void testBoundingBox_Chain() throws StructureException {
        System.out.println("boundingBox");
        Chain chain = new ChainImpl();
        Group group=new AminoAcidImpl();
        chain.addGroup(group);

        Atom at=new AtomImpl();
        at.setCoords(new double [] {0.2,-0.7,1});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {-1,1,0.56});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {0,-1,0});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {1,-0.9,-1});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {0,-0.97,-1});
        group.addAtom(at);

        Atom[] expResult = new Atom [2];
        expResult[0]=new AtomImpl();
        expResult[1]=new AtomImpl();
        expResult[0].setCoords(new double [] {-1,-1,-1});
        expResult[1].setCoords(new double [] {1,1,1});
        Atom[] result = BSPTree.boundingBox(chain);
        assertTrue( Math.abs(Calc.getDistance(result[0],expResult[0]))<1e-12);
        assertTrue( Math.abs(Calc.getDistance(result[1],expResult[1]))<1e-12);
    }

    /**
     * Test of getBBox method, of class BSPTree.
     */
    @Test
    public void testGetBBox() throws StructureException {
        System.out.println("getBBox");
        
        Chain chain = new ChainImpl();
        Group group=new AminoAcidImpl();
        chain.addGroup(group);

        Atom at=new AtomImpl();
        at.setCoords(new double [] {0.2,-0.7,1});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {-1,1,0.56});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {0,-1,0});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {1,-0.9,-1});
        group.addAtom(at);

        at=new AtomImpl();
        at.setCoords(new double [] {0,-0.97,-1});
        group.addAtom(at);

        BSPTree linstance = new BSPTree(chain);

        Atom[] expResult = new Atom [2];
        expResult[0]=new AtomImpl();
        expResult[1]=new AtomImpl();
        expResult[0].setCoords(new double [] {-1,-1,-1});
        expResult[1].setCoords(new double [] {1,1,1});
        Atom[] result = linstance.getBBox();
        assertTrue( Math.abs(Calc.getDistance(result[0],expResult[0]))<1e-12);
        assertTrue( Math.abs(Calc.getDistance(result[1],expResult[1]))<1e-12);
    }

    /**
     * Test of isInsideVolume method, of class BSPTree.
     */
    @Test
    public void testIsInsideVolume() {
        System.out.println("isInsideVolume");
        Atom atom = new AtomImpl();
        IForceField ffield = new TestForceField();
        double p_radius = 0.0;

        Atom[] bbox = instance.getBBox();

        atom.setCoords(new double [] {79,12,12});
        boolean expResult = false;
        boolean result = instance.isInsideVolume(atom, ffield, p_radius);
        assertEquals(expResult, result);

        atom.setCoords(new double [] {65, 41 ,-3});
        expResult = true;
        result = instance.isInsideVolume(atom, ffield, p_radius);
        assertEquals(expResult, result);

        atom.setCoords(new double [] {65, bbox[0].getY()-p_radius-0.1 ,-3});
        expResult = false;
        result = instance.isInsideVolume(atom, ffield, p_radius);
        assertEquals(expResult, result);

        atom.setCoords(new double [] {39, 35, 24.5});
        expResult = true;
        result = instance.isInsideVolume(atom, ffield, p_radius);
        assertEquals(expResult, result);
        
    }

    /**
     * Test of neighbours method, of class BSPTree.
     */
    @Test
    public void testNeighbours() {
        System.out.println("neighbours");
        Atom atom = new AtomImpl();
        
        List<Atom> list = new ArrayList<Atom>();

        //for eif6_swiss PDB

        double distance = 1.55;// a CB atom
        atom.setCoords(new double [] {44.440,38.782,19.880});
        instance.neighbours(atom, distance, list);
        int expResult=3;
        assertEquals(expResult, list.size());

        list.clear();
        distance = 1.40;// a CZ atom of TYR
        atom.setCoords(new double [] {42.319, 33.031, 6.454});
        instance.neighbours(atom, distance, list);
        expResult=4;
        assertEquals(expResult, list.size());

        list.clear();
        distance = 1.20;// a C atom of GLY
        atom.setCoords(new double [] {54.829,32.116,9.599});
        instance.neighbours(atom, distance, list);
        expResult=1;
        assertEquals(expResult, list.size());
    }

    /**
     * Test of nearestNeighbour method, of class BSPTree.
     */
    @Test
    public void testNearestNeighbour() throws StructureException {
        System.out.println("nearestNeighbour");
        Atom point = new AtomImpl();
        point.setCoords(new double [] {54.7,32.05,9.4});
        
        Atom expResult = new AtomImpl();
        expResult.setCoords(new double [] {54.829,32.116,9.599});
        Atom result = instance.nearestNeighbour(point);
        assertTrue( Math.abs(Calc.getDistance(result,expResult) )<1e-12);
    }

}