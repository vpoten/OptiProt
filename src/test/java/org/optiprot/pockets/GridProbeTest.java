/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.pockets;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.maths.AtomGrid;

/**
 *
 * @author victor
 */
public class GridProbeTest {

    public GridProbeTest() {
    }

    
    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    
    /**
     * Test of incBuriedIdx method, of class GridProbe.
     */
    @Test
    public void testIncBuriedIdx() {
        System.out.println("incBuriedIdx");
        
        GridProbe instance = new GridProbe(0,null,0);
        
        instance.incBuriedIdx(1);
        assertEquals(1, instance.getBuriedIdx());

        instance.incBuriedIdx(2);
        assertEquals(3, instance.getBuriedIdx());
    }

    /**
     * Test of getNeighbors method, of class GridProbe.
     */
    @Test
    public void testGetNeighbors() {
        System.out.println("getNeighbors");

        AtomGrid grid = ConfigOptiProtTest.getParameters().getGrid();
        grid.setResol(1.0);
        grid.setDimension(10, 10, 10);

        int[] neighbors = new int [GridProbe.NUM_NEIGHBORS];

        GridProbe instance = new GridProbe(0,null,0);
        instance.getNeighbors(grid, neighbors);

        assertEquals( neighbors[0], -1);
        assertEquals( neighbors[8], -1);

        
        instance.setIndex(999);
        instance.getNeighbors(grid, neighbors);

        assertEquals( neighbors[25], -1);
        assertEquals( neighbors[16], -1);

        instance.setIndex( (int)Math.floor(Math.random()*1000) );
        instance.getNeighbors(grid, neighbors);

        //assertEquals( neighbors[13], instance.getIndex());
    }

}