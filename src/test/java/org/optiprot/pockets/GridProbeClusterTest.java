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
public class GridProbeClusterTest {

    public GridProbeClusterTest() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

   
    /**
     * Test of isAdjoin method, of class GridProbeCluster.
     */
    @Test
    public void testIsAdjoin_GridProbe_AtomGrid() {
        System.out.println("isAdjoin");
        
        GridProbe probe = new GridProbe(9,null,0);
        AtomGrid grid = ConfigOptiProtTest.getParameters().getGrid();
        grid.setResol(1.0);
        grid.setDimension(10, 10, 10);

        GridProbeCluster instance = new GridProbeCluster();
        instance.getProbes().add( new GridProbe(0,null,0) );
        instance.getProbes().add( new GridProbe(10,null,0) );

        boolean expResult = false;
        boolean result = instance.isAdjoin(probe, grid);

        assertEquals(expResult, result);

        probe.setIndex(101);
        expResult = true;
        result = instance.isAdjoin(probe, grid);

        assertEquals(expResult, result);
    }

    /**
     * Test of isAdjoin method, of class GridProbeCluster.
     */
    @Test
    public void testIsAdjoin_GridProbeCluster_AtomGrid() {
        System.out.println("isAdjoin");

        AtomGrid grid = ConfigOptiProtTest.getParameters().getGrid();
        grid.setResol(1.0);
        grid.setDimension(10, 10, 10);

        GridProbeCluster cluster = new GridProbeCluster();
        cluster.getProbes().add( new GridProbe(1,null,0) );
        cluster.getProbes().add( new GridProbe(11,null,0) );

        GridProbeCluster instance = new GridProbeCluster();
        instance.getProbes().add( new GridProbe(0,null,0) );
        instance.getProbes().add( new GridProbe(10,null,0) );

        boolean expResult = true;
        boolean result = instance.isAdjoin(cluster, grid);

        assertEquals(expResult, result);

        cluster.getProbes().clear();
        cluster.getProbes().add( new GridProbe(9,null,0) );
        cluster.getProbes().add( new GridProbe(19,null,0) );

        expResult = false;
        result = instance.isAdjoin(cluster, grid);

        assertEquals(expResult, result);
    }

    /**
     * Test of join method, of class GridProbeCluster.
     */
    @Test
    public void testJoin() {
        System.out.println("join");

        GridProbeCluster cluster = new GridProbeCluster();
        cluster.getProbes().add( new GridProbe(1,null,0) );
        cluster.getProbes().add( new GridProbe(11,null,0) );

        GridProbeCluster instance = new GridProbeCluster();
        instance.getProbes().add( new GridProbe(0,null,0) );
        instance.getProbes().add( new GridProbe(10,null,0) );


        instance.join(cluster);

        assertTrue( instance.getProbes().size()==4 );
        assertTrue( cluster.getProbes().isEmpty() );
    }


}