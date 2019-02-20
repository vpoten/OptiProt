/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.Iterator;
import org.biojava.bio.structure.*;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class CalcGeomTest {


    public CalcGeomTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        ConfigOptiProtTest.setName("CalcGeomTest");

        //rotate chain1 randomly

        double angle=Math.random()*Math.PI;
        Atom vaxis=new AtomImpl();
        vaxis.setX(0.5);
        vaxis.setY(Math.random()-0.5 );
        vaxis.setZ(Math.random()-0.5 );

        ConfigOptiProtTest.setChain2( (Chain) ConfigOptiProtTest.getChain1().clone());

        Atom center= CalcGeom.getCentroid(ConfigOptiProtTest.getChain1());

        CalcTransform.rotateSubChain( center, vaxis, angle,
            ConfigOptiProtTest.getChain1(), 0,
            ConfigOptiProtTest.getChain1().getAtomGroups().size() );
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of squareDistance method, of class CalcGeom.
     */
    @Test
    public void testSquareDistance() {
        System.out.println("squareDistance");
        
        Atom at1 = new AtomImpl();
        Atom at2 = new AtomImpl();

        at1.setCoords(new double [] {1,1,1});
        at2.setCoords(new double [] {-1,-1,-1});
        
        double expResult = 12.0;
        double result = CalcGeom.squareDistance(at1, at2);

        assertTrue( Math.abs(result-expResult)<1e-12);
    }

    

    /**
     * Test of getCentroid method, of class CalcRmsd.
     */
    @Test
    public void testGetCentroid1() throws Exception {
        System.out.println("Test Calc.getCentroid");

        Atom center1=CalcGeom.getCentroid(ConfigOptiProtTest.getChain1());
        Atom center2=CalcGeom.getCentroid(ConfigOptiProtTest.getChain2());

        Atom shiftVect1=new AtomImpl();
        shiftVect1.setCoords( new double [] {-center1.getX(),-center1.getY(),-center1.getZ()});
        Atom shiftVect2=new AtomImpl();
        shiftVect2.setCoords( new double [] {-center2.getX(),-center2.getY(),-center2.getZ()});

        // center chain1 & chain2
        Iterator<Group> it1=ConfigOptiProtTest.getChain1().getAtomGroups().iterator();
        Iterator<Group> it2=ConfigOptiProtTest.getChain2().getAtomGroups().iterator();

        while(it1.hasNext()){

            AminoAcid aa1=(AminoAcid) it1.next();
            AminoAcid aa2=(AminoAcid) it2.next();

            Calc.shift(aa1, shiftVect1);
            Calc.shift(aa2, shiftVect2);
        }

        center1=CalcGeom.getCentroid(ConfigOptiProtTest.getChain1());
        center2=CalcGeom.getCentroid(ConfigOptiProtTest.getChain2());

        double dist=Calc.getDistance(center1, center2);

        assertTrue( dist<0.001 );
    }



    /**
     * Test of rayPlaneIntersec method, of class CalcGeom.
     */
    @Test
    public void testRayPlaneIntersec() throws StructureException {
        System.out.println("rayPlaneIntersec");
        
        Atom org = new AtomImpl();
        Atom vector = new AtomImpl();
        Atom normal = new AtomImpl();
        Atom x0 = new AtomImpl();

        org.setCoords(new double [] {1,1,1});
        vector.setCoords(new double [] {-1,-1,-1});
        vector=Calc.unitVector(vector);
        normal.setCoords(new double [] {1,0,0});
        x0.setCoords(new double [] {0,-1,-1});
        
        Atom expResult = new AtomImpl();
        expResult.setCoords(new double [] {0,0,0});
        Atom result = CalcGeom.rayPlaneIntersec(org, vector, normal, x0);

        assertTrue( Calc.getDistance(result,expResult)<1e-12);
        
    }

}