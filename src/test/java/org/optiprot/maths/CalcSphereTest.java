/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.*;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class CalcSphereTest {

    public CalcSphereTest() {
    }


    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of volume method, of class CalcSphere.
     */
    @Test
    public void testVolume() {
        System.out.println("volume");
        double radius = 1.0;
        double expResult = 4.0*Math.PI/3.0;
        double result = CalcSphere.volume(radius);
        assertTrue( Math.abs(expResult-result)<1e-20 );
    }

    /**
     * Test of area method, of class CalcSphere.
     */
    @Test
    public void testArea() {
        System.out.println("area");
        double radius = 1.0;
        double expResult = 4.0*Math.PI;
        double result = CalcSphere.area(radius);
        assertTrue( Math.abs(expResult-result)<1e-20 );
    }

    /**
     * Test of intersecVol method, of class CalcSphere.
     */
    @Test
    public void testIntersecVol() throws StructureException {
        System.out.println("intersecVol");
        Atom c1 = new AtomImpl();
        c1.setCoords(new double [] {0,0,0} );
        double r1 = 1.0;
        Atom c2 = new AtomImpl();
        c2.setCoords(new double [] {0,0,2} );
        double r2 = 2.0;

        double dist=Calc.getDistance(c1, c2);
        double a=(1.0/(2*dist))*Math.sqrt( (-dist+r1-r2)*
                (-dist-r1+r2)*(-dist+r1+r2)*(dist+r1+r2) );
        double h1=r1-Math.sqrt(r1*r1-a*a);
        double h2=r2-Math.sqrt(r2*r2-a*a);

        double expResult = (Math.PI/6.0)*( h1*(3*a*a+h1*h1)+h2*(3*a*a+h2*h2) );
        double result = CalcSphere.intersecVol(c1, r1, c2, r2);

        assertTrue( Math.abs(expResult-result)<1e-15 );

        c2.setCoords(new double [] {0,0,3} );
        result = CalcSphere.intersecVol(c1, r1, c2, r2);
        assertTrue( Math.abs(result)<1e-15 );
    }

    /**
     * Test of intersecArea method, of class CalcSphere.
     */
    @Test
    public void testIntersecArea_4args() throws StructureException {
        System.out.println("intersecArea");
        Atom c1 = new AtomImpl();
        c1.setCoords(new double [] {0,0,0} );
        double r1 = 1.0;
        Atom c2 = new AtomImpl();
        c2.setCoords(new double [] {0,0,2} );
        double r2 = 2.0;

        double dist=Calc.getDistance(c1, c2);
        double a=(1.0/(2*dist))*Math.sqrt( (-dist+r1-r2)*
                (-dist-r1+r2)*(-dist+r1+r2)*(dist+r1+r2) );
        double h1=r1-Math.sqrt(r1*r1-a*a);
        double h2=r2-Math.sqrt(r2*r2-a*a);

        double expResult = 2.0*Math.PI*(r1*h1+r2*h2);
        double result = CalcSphere.intersecArea(c1, r1, c2, r2);
        
        assertTrue( Math.abs(expResult-result)<1e-12 );
    }

    /**
     * Test of sphericCapArea method, of class CalcSphere.
     */
    @Test
    public void testSphericCapArea() {
        System.out.println("sphericCapArea");
        double radius = 2.0;
        double h = 2.0;
        double expResult = 25.1327412287183459077;
        double result = CalcSphere.sphericCapArea(radius, h);

        assertTrue( Math.abs(expResult-result)<1e-20 );
    }

    /**
     * Test of sphericCapVol method, of class CalcSphere.
     */
    @Test
    public void testSphericCapVol() {
        System.out.println("sphericCapVol");
        double radius = 2.0;
        double h = 2.0;
        double expResult = 16.75516081914556393847;
        double result = CalcSphere.sphericCapVol(radius, h);

        assertTrue( Math.abs(expResult-result)<1e-12 );
    }

    
    /**
     * Test of sphericalExcess method, of class CalcSphere.
     */
    @Test
    public void testSphericalExcess() {
        System.out.println("sphericalExcess");
        Atom vA = new AtomImpl();
        vA.setCoords(new double [] {1,0,0} );
        Atom vB = new AtomImpl();
        vB.setCoords(new double [] {0,1,0} );
        Atom vC = new AtomImpl();
        vC.setCoords(new double [] {0,0,1} );

        double expResult = Math.PI*0.5;
        double result = CalcSphere.sphericalExcess(vA, vB, vC);

        assertTrue( Math.abs(expResult-result)<1e-12 );

        result = CalcSphere.sphericalExcess(vB, vA, vC);
        assertTrue( Math.abs(expResult-result)<1e-12 );
        
        vA.setCoords(new double [] {0,1,0} );
        vB.setCoords(new double [] {-1,0,0} );
        vB=Calc.unitVector(vB);
        vC.setCoords(new double [] {1,0,0} );
        vC=Calc.unitVector(vC);

        expResult = 2.0*Math.PI;
        result = CalcSphere.sphericalExcess(vA, vB, vC);

        assertTrue( Math.abs(expResult-result)<1e-12 );

        vA.setCoords(new double [] {0,0,1} );
        vB.setCoords(new double [] {0,1,0} );
        vC.setCoords(new double [] {0,1,0} );

        expResult = 0.0;
        result = CalcSphere.sphericalExcess(vA, vB, vC);

        assertTrue( Math.abs(expResult-result)<1e-12 );

    }

    /**
     * Test of intersecArea method, of class CalcSphere.
     */
    @Test
    public void testIntersecArea_6args() throws Exception {
        System.out.println("intersecArea");
        Atom c1 = new AtomImpl();
        c1.setCoords(new double [] {0,0,0} );
        double r1 = 1.0;
        Atom c2 = new AtomImpl();
        c2.setCoords(new double [] {2,0,0} );
        double r2 = 2.0;
        Atom c3 = new AtomImpl();
        c3.setCoords(new double [] {0,2.5,0} );
        double r3 = 1.5;
        double expResult = 0.0;
        double result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result)<1e-12 );


        //topological situation 1 : intersection without contact
        c1.setCoords(new double [] {0,0,0} );
        c2.setCoords(new double [] {0,0,0} );
        c3.setCoords(new double [] {0,0,0} );

        expResult = CalcSphere.area(r1);
        result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result)<1e-12 );

        //empty simplicia1 topology
        c1.setCoords(new double [] {0,0,0} );
        c2.setCoords(new double [] {1.8,0,0} );
        c3.setCoords(new double [] {0.9,1.6,0} );
        r1=1;
        r2=1;
        r3=1;

        expResult = 0;
        result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result)<1e-12 );

        //general situation
        c1.setCoords(new double [] {1,1,1} );
        c2.setCoords(new double [] {-1,-1,1} );
        c3.setCoords(new double [] {-1,1,-1} );
        r1=2.0*Math.sqrt(2.0)+1.4;
        r2=2.0*Math.sqrt(2.0)+1.4;
        r3=2.0*Math.sqrt(2.0)+1.4;

        expResult = CalcSphere.monteCarloSurf(c1, r1, c2, r2, c3, r3, 100000);
        result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result) < expResult*0.01 );
        result = CalcSphere.intersecArea(c1, r1, c3, r3, c2, r2, 4);
        assertTrue( Math.abs(expResult-result) < expResult*0.01 );

        //if a sphere i is inside a sphere k and the intersection i,j is inside k
        c1.setCoords(new double [] {0,0,0} );
        c2.setCoords(new double [] {2,0,0} );
        c3.setCoords(new double [] {2,0,0} );
        r1=2;
        r2=1;
        r3=2;

        expResult = CalcSphere.intersecArea(c1,r1,c2,r2);
        result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result)<1e-12 );
        result = CalcSphere.intersecArea(c1, r1, c3, r3, c2, r2, 4);
        assertTrue( Math.abs(expResult-result)<1e-12 );

        //if an intersection of order 2 (i,j) is inside the remaining sphere (k)
        // and the union of i,j is inside k
        c1.setCoords(new double [] {0,0,0} );
        c2.setCoords(new double [] {-0.6,0,0} );
        c3.setCoords(new double [] {0.6,0,0} );
        r1=3;
        r2=1;
        r3=1;

        expResult = CalcSphere.intersecArea(c3,r3,c2,r2);
        result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result)<1e-12 );
        result = CalcSphere.intersecArea(c2, r2, c1, r1, c3, r3, 4);
        assertTrue( Math.abs(expResult-result)<1e-12 );

        //if the intersection i,j is inside k and k intersects with i and j
        c1.setCoords(new double [] {0,0,0} );
        c2.setCoords(new double [] {-0.6,0,0} );
        c3.setCoords(new double [] {0.6,0,0} );
        r1=1.1;
        r2=1;
        r3=1;

        expResult = CalcSphere.intersecArea(c3,r3,c2,r2);
        result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result)<1e-12 );
        result = CalcSphere.intersecArea(c2, r2, c1, r1, c3, r3, 4);
        assertTrue( Math.abs(expResult-result)<1e-12 );

        //if a sphere k intersects the i,j intersection where i,j intersection
        // is not inside k
        c1.setCoords(new double [] {0,0,0} );
        c2.setCoords(new double [] {-3.7,0,0} );
        c3.setCoords(new double [] {3.6,0,0} );
        r1=1;
        r2=4;
        r3=4;

        expResult = CalcSphere.intersecArea(c1,r1,c2,r2) +
                CalcSphere.intersecArea(c3,r3,c1,r1) - CalcSphere.area(r1);
        result = CalcSphere.intersecArea(c1, r1, c2, r2, c3, r3, 4);

        assertTrue( Math.abs(expResult-result)<1e-12 );
        result = CalcSphere.intersecArea(c2, r2, c1, r1, c3, r3, 4);
        assertTrue( Math.abs(expResult-result)<1e-12 );

    }

}