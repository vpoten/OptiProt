/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.cgal;

import org.biojava.bio.structure.*;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.maths.CalcTransform;
import org.optiprot.potential.IForceField;
import org.optiprot.potential.TestForceField;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class AlphaShapeTest {

    public AlphaShapeTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        ConfigOptiProtTest.setName("AlphaShapeTest");
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of getSimplices method, of class AlphaShape.
     */
    @Test
    public void testGetSimplices() {
        System.out.println("getSimplices");

        Chain chain1=ConfigOptiProtTest.getChain1();

        Atom [] molecule=CalcTransform.toAtomArray(chain1);
        double [] radii=new double [molecule.length];

        for(int i=0;i<radii.length;i++){
            radii[i]=1.6;
        }

        AlphaShape instance = new AlphaShape(molecule, radii, 1.4);

        //get the simplices
        BaseSimplex [] simplices = instance.getSimplices();

        testSimplices(simplices);

        assertTrue(simplices!=null);
    }

    @Test
    public void testFigureProbe() {
        System.out.println("testFigureProbe");

        int natom=15;

        Atom [] molecule=new Atom [natom];
        double [] radii=new double [natom];

        molecule[0]=new AtomImpl();
        molecule[0].setCoords(new double [] {0,0,0});
        molecule[1]=new AtomImpl();
        molecule[1].setCoords(new double [] {4,0,0});
        molecule[2]=new AtomImpl();
        molecule[2].setCoords(new double [] {8,0,0});
        molecule[3]=new AtomImpl();
        molecule[3].setCoords(new double [] {0,0,4});
        molecule[4]=new AtomImpl();
        molecule[4].setCoords(new double [] {4,0,4});
        molecule[5]=new AtomImpl();
        molecule[5].setCoords(new double [] {8,0,4});
        molecule[6]=new AtomImpl();
        molecule[6].setCoords(new double [] {0,4,0});
        molecule[7]=new AtomImpl();
        molecule[7].setCoords(new double [] {4,4,0});
        molecule[8]=new AtomImpl();
        molecule[8].setCoords(new double [] {8,4,0});
        molecule[9]=new AtomImpl();
        molecule[9].setCoords(new double [] {0,4,4});
        molecule[10]=new AtomImpl();
        molecule[10].setCoords(new double [] {4,4,4});
        molecule[11]=new AtomImpl();
        molecule[11].setCoords(new double [] {8,4,4});
        molecule[12]=new AtomImpl();
        molecule[12].setCoords(new double [] {2,6,2});
        molecule[13]=new AtomImpl();
        molecule[13].setCoords(new double [] {6,6,2});
        molecule[14]=new AtomImpl();
        molecule[14].setCoords(new double [] {10,2,2});
//        molecule[15]=new AtomImpl();
//        molecule[15].setCoords(new double [] {6,2,2});

        for(int i=0;i<natom;i++){
            radii[i]=3.9;
        }

        radii[12]=4.2;
        radii[13]=4.2;
        radii[14]=4.2;

        AlphaShape instance = new AlphaShape(molecule, radii, 2.8);

        //get the simplices
        BaseSimplex [] simplices = instance.getSimplices();
        //printSimplices(simplices);
        
        testSimplices(simplices);

    }

    private void printSimplices(BaseSimplex [] simplices){
        for(BaseSimplex simplex : simplices){
            if( simplex==null)
                break;
            
            simplex.printSimplex();
        }
    }

    private void testSimplices(BaseSimplex [] simplices){

        for(BaseSimplex simplex : simplices){

            if( simplex==null)
                break;

            if( simplex.getOutAngle()>1.0){
                fail( "Out angle >1");
            }

            if( simplex.getOutAngle()<0.0){
                fail( "Out angle <0");
            }
        }

    }

    @Test
    public void testCalcAreaVolume(){

        System.out.println("CalcAreaVolume");

        Chain chain1=ConfigOptiProtTest.getChain1();
        Atom [] molecule=CalcTransform.toAtomArray(chain1);
        IForceField ff = new TestForceField();

        double [] radii=new double [molecule.length];

        for(int i=0;i<radii.length;i++){
            radii[i]=ff.getVDWRadius(molecule[i]);
        }

        AlphaShape instance = new AlphaShape(molecule, radii, 1.4);

        double result=0.0;
        double expResult=0.0;
        double error=0.0;

        expResult = 9561.13346695531;//for eif6_swiss
        result = instance.getArea(4);
        error=expResult*0.15;
        assertTrue( Math.abs(result-expResult)<error);

        
        chain1=ConfigOptiProtTest.readChain("1CAU.pdb", new int [] {0,1} );
        molecule=CalcTransform.toAtomArray(chain1);

        radii=new double [molecule.length];

        for(int i=0;i<radii.length;i++){
            radii[i]=ff.getVDWRadius(molecule[i]);
        }

        instance = new AlphaShape(molecule, radii, 0);
        expResult = 37453.0;//1CAU
        result = instance.getArea(4);
        error=expResult*0.1;
        assertTrue( Math.abs(result-expResult)<error);


        chain1=ConfigOptiProtTest.readChain("1ECA.pdb", new int [] {0} );
        molecule=CalcTransform.toAtomArray(chain1);

        radii=new double [molecule.length];

        for(int i=0;i<radii.length;i++){
            radii[i]=ff.getVDWRadius(molecule[i]);
        }

        instance = new AlphaShape(molecule, radii, 1.4);
        expResult = 25000.0;//1ECA
        result = instance.getVolume(50000);
        error=expResult*0.1;
        assertTrue( Math.abs(result-expResult)<error);
        
    }


}