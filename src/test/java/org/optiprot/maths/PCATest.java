/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import Jama.EigenvalueDecomposition;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.math.array.LinearAlgebra;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class PCATest {

    public PCATest() {
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
     * Test of transformPCA method, of class PCA.
     */
    @Test
    public void testTransformPCA() throws Exception {
        System.out.println("transformPCA");
        
        Chain chain = ConfigOptiProtTest.readChain("eif6_swiss.pdb", new int [] {0});
        CalcTransform.center(chain);

        double[][] expResult = null;
        double[][] result = PCA.transformPCA(chain);

        assertTrue( result!=null );
    }

    @Test
    public void testTransformPCA2() throws Exception {
        System.out.println("transformPCA2");

        Chain chain = ConfigOptiProtTest.readChain("eif6_swiss.pdb", new int [] {0});
        CalcTransform.center(chain);
        
        Chain chain2=(Chain) chain.clone();


        double[][] result = PCA.transformPCA(chain);

        double rmsd = CalcRmsd.rmsd( chain, chain2);

        assertTrue(rmsd<0.01);
    }
    
    /**
     * Test of inversePCA method, of class PCA.
     */
    @Test
    public void testInversePCA() throws StructureException {
        System.out.println("inversePCA");

        Chain chain = ConfigOptiProtTest.readChain("eif6_swiss.pdb", new int [] {0});
        CalcTransform.center(chain);

        List<Atom> list=BSPTree.getListAtoms(chain);
        Atom [] atomSet=list.toArray(new Atom [list.size()]);
        
        int natoms=30;

        Atom[] atoms = new Atom [natoms];
        Atom[] expResult = new Atom [natoms];

        // pick natoms randomly from the chain
        for( int i=0;i<natoms;i++ ){
            atoms[i] = atomSet[ (int)Math.floor(Math.random()*atomSet.length) ];
            expResult[i]=(Atom) atoms[i].clone();
        }

        double[][] pca = PCA.transformPCA(chain);

        Atom[] result = PCA.inversePCA(atoms, pca);


        for( int i=0;i<natoms;i++ ){
            assertTrue( Calc.getDistance(result[i], expResult[i])<1e-6 );
        }

    }

    
    /**
     * Test of covariance method, of class PCA.
     */
    @Test
    public void testCovariance() {
        System.out.println("covariance");

        // a test case with mean 0
        double[][] v = new double [][] {
            {0.69, 0.49},
            {-1.31, -1.21},
            {0.39, 0.99},
            {0.09, 0.29},
            {1.29, 1.09},
            {0.49, 0.79},
            {0.19, -0.31},
            {-0.81, -0.81},
            {-0.31, -0.31},
            {-0.71, -1.01}
        };


        double[][] expResult = new double [][] {
            {0.616555, 0.615444},
            {0.615444, 0.716555}
        };

        double[][] result = PCA.covariance(v);
       
        assertTrue( Math.abs(expResult[0][0]-result[0][0])<1e-5 );
        assertTrue( Math.abs(expResult[0][1]-result[0][1])<1e-5 );
        assertTrue( Math.abs(expResult[1][0]-result[1][0])<1e-5 );
        assertTrue( Math.abs(expResult[1][1]-result[1][1])<1e-5 );

        EigenvalueDecomposition e = LinearAlgebra.eigen(result);
		result = e.getV().transpose().getArray();

        expResult = new double [][] {
            {-0.735178, 0.677873},
            {-0.677873, -0.735178}
        };

        assertTrue( Math.abs(expResult[0][0]-result[0][0])<1e-5 );
        assertTrue( Math.abs(expResult[0][1]-result[0][1])<1e-5 );

    }

}