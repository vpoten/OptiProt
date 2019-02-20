/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.pockets;

import java.util.ArrayList;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.PDBParseException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import org.optiprot.jmol.JmolPanel;
import static org.junit.Assert.*;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.BSPTree;

/**
 *
 * @author victor
 */
public class CalcPocketsTest {

    public CalcPocketsTest() {
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
     * Test of calc method, of class CalcPockets.
     */
    @Test
    public void testCalc() throws Exception {
        System.out.println("calc");

        //Chain chain = ConfigOptiProtTest.readChainPDB("3tms", new int [] {0});
        Chain chain = ConfigOptiProtTest.readChainPDB("1bid", new int [] {0});

        //UMP coords
        HetatomImpl ligand =
                ConfigOptiProtTest.readHetatomPDB("1bid", "UMP", 0);
 
        ////CalcTransform.center(chain);
        ////double [][] pca = PCA.transformPCA( chain );

        AtomGrid grid = ConfigOptiProtTest.getParameters().getGrid();
        int npockets = 3;


        // calc pockets
        GridProbeCluster[] result = CalcPockets.calc(chain, grid, npockets);

        //drawChainPockets( chain , result ) ;

        //check the distance between the pocket center and the ligand
        assertTrue( getMinDistance(result[0],ligand)<4.0 );


        // test 2ifb
        chain = ConfigOptiProtTest.readChainPDB("2ifb", new int [] {0});

        ligand = ConfigOptiProtTest.readHetatomPDB("2ifb", "PLM", 0);

        result = CalcPockets.calc(chain, grid, npockets);

        //check the distance between the pocket center and the ligand
        assertTrue( getMinDistance(result[0],ligand)<4.0 );


    }


    static public void drawChainPockets( Chain chain , GridProbeCluster[] pockets )
            throws PDBParseException {

        //draw the chain and the pockets (jmol)
        JmolPanel panel=new JmolPanel();
        panel.addChain(chain);

        for( GridProbeCluster cluster : pockets){
            HetatomImpl het=new HetatomImpl();
            ArrayList<Atom> list=new ArrayList<Atom>();
            cluster.getAtoms(list);

            for( Atom at : list ){
                het.addAtom(at);
            }
            ///het.addAtom( center );
            het.setPDBCode("");
            het.setPDBName("");
            panel.addHetatm(0, het);
        }


        panel.setStruct();

        panel.evalCommand("select * ; colour chain;");
        panel.evalCommand("select hetero ; colour red;");
        panel.evalCommand("select protein ; isosurface ignore(hetero) molecular");
        ///panel.evalCommand("select *; spacefill off; wireframe off; backbone 0.4;");

    }
    
    private double getMinDistance( GridProbeCluster cluster, HetatomImpl ligand )
            throws StructureException{

        Atom center=cluster.getCentroid();

        double minDist=java.lang.Double.MAX_VALUE;

        for( Atom at : ligand.getAtoms() ){
            double dist=Calc.getDistance(at, center);

            if( dist<minDist ){
                minDist=dist;
            }
        }

        return minDist;
    }

    /**
     * Test of calcActiveSite method, of class CalcPockets.
     */
    @Test
    public void testCalcActiveSite() throws Exception {
        System.out.println("calcActiveSite");

        String pdbcode="1bid";

        Chain chain = ConfigOptiProtTest.readChainPDB(pdbcode, new int [] {0});

        //UMP coords
        HetatomImpl ligand =
                ConfigOptiProtTest.readHetatomPDB(pdbcode, "UMP", 0);

        ///chain = ConfigOptiProtTest.readChainPDB("2ifb", new int [] {0});
        ///ligand = ConfigOptiProtTest.readHetatomPDB("2ifb", "PLM", 0);

        Atom actSiteCenter =
                Calc.getCentroid( ligand.getAtoms().toArray(new Atom [ligand.getAtoms().size()]));

        AtomGrid grid = ConfigOptiProtTest.getParameters().getGrid();

        double side_length=18;
        BSPTree btree=new BSPTree( chain );

        GridProbeCluster result = CalcPockets.calcActiveSite(btree, actSiteCenter, side_length, grid);

        //drawChainPockets( chain , new GridProbeCluster [] {result} );

        assertTrue( Calc.getDistance(actSiteCenter, result.getCentroid())<4.0 );
    }

    @Test
    public void testGenerateActiveSite() throws Exception {
        System.out.println("generateActiveSite");

        String pdbcode="1bid";

        Chain chain = ConfigOptiProtTest.readChainPDB(pdbcode, new int [] {0});

        //UMP coords
        HetatomImpl ligand =
                ConfigOptiProtTest.readHetatomPDB(pdbcode, "UMP", 0);

        ///chain = ConfigOptiProtTest.readChainPDB("2ifb", new int [] {0});
        ///ligand = ConfigOptiProtTest.readHetatomPDB("2ifb", "PLM", 0);

        Atom actSiteCenter =
                Calc.getCentroid( ligand.getAtoms().toArray(new Atom [ligand.getAtoms().size()]));

        AtomGrid grid = ConfigOptiProtTest.getParameters().getGrid();


        BSPTree btree=new BSPTree( chain );
        double radii=5;

        GridProbeCluster result = CalcPockets.generateActiveSite(btree, actSiteCenter, radii, grid);

        //drawChainPockets( chain , new GridProbeCluster [] {result} );

        assertTrue( Calc.getDistance(actSiteCenter, result.getCentroid())<4 );
    }

}