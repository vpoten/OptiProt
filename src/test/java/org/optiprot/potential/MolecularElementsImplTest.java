/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import org.biojava.bio.structure.Chain;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class MolecularElementsImplTest {

    static MolecularElementsImpl instance = null;

    public MolecularElementsImplTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        ConfigOptiProtTest.setName("MolecularElementsImplTest");
        ConfigOptiProtTest.loadRotamerLibrary();
    }
    
    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    @Test
    public void testConstructor() {
        System.out.println("constructor");

        Chain chain=ConfigOptiProtTest.getChain1();
        instance=new MolecularElementsImpl( chain, ConfigOptiProtTest.getParameters().getGrid(),
                ConfigOptiProtTest.getParameters() );

        assertTrue(instance!=null);
    }

    @Test
    public void testCalcEnergy() {
        System.out.println("calcEnergy");

        Chain chain=ConfigOptiProtTest.getChain1();
        Double energ =
                MolecularElementsImpl.calcEnergy(chain, ConfigOptiProtTest.getParameters().getGrid(),
                ConfigOptiProtTest.getParameters(),true,true,true);

        assertTrue( energ!=null );
    }

    @Test
    public void testCalcEnergy2() {
        System.out.println("calcEnergy2");

        Chain chain=ConfigOptiProtTest.getChain1();
        Double energ =
                MolecularElementsImpl.calcEnergy(chain, ConfigOptiProtTest.getParameters().getGrid(),
                ConfigOptiProtTest.getParameters(),true,false,true);

        assertTrue( energ!=null );
    }

    @Test
    public void testCalcEnergy3() {
        System.out.println("calcEnergy3");

        Chain chain=ConfigOptiProtTest.getChain1();
        Double energ =
                MolecularElementsImpl.calcEnergy(chain, ConfigOptiProtTest.getParameters().getGrid(),
                ConfigOptiProtTest.getParameters(),true,false,false);

        assertTrue( energ!=null );
    }

    @Test
    public void testCalcEnergy4() {
        System.out.println("calcEnergy4");

        Chain chain=ConfigOptiProtTest.getChain1();
        Double energ =
                MolecularElementsImpl.calcEnergy(chain, ConfigOptiProtTest.getParameters().getGrid(),
                ConfigOptiProtTest.getParameters(),false,true,false);

        assertTrue( energ!=null );
    }

    /**
     * Test of getChain method, of class MolecularElementsImpl.
     */
    @Test
    public void testGetChain() {
        System.out.println("getChain");

        Chain result = instance.getChain();
        assertTrue( result!=null );
    }

    /**
     * Test of hasAngles method, of class MolecularElementsImpl.
     */
    @Test
    public void testHasAngles() {
        System.out.println("hasAngles");

        boolean expResult = false;
        boolean result = instance.hasAngles();
        assertEquals(expResult, result);
        
    }

    /**
     * Test of hasBonds method, of class MolecularElementsImpl.
     */
    @Test
    public void testHasBonds() {
        System.out.println("hasBonds");

        boolean expResult = true;
        boolean result = instance.hasBonds();
        assertEquals(expResult, result);

    }

    /**
     * Test of hasDihedrals method, of class MolecularElementsImpl.
     */
    @Test
    public void testHasDihedrals() {
        System.out.println("hasDihedrals");

        boolean expResult = true;
        boolean result = instance.hasDihedrals();
        assertEquals(expResult, result);
    }

    /**
     * Test of hasImpropers method, of class MolecularElementsImpl.
     */
    @Test
    public void testHasImpropers() {
        System.out.println("hasImpropers");

        boolean expResult = true;
        boolean result = instance.hasImpropers();
        assertEquals(expResult, result);
    }

    /**
     * Test of hasNonbondeds method, of class MolecularElementsImpl.
     */
    @Test
    public void testHasNonbondeds() {
        System.out.println("hasNonbondeds");

        boolean expResult = true;
        boolean result = instance.hasNonbondeds();
        assertEquals(expResult, result);
    }

    /**
     * Test of hasPairs method, of class MolecularElementsImpl.
     */
    @Test
    public void testHasPairs() {
        System.out.println("hasPairs");

        boolean expResult = true;
        boolean result = instance.hasPairs();
        assertEquals(expResult, result);
    }


    @Test
    public void testFitnessOrder() {

        System.out.println("FitnessOrder");

        int chains=3;
        int step=3;

        double energ[]=new double [chains*step];

//        String names [] = new String [] {"romo1_itasser.pdb", "optiprot_01.pdb",
//                "optiprot_06.pdb" };
        String names [] = new String [] {"3TRX.pdb", "optiprot_3trx_01.pdb",
                "optiprot_3trx_02.pdb" };

        for(int i=0;i<chains;i++){

            Chain chain=ConfigOptiProtTest.readChain( names[i], new int [] {0});

            energ[i*step+0] = MolecularElementsImpl.calcEnergy(chain, ConfigOptiProtTest.getParameters().getGrid(),
                       ConfigOptiProtTest.getParameters(),true,false,false);

            energ[i*step+1] = MolecularElementsImpl.calcEnergy(chain, ConfigOptiProtTest.getParameters().getGrid(),
                       ConfigOptiProtTest.getParameters(),false,true,false);

            energ[i*step+2] = MolecularElementsImpl.calcEnergy(chain, ConfigOptiProtTest.getParameters().getGrid(),
                       ConfigOptiProtTest.getParameters(),false,false,true);

        }
       

        assertTrue( energ[0]<energ[step] );

    }
  

}