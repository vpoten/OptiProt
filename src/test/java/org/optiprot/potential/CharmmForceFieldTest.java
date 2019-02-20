/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.optiprot.configtests.ConfigOptiProtTest;
import static org.junit.Assert.*;
import org.optiprot.potential.element.MolecularAngle;
import org.optiprot.potential.element.MolecularBond;
import org.optiprot.potential.element.MolecularDihedral;
import org.optiprot.potential.element.MolecularImproper;
import org.optiprot.potential.element.MolecularNonbonded;
import org.optiprot.potential.element.MolecularPair;

/**
 *
 * @author victor
 */
public class CharmmForceFieldTest {

    static CharmmForceField instance =
            (CharmmForceField) ConfigOptiProtTest.getParameters().getForceField();

    public CharmmForceFieldTest() {
    }

    
    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of getKBond method, of class CharmmForceField.
     */
    @Test
    public void testGetKBond() {
        System.out.println("getKBond");
        MolecularBond bond = new MolecularBond();

        bond.setAtomTypeA( instance.getAtomType("CP2"));
        bond.setAtomTypeB( instance.getAtomType("CP3"));
        double expResult = 222.50;
        double result = instance.getKBond(bond);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        bond.setAtomTypeA(instance.getAtomType("HB"));
        bond.setAtomTypeB(instance.getAtomType("CT3"));
        expResult = 330;
        result = instance.getKBond(bond);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        bond.setAtomTypeA(instance.getAtomType("SS"));
        bond.setAtomTypeB(instance.getAtomType("CS"));
        expResult = 205;
        result = instance.getKBond(bond);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getEqDistBond method, of class CharmmForceField.
     */
    @Test
    public void testGetEqDistBond() {
        System.out.println("getEqDistBond");
        MolecularBond bond = new MolecularBond();

        bond.setAtomTypeA(instance.getAtomType("CA"));
        bond.setAtomTypeB(instance.getAtomType("CT3"));
        double expResult = 1.49;
        double result = instance.getEqDistBond(bond);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        bond.setAtomTypeA(instance.getAtomType("HA"));
        bond.setAtomTypeB(instance.getAtomType("CPM"));
        expResult = 1.09;
        result = instance.getEqDistBond(bond);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getKUreyBradley method, of class CharmmForceField.
     */
    @Test
    public void testGetKUreyBradley() {
        System.out.println("getKUreyBradley");
        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();

        MolecularAngle angle = new MolecularAngle();
        at1.setName("H1");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("HA"));
        at2.setName("C");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("CP2"));
        at3.setName("H2");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("HA"));
        double expResult = 5.40;
        double result = instance.getKUreyBradley(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        angle = new MolecularAngle();
        at1.setName("C1");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("CC"));
        at2.setName("C2");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("CT2"));
        at3.setName("H");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("HA"));
        expResult = 30;
        result = instance.getKUreyBradley(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getEqDistUreyBradley method, of class CharmmForceField.
     */
    @Test
    public void testGetEqDistUreyBradley() {
        System.out.println("getEqDistUreyBradley");

        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();

        MolecularAngle angle = new MolecularAngle();
        at1.setName("H1");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("HA"));
        at2.setName("C");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("CT1"));
        at3.setName("C2");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("CT2"));
        double expResult = 2.179;
        double result = instance.getEqDistUreyBradley(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        angle = new MolecularAngle();
        at1.setName("H");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("HA"));
        at2.setName("C1");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("CPM"));
        at3.setName("C2");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("CPA"));
        expResult = 0;
        result = instance.getEqDistUreyBradley(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getKBondAngle method, of class CharmmForceField.
     */
    @Test
    public void testGetKBondAngle() {
        System.out.println("getKBondAngle");

        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();

        MolecularAngle angle = new MolecularAngle();
        at1.setName("H1");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("H"));
        at2.setName("N");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("NY"));
        at3.setName("C2");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("CA"));
        double expResult = 28;
        double result = instance.getKBondAngle(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        angle = new MolecularAngle();
        at1.setName("H");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("HA"));
        at2.setName("C1");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("CPM"));
        at3.setName("C2");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("CPA"));
        expResult = 12.7;
        result = instance.getKBondAngle(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getEqDistBondAngle method, of class CharmmForceField.
     */
    @Test
    public void testGetEqDistBondAngle() {
        System.out.println("getEqDistBondAngle");

        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();

        MolecularAngle angle = new MolecularAngle();
        at1.setName("H1");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("H"));
        at2.setName("N");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("NY"));
        at3.setName("C2");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("CA"));
        double expResult = 126;
        double result = instance.getEqDistBondAngle(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        angle = new MolecularAngle();
        at1.setName("H");
        angle.setAtom( 0, at1, 0.0, instance.getAtomType("HB"));
        at2.setName("C1");
        angle.setAtom( 1, at2, 0.0, instance.getAtomType("CT1"));
        at3.setName("C2");
        angle.setAtom( 2, at3, 0.0, instance.getAtomType("C"));
        expResult = 109.5;
        result = instance.getEqDistBondAngle(angle);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getKDihedral method, of class CharmmForceField.
     */
    @Test
    public void testGetKDihedral() {
        System.out.println("getKDihedral");
        
        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();
        Atom at4=new AtomImpl();

        MolecularDihedral dihedral = new MolecularDihedral();
        at1.setName("H");
        dihedral.setAtom( 0, at1, 0.0, instance.getAtomType("H"));
        at2.setName("N");
        dihedral.setAtom( 1, at2, 0.0, instance.getAtomType("NH2"));
        at3.setName("C2");
        dihedral.setAtom( 2, at3, 0.0, instance.getAtomType("CC"));
        at4.setName("C3");
        dihedral.setAtom( 3, at4, 0.0, instance.getAtomType("CT2"));
        double expResult = 1.40;
        double result = instance.getKDihedral(dihedral);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        dihedral = new MolecularDihedral();
        at1.setName("H");
        dihedral.setAtom( 0, at1, 0.0, instance.getAtomType("SH"));
        at2.setName("N");
        dihedral.setAtom( 1, at2, 0.0, instance.getAtomType("CC"));
        at3.setName("C2");
        dihedral.setAtom( 2, at3, 0.0, instance.getAtomType("CT1"));
        at4.setName("C3");
        dihedral.setAtom( 3, at4, 0.0, instance.getAtomType("ZN"));
        expResult = 0.05;
        result = instance.getKDihedral(dihedral);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getNDihedral method, of class CharmmForceField.
     */
    @Test
    public void testGetNDihedral() {
        System.out.println("getNDihedral");

        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();
        Atom at4=new AtomImpl();

        MolecularDihedral dihedral = new MolecularDihedral();
        at1.setName("C11");
        dihedral.setAtom( 0, at1, 0.0, instance.getAtomType("CT3"));
        at2.setName("C12");
        dihedral.setAtom( 1, at2, 0.0, instance.getAtomType("CT2"));
        at3.setName("C2");
        dihedral.setAtom( 2, at3, 0.0, instance.getAtomType("CY"));
        at4.setName("C3");
        dihedral.setAtom( 3, at4, 0.0, instance.getAtomType("CPT"));
        double expResult = 2;
        double result = instance.getNDihedral(dihedral);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        dihedral = new MolecularDihedral();
        at1.setName("SH");
        dihedral.setAtom( 0, at1, 0.0, instance.getAtomType("SH"));
        at2.setName("C1");
        dihedral.setAtom( 1, at2, 0.0, instance.getAtomType("CT2"));
        at3.setName("C2");
        dihedral.setAtom( 2, at3, 0.0, instance.getAtomType("CD"));
        at4.setName("Z3");
        dihedral.setAtom( 3, at4, 0.0, instance.getAtomType("ZN"));
        expResult = 6;
        result = instance.getNDihedral(dihedral);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getDDihedral method, of class CharmmForceField.
     */
    @Test
    public void testGetDDihedral() {
        System.out.println("getDDihedral");

        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();
        Atom at4=new AtomImpl();

        MolecularDihedral dihedral = new MolecularDihedral();
        at1.setName("H");
        dihedral.setAtom( 0, at1, 0.0, instance.getAtomType("H"));
        at2.setName("N");
        dihedral.setAtom( 1, at2, 0.0, instance.getAtomType("NH1"));
        at3.setName("C2");
        dihedral.setAtom( 2, at3, 0.0, instance.getAtomType("C"));
        at4.setName("C3");
        dihedral.setAtom( 3, at4, 0.0, instance.getAtomType("CP1"));
        double expResult = 180;
        double result = instance.getDDihedral(dihedral);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        dihedral = new MolecularDihedral();
        at1.setName("SH");
        dihedral.setAtom( 0, at1, 0.0, instance.getAtomType("SH"));
        at2.setName("C1");
        dihedral.setAtom( 1, at2, 0.0, instance.getAtomType("CPB"));
        at3.setName("C2");
        dihedral.setAtom( 2, at3, 0.0, instance.getAtomType("CT3"));
        at4.setName("C3");
        dihedral.setAtom( 3, at4, 0.0, instance.getAtomType("ZN"));
        expResult = 0;
        result = instance.getDDihedral(dihedral);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getKTorsionAngle method, of class CharmmForceField.
     */
    @Test
    public void testGetKTorsionAngle() {
        System.out.println("getKTorsionAngle");

        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();
        Atom at4=new AtomImpl();

        MolecularImproper improper = new MolecularImproper();
        at1.setName("C1");
        improper.setAtom( 0, at1, 0.0, instance.getAtomType("CPB"));
        at2.setName("S");
        improper.setAtom( 1, at2, 0.0, instance.getAtomType("SH"));
        at3.setName("M");
        improper.setAtom( 2, at3, 0.0, instance.getAtomType("MG"));
        at4.setName("C3");
        improper.setAtom( 3, at4, 0.0, instance.getAtomType("CE1"));
        double expResult = 90;
        double result = instance.getKTorsionAngle(improper);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        improper = new MolecularImproper();
        at1.setName("N");
        improper.setAtom( 0, at1, 0.0, instance.getAtomType("NR3"));
        at2.setName("C1");
        improper.setAtom( 1, at2, 0.0, instance.getAtomType("CPH1"));
        at3.setName("C2");
        improper.setAtom( 2, at3, 0.0, instance.getAtomType("CPH2"));
        at4.setName("H3");
        improper.setAtom( 3, at4, 0.0, instance.getAtomType("H"));
        expResult = 1.20;
        result = instance.getKTorsionAngle(improper);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getEqDistTorsionAngle method, of class CharmmForceField.
     */
    @Test
    public void testGetEqDistTorsionAngle() {
        System.out.println("getEqDistTorsionAngle");

        Atom at1=new AtomImpl();
        Atom at2=new AtomImpl();
        Atom at3=new AtomImpl();
        Atom at4=new AtomImpl();

        MolecularImproper improper = new MolecularImproper();
        at1.setName("C1");
        improper.setAtom( 0, at1, 0.0, instance.getAtomType("CPB"));
        at2.setName("S");
        improper.setAtom( 1, at2, 0.0, instance.getAtomType("SH"));
        at3.setName("M");
        improper.setAtom( 2, at3, 0.0, instance.getAtomType("MG"));
        at4.setName("C3");
        improper.setAtom( 3, at4, 0.0, instance.getAtomType("CE1"));
        double expResult = 0;
        double result = instance.getEqDistTorsionAngle(improper);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        improper = new MolecularImproper();
        at1.setName("O");
        improper.setAtom( 0, at1, 0.0, instance.getAtomType("O"));
        at2.setName("C1");
        improper.setAtom( 1, at2, 0.0, instance.getAtomType("CT3"));
        at3.setName("N2");
        improper.setAtom( 2, at3, 0.0, instance.getAtomType("NH2"));
        at4.setName("C");
        improper.setAtom( 3, at4, 0.0, instance.getAtomType("CC"));
        expResult = 0;
        result = instance.getEqDistTorsionAngle(improper);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getVDWRmin method, of class CharmmForceField.
     */
    @Test
    public void testGetVDWRmin() {
        System.out.println("getVDWRmin");
        
        MolecularNonbonded nonbond = new MolecularNonbonded();
        nonbond.setAtomTypeA(instance.getAtomType("CA"));
        nonbond.setAtomTypeB(instance.getAtomType("CE1"));
        double expResult = 1.992400+2.09;
        double result = instance.getVDWRmin(nonbond);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        nonbond = new MolecularNonbonded();
        nonbond.setAtomTypeA(instance.getAtomType("NY"));
        nonbond.setAtomTypeB(instance.getAtomType("HT"));
        expResult = 1.85+0.224500;
        result = instance.getVDWRmin(nonbond);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getVDWEner method, of class CharmmForceField.
     */
    @Test
    public void testGetVDWEner() {
        System.out.println("getVDWEner");

        MolecularNonbonded nonbond = new MolecularNonbonded();
        nonbond.setAtomTypeA(instance.getAtomType("CP3"));
        nonbond.setAtomTypeB(instance.getAtomType("CPB"));
        double expResult = Math.sqrt(-0.055*-0.09);
        double result = instance.getVDWEner(nonbond);
        assertTrue( Math.abs(result-expResult)<1e-12 );

        nonbond = new MolecularNonbonded();
        nonbond.setAtomTypeA(instance.getAtomType("NY"));
        nonbond.setAtomTypeB(instance.getAtomType("HA3"));
        expResult = Math.sqrt(-0.20*-0.024);
        result = instance.getVDWEner(nonbond);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    
    /**
     * Test of getChargeAB method, of class CharmmForceField.
     */
    @Test
    public void testGetChargeAB() {
        System.out.println("getChargeAB");
        MolecularPair pair = new MolecularPair();
        pair.setChargeA(0.5);
        pair.setChargeB(0.007);
        double expResult = 0.5*0.007;
        double result = instance.getChargeAB(pair);
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getKProtLigDielectric method, of class CharmmForceField.
     */
    @Test
    public void testGetKProtLigDielectric() {
        System.out.println("getKProtLigDielectric");
        double expResult = 1.0;
        double result = instance.getKProtLigDielectric();
       assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getKWaterDielectric method, of class CharmmForceField.
     */
    @Test
    public void testGetKWaterDielectric() {
        System.out.println("getKWaterDielectric");
        double expResult = 78.4;
        double result = instance.getKWaterDielectric();
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

    /**
     * Test of getKSurfTensionWater method, of class CharmmForceField.
     */
    @Test
    public void testGetKSurfTensionWater() {
        System.out.println("getKSurfTensionWater");
        double expResult = 0.0072;
        double result = instance.getKSurfTensionWater();
        assertTrue( Math.abs(result-expResult)<1e-12 );
    }

   
  
}