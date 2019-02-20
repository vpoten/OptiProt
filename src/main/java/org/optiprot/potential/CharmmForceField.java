/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import java.util.HashMap;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.optiprot.potential.element.CharmmElement;
import org.optiprot.potential.element.MolecularAngle;
import org.optiprot.potential.element.MolecularBond;
import org.optiprot.potential.element.MolecularDihedral;
import org.optiprot.potential.element.MolecularImproper;
import org.optiprot.potential.element.MolecularNonbonded;
import org.optiprot.potential.element.MolecularPair;
import org.optiprot.rotamer.RotamerLibrary;

/**
 *
 * @author victor
 */
public class CharmmForceField implements IForceField, IForceFieldTable {

    private float loadFactor=0.3f;

    //bonds params table
    private HashMap<String,double []> m_bonds=new HashMap<String,double []>(50,loadFactor);

    //angles params table
    private HashMap<String,double []> m_angles=new HashMap<String,double []>(50,loadFactor);

    //dihedrals params table
    private HashMap<String,double []> m_dihedrals=new HashMap<String,double []>(50,loadFactor);

    //impropers params table
    private HashMap<String,double []> m_impropers=new HashMap<String,double []>(50,loadFactor);

    //nonbondend params array
    private double [][] m_nonbondeds=null;

    static final private String TABLE_SEP="_";
    static final private String WILDCARD="X";

    //elements indexed by atom type
    CharmmElement [] elements=null;
    HashMap<String,CharmmElement> elemTable=new HashMap<String,CharmmElement>(50,loadFactor);



    ////////////////////////

    public CharmmForceField(List<CharmmElement> elements) {

        this.elements=elements.toArray(new CharmmElement [elements.size()] );
        this.m_nonbondeds=new double [elements.size()][4];

        //fill the table indexed by elem type
        
        for( CharmmElement ele : this.elements ){
            if(ele!=null)
                this.elemTable.put( ele.getType(), ele);
        }
    }

    protected String getAtomType( int id ){
        return this.elements[id].getType();
    }

    public Integer getAtomType( String atomType ){

        if( this.elemTable.get(atomType)!=null )
            return this.elemTable.get(atomType).getIndex();

        return CharmmElement.WILDCARD;
    }
    
    public double getKBond(MolecularBond bond) {

        double [] params=this.getBond( bond.getAtomTypeA(),
                bond.getAtomTypeB() );

        if( params==null )
            return 0;

        return params[0];
    }

    public double getEqDistBond(MolecularBond bond) {

        double [] params=this.getBond( bond.getAtomTypeA(),
                bond.getAtomTypeB() );

        if( params==null )
            return 0;

        return params[1];
    }

    public double getKUreyBradley(MolecularAngle angle) {
        
        double [] params=this.getAngle(
                angle.getAtomType(0),
                angle.getAtomType(1),
                angle.getAtomType(2) );

        if( params==null )
            return 0;

        return params[2];
    }

    public double getEqDistUreyBradley(MolecularAngle angle) {

        double [] params=this.getAngle( 
                angle.getAtomType(0),
                angle.getAtomType(1),
                angle.getAtomType(2) );

        if( params==null )
            return 0;

        return params[3];
    }

    public double getKBondAngle(MolecularAngle angle) {

        double [] params=this.getAngle( 
                angle.getAtomType(0),
                angle.getAtomType(1),
                angle.getAtomType(2) );

        if( params==null )
            return 0;

        return params[0];
    }

    public double getEqDistBondAngle(MolecularAngle angle) {

        double [] params=this.getAngle(
                angle.getAtomType(0),
                angle.getAtomType(1),
                angle.getAtomType(2) );

        if( params==null )
            return 0;

        return params[1];
    }

    public double getKDihedral(MolecularDihedral dihedral) {

        double [] params=this.getDihedral(
                dihedral.getAtomType(0),
                dihedral.getAtomType(1),
                dihedral.getAtomType(2),
                dihedral.getAtomType(3) );

        if( params==null )
            return 0;

        return params[0];
    }

    public double getNDihedral(MolecularDihedral dihedral) {

        double [] params=this.getDihedral(
                dihedral.getAtomType(0),
                dihedral.getAtomType(1),
                dihedral.getAtomType(2),
                dihedral.getAtomType(3) );

        if( params==null )
            return 0;

        return params[1];
    }

    public double getDDihedral(MolecularDihedral dihedral) {

        double [] params=this.getDihedral(
                dihedral.getAtomType(0),
                dihedral.getAtomType(1),
                dihedral.getAtomType(2),
                dihedral.getAtomType(3) );

        if( params==null )
            return 0;

        return params[2];
    }

    public double getKTorsionAngle(MolecularImproper improper) {

        double [] params=this.getImproper(
                improper.getAtomType(0),
                improper.getAtomType(1),
                improper.getAtomType(2),
                improper.getAtomType(3) );

        if( params==null )
            return 0;

        return params[0];
    }

    public double getEqDistTorsionAngle(MolecularImproper improper) {

        double [] params=this.getImproper(
                improper.getAtomType(0),
                improper.getAtomType(1),
                improper.getAtomType(2),
                improper.getAtomType(3) );

        if( params==null )
            return 0;

        return params[1];
    }

    public double getVDWRmin(MolecularNonbonded nonbond) {

        double r1=0;
        double r2=0;

        double [] params=this.getNonbonded(
                nonbond.getAtomTypeA() );

        if(params!=null)
            r1=params[1];

        params =this.getNonbonded(
                nonbond.getAtomTypeB() );

        if(params!=null)
            r2=params[1];

        return r1+r2;
    }

    public double getVDWEner(MolecularNonbonded nonbond) {

        double eps1=0;
        double eps2=0;

        double [] params=this.getNonbonded(
                nonbond.getAtomTypeA() );

        if(params!=null)
            eps1=params[0];

        params =this.getNonbonded(
                nonbond.getAtomTypeB() );

        if(params!=null)
            eps2=params[0];

        return Math.sqrt(eps1*eps2);
    }

    
    public double getVDWRadius(Atom a) {

        char symbol='\0';

        if( RotamerLibrary.isHidrogen(a) ){
            symbol='H';
        }
        else{
            symbol=a.getName().charAt(0);
        }

        switch(symbol) {
           case 'C':
                return 1.70;
           case 'N':
                return 1.55;
           case 'O':
                return 1.52;
           case 'H':
                return 1.09;
           case 'S':
                return 1.80;
           default:
                return 2.00;
        }

    }

    public double getInvVDWRadius(Atom a) {

        char symbol='\0';

        if( RotamerLibrary.isHidrogen(a) ){
            symbol='H';
        }
        else{
            symbol=a.getName().charAt(0);
        }

        switch(symbol) {
           case 'C':
                return 0.58823;
           case 'N':
                return 0.64516;
           case 'O':
                return 0.65789;
           case 'H':
                return 0.91743;
           case 'S':
                return 0.55555;
           default:
                return 0.5;
        }

    }

    public double getChargeA(MolecularNonbonded nonbond) {
        return nonbond.getChargeA();
    }

    public double getChargeB(MolecularNonbonded nonbond) {
        return nonbond.getChargeB();
    }

//    public double getChargeAB(MolecularNonbonded nonbond) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }

    public double getChargeAB(MolecularPair pair) {
        return pair.getChargeA()*pair.getChargeB();
    }


    public double getKProtLigDielectric() {
        return 1.0;
    }


    public double getKWaterDielectric() {
        return 78.4;
    }


    public double getKSurfTensionWater() {
        return 0.0072;
    }


    public void addBond(Integer idat1, Integer idat2, double kb, double b0) {

        String at1=getAtomType(idat1);
        String at2=getAtomType(idat2);

        this.getBonds().put(at1+TABLE_SEP+at2, new double [] {kb,b0} );
    }


    public void addAngle(Integer idat1, Integer idat2, Integer idat3, double ktheta, double theta0, double kub, double s0) {

        String at1=getAtomType(idat1);
        String at2=getAtomType(idat2);
        String at3=getAtomType(idat3);

        this.getAngles().put(at1+TABLE_SEP+at2+TABLE_SEP+at3,
                new double [] {ktheta,theta0,kub,s0} );
    }


    /**
     *
     * @param idat1 : id or CharmmElement.WILDCARD
     * @param idat2
     * @param idat3
     * @param idat4 : id or CharmmElement.WILDCARD
     * @param kchi
     * @param n
     * @param delta
     */
    public void addDihedral(Integer idat1, Integer idat2, Integer idat3, Integer idat4, double kchi, double n, double delta) {

        String at1= (idat1==CharmmElement.WILDCARD) ? WILDCARD : getAtomType(idat1);
        String at2=getAtomType(idat2);
        String at3=getAtomType(idat3);
        String at4= (idat4==CharmmElement.WILDCARD) ? WILDCARD : getAtomType(idat4);

        this.getDihedrals().put(at1+TABLE_SEP+at2+TABLE_SEP+at3+TABLE_SEP+at4,
                new double [] {kchi,n,delta} );
    }


    /**
     *
     * @param idat1
     * @param idat2 : id or CharmmElement.WILDCARD
     * @param idat3 : id or CharmmElement.WILDCARD
     * @param idat4
     * @param kpsi
     * @param psi0
     */
    public void addImproper(Integer idat1, Integer idat2, Integer idat3, Integer idat4, double kpsi, double psi0) {

        String at1=getAtomType(idat1);
        String at2= (idat2==CharmmElement.WILDCARD) ? WILDCARD : getAtomType(idat2);
        String at3= (idat3==CharmmElement.WILDCARD) ? WILDCARD : getAtomType(idat3);
        String at4=getAtomType(idat4);

        this.getImpropers().put(at1+TABLE_SEP+at2+TABLE_SEP+at3+TABLE_SEP+at4,
                new double [] {kpsi,psi0} );
    }


    public void addNonbonded(Integer idat1, double eps, double rmin2, double eps14, double rmin2_14) {

        this.getNonbondeds()[idat1][0]=eps;
        this.getNonbondeds()[idat1][1]=rmin2;
        this.getNonbondeds()[idat1][2]=eps14;
        this.getNonbondeds()[idat1][3]=rmin2_14;
    }


    public double[] getBond(Integer idat1, Integer idat2) {

        String at1=getAtomType(idat1);
        String at2=getAtomType(idat2);

        double [] val=this.getBonds().get(at1+TABLE_SEP+at2);

        if( val==null )
            return this.getBonds().get(at2+TABLE_SEP+at1);

        return val;
    }


    public double[] getAngle(Integer idat1, Integer idat2, Integer idat3) {

        String at1=getAtomType(idat1);
        String at2=getAtomType(idat2);
        String at3=getAtomType(idat3);

        double [] val=this.getAngles().get(at1+TABLE_SEP+at2+TABLE_SEP+at3);

        if( val==null )
            return this.getAngles().get(at3+TABLE_SEP+at2+TABLE_SEP+at1);

        return val;
    }


    public double[] getDihedral(Integer idat1, Integer idat2, Integer idat3, Integer idat4) {

        String at1= (idat1==CharmmElement.WILDCARD) ? CharmmForceField.WILDCARD : getAtomType(idat1);
        String at2=getAtomType(idat2);
        String at3=getAtomType(idat3);
        String at4= (idat4==CharmmElement.WILDCARD) ? CharmmForceField.WILDCARD : getAtomType(idat4);

        double [] val=this.getDihedrals().get(at1+TABLE_SEP+at2+TABLE_SEP+at3+TABLE_SEP+at4);

        if( val!=null )
            return val;

        val=this.getDihedrals().get(at4+TABLE_SEP+at3+TABLE_SEP+at2+TABLE_SEP+at1);

        if( val!=null )
            return val;

        val=this.getDihedrals().get(WILDCARD+TABLE_SEP+at2+TABLE_SEP+at3+TABLE_SEP+WILDCARD);

        if( val!=null )
            return val;

        val=this.getDihedrals().get(WILDCARD+TABLE_SEP+at3+TABLE_SEP+at2+TABLE_SEP+WILDCARD);

        return val;
    }


    public double[] getImproper(Integer idat1, Integer idat2, Integer idat3, Integer idat4) {

        String at1=getAtomType(idat1);
        String at2= (idat2==CharmmElement.WILDCARD) ? CharmmForceField.WILDCARD : getAtomType(idat2);
        String at3= (idat3==CharmmElement.WILDCARD) ? CharmmForceField.WILDCARD : getAtomType(idat3);
        String at4=getAtomType(idat4);

        double [] val=this.getImpropers().get(at1+TABLE_SEP+at2+TABLE_SEP+at3+TABLE_SEP+at4);

        if( val!=null )
            return val;

        val=this.getImpropers().get(at1+TABLE_SEP+WILDCARD+TABLE_SEP+WILDCARD+TABLE_SEP+at4);

        return val;
    }


    public double[] getNonbonded( Integer idat1) {

        return this.getNonbondeds()[idat1];
    }

    /**
     * @return the m_nonbondeds
     */
    private double[][] getNonbondeds() {
        return m_nonbondeds;
    }

    /**
     * @return the m_bonds
     */
    private HashMap<String, double[]> getBonds() {
        return m_bonds;
    }

    /**
     * @param m_bonds the m_bonds to set
     */
    private void setBonds(HashMap<String, double[]> m_bonds) {
        this.m_bonds = m_bonds;
    }

    /**
     * @return the m_angles
     */
    private HashMap<String, double[]> getAngles() {
        return m_angles;
    }

    /**
     * @param m_angles the m_angles to set
     */
    private void setAngles(HashMap<String, double[]> m_angles) {
        this.m_angles = m_angles;
    }

    /**
     * @return the m_dihedrals
     */
    private HashMap<String, double[]> getDihedrals() {
        return m_dihedrals;
    }

    /**
     * @param m_dihedrals the m_dihedrals to set
     */
    private void setDihedrals(HashMap<String, double[]> m_dihedrals) {
        this.m_dihedrals = m_dihedrals;
    }

    /**
     * @return the m_impropers
     */
    private HashMap<String, double[]> getImpropers() {
        return m_impropers;
    }

    /**
     * @param m_impropers the m_impropers to set
     */
    private void setImpropers(HashMap<String, double[]> m_impropers) {
        this.m_impropers = m_impropers;
    }

    public double[] getParBond(MolecularBond bond) {
        return this.getBond( bond.getAtomTypeA(), bond.getAtomTypeB() );
    }

    public double[] getParAngle(MolecularAngle angle) {
        return this.getAngle(
                angle.getAtomType(0),
                angle.getAtomType(1),
                angle.getAtomType(2) );
    }

    public double[] getParDihedral(MolecularDihedral dihedral) {
        return this.getDihedral(
                dihedral.getAtomType(0),
                dihedral.getAtomType(1),
                dihedral.getAtomType(2),
                dihedral.getAtomType(3) );
    }

    public double[] getParImproper(MolecularImproper improper) {
        return this.getImproper(
                improper.getAtomType(0),
                improper.getAtomType(1),
                improper.getAtomType(2),
                improper.getAtomType(3) );
    }

    public double[] getParVDW(MolecularNonbonded nonbond) {

        double r1=0;
        double r2=0;
        double eps1=0;
        double eps2=0;

        double [] params=this.getNonbonded(nonbond.getAtomTypeA());

        if(params!=null){
            r1=params[1];
            eps1=params[0];
        }

        params =this.getNonbonded(nonbond.getAtomTypeB());

        if(params!=null){
            r2=params[1];
            eps2=params[0];
        }

        return new double [] { Math.sqrt(eps1*eps2), r1+r2 };
    }

    
}
