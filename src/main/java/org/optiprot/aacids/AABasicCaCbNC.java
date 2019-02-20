/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.aacids;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.io.PDBParseException;
import org.optiprot.OptiProtParameters;
import org.optiprot.rotamer.RotamerLibrary;
import org.biojava.bio.structure.*;
import org.biojava.bio.structure.jama.Matrix;


/**
 * Implementas a reduced representation of an aminoacid, with Calpha, Cbeta and
 * the sidechain's center of mass
 *
 * @author victor
 */
public class AABasicCaCbNC implements java.io.Serializable, IAABasic {

    private OptiProtParameters parameters=null;

    static private final int NUM_ATOMS = 5;

    private Atom [] atoms = new Atom[this.getNumAtoms()];

    static private final int CA = 0;
    static private final int CB = 1;
    static private final int CoM = 2;
    static private final int N = 3;
    static private final int C = 4;

    private String name="";

    private Integer rotIdx=null;

    


    public AABasicCaCbNC(){

    }

    public AABasicCaCbNC(Group aminoacid, OptiProtParameters par )
            throws StructureException {

        if( aminoacid instanceof AminoAcid ){
            
            AminoAcid aa=(AminoAcid)aminoacid;

            this.setParameters(par);

            this.setCa( (Atom)aa.getCA().clone() );

            try{
                this.setCb( (Atom)aa.getCB().clone() );
            }
            catch(StructureException ex){
                //case of glycine
                this.setCb( Calc.createVirtualCBAtom(aa) );
            }

            this.setN( (Atom)aa.getN().clone() );
            this.setC( (Atom)aa.getC().clone() );

            this.setName( aa.getPDBName() );

            this.setRotIdx( par.getRotLib().getRotamerIdx(aa) );

            this.setCoM(
                    Calc.getCentroid(  RotamerLibrary.getSideChainAtoms(aa) ) );
            this.getCoM().setName( "CoM" );


        }
        else{
            throw new ClassCastException();
        }

    }

    

    /**
     * @return the cb
     */
    public Atom getCb() {
        return this.getAtoms()[AABasicCaCbNC.CB];
    }

    /**
     * @param cb the cb to set
     */
    public void setCb(Atom cb) {
        this.getAtoms()[AABasicCaCbNC.CB] = cb;
    }

    /**
     * @return the coM
     */
    public Atom getCoM() {
        return this.getAtoms()[AABasicCaCbNC.CoM];
    }

    /**
     * @param coM the coM to set
     */
    public void setCoM(Atom coM) {
        this.getAtoms()[AABasicCaCbNC.CoM] = coM;
    }

    /**
     * @return the rotIdx
     */
    public Integer getRotIdx() {
        return rotIdx;
    }

    /**
     * @param rotIdx the rotIdx to set
     */
    private void setRotIdx(Integer rotIdx) {
        this.rotIdx = rotIdx;
    }

    

    /**
     * changes the rotamer index randomly
     *
     * @param i
     * @param v
     */
    public void applyMutationChangeRotamer() {

        int nrotamer=this.getParameters().getRotLib().getNumAngles(this.getName());

        if(nrotamer<=1)
            return;

        int newIdx=(int)Math.floor( Math.random()*(double)nrotamer );

        while( this.getRotIdx()==newIdx ){
            newIdx=(int)Math.floor( Math.random()*(double)nrotamer );
        }

        this.changeRotIdx(newIdx);
    }

    public void applyMutation(int i, double v) {

        applyMutationChangeRotamer();
    }

    
    public Atom[] getAtoms() {

        return atoms;
    }

    public int getNumAtoms() {
        return AABasicCaCbNC.NUM_ATOMS;
    }

    /**
     * @return the N
     */
    public Atom getN() {
        return this.getAtoms()[AABasicCaCbNC.N];
    }

    /**
     * @param N the N to set
     */
    public void setN(Atom N) {
        this.getAtoms()[AABasicCaCbNC.N] = N;
    }


    /**
     * @return the C
     */
    public Atom getC() {
        return this.getAtoms()[AABasicCaCbNC.C];
    }

    /**
     * @param C the C to set
     */
    public void setC(Atom C) {
        this.getAtoms()[AABasicCaCbNC.C] = C;
    }

    public Atom getVectorNCa() {
        try {
            return Calc.substract(this.getCa(), this.getN());
        } catch (StructureException ex) {
            return null;
        }
    }

    public Atom getVectorCaC() {
        try {
            return Calc.substract(this.getC(), this.getCa());
        } catch (StructureException ex) {
            return null;
        }
    }

    /**
     * @return the parameters
     */
    public OptiProtParameters getParameters() {
        return parameters;
    }

    /**
     * @param parameters the parameters to set
     */
    public void setParameters(OptiProtParameters parameters) {
        this.parameters = parameters;
    }


    /**
     * @return the ca
     */
    public Atom getCa() {
        return this.getAtoms()[AABasicCaCbNC.CA];
    }

    /**
     * @param ca the ca to set
     */
    public void setCa(Atom ca) {
        this.getAtoms()[AABasicCaCbNC.CA] = ca;
    }


    /**
     * @param atoms the atoms to set
     */
    private void setAtoms(Atom[] atoms) {
        this.atoms = atoms;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    @Override
    public Object clone(){

        AABasicCaCbNC obj=new AABasicCaCbNC();

        obj.setC( (Atom)this.getC().clone() );
        obj.setN( (Atom)this.getN().clone() );
        obj.setCa( (Atom)this.getCa().clone() );
        obj.setCoM( (Atom)this.getCoM().clone() );
        obj.setCb( (Atom)this.getCb().clone() );

        obj.setName( this.getName() );
        obj.setRotIdx( this.getRotIdx() );
        obj.setParameters( this.getParameters() );

        return obj;
    }

    public AminoAcid toAminoAcid(boolean generatesH) {

        AminoAcid aa = new AminoAcidImpl();

        try {
            aa.setPDBName( this.getName() );
        } catch (PDBParseException ex) {
            Logger.getLogger(AABasicCaCbNC.class.getName()).log(Level.SEVERE, null, ex);
        }

        aa.addAtom((Atom) this.getN().clone());
        aa.addAtom((Atom) this.getCa().clone());
        aa.addAtom((Atom) this.getC().clone());

        // nots glycine -> add  Cb
        if( !aa.getPDBName().toUpperCase().equals("GLY") ){
            aa.addAtom((Atom) this.getCb().clone());
        }

        Atom[] arr1 = new Atom[3];
        Atom[] arr2 = new Atom[3];

        arr2[0] = this.getN();
        arr2[1] = this.getCa();
        arr2[2] = this.getC();

        AminoAcid rotamer =
                (AminoAcid) getParameters().getRotLib().getRotamer( getName(), getRotIdx());

        //calculate HA & sidechain atoms
        try {

            arr1[0] = rotamer.getN();
            arr1[1] = rotamer.getCA();
            arr1[2] = rotamer.getC();

            List<Atom> atomsToTransform=new ArrayList<Atom>();

            // get O
            atomsToTransform.add( (Atom)rotamer.getO().clone() );

            if( !aa.getPDBName().toUpperCase().equals("GLY") ){
                //nots Glycine

                //add HA
                if( generatesH ){
                    atomsToTransform.add( (Atom)rotamer.getAtom("HA").clone() );
                }

                //add sidechain
                Atom [] sidechain = this.getParameters().getRotLib().getSideChainAtoms(
                        aa.getPDBName(),
                        this.getRotIdx() );

                for(int i=0;i<sidechain.length;i++){

                    if( !generatesH && RotamerLibrary.isHidrogen(sidechain[i]) ) {
                            continue;
                    }

                    atomsToTransform.add( (Atom)sidechain[i].clone() );
                }
            }

            superImposer(atomsToTransform,arr1,arr2);

            aa.getAtoms().addAll(atomsToTransform);


        } catch (StructureException ex) {
            Logger.getLogger(AABasicCaCbNC.class.getName()).log(Level.SEVERE, null, ex);
        }

        return aa;
    }

    /**
     * Transforms the atoms in src given in coordinates of arr1 to coordinates
     * of arr2
     *
     * @param src
     * @param arr1 src base
     * @param arr2 target base
     * @throws org.biojava.bio.structure.StructureException
     */
    private void superImposer( List<Atom> src, Atom[] arr1, Atom[] arr2)
            throws StructureException{

        // ok now we got the two arrays, do a SVD:
        SVDSuperimposer svd = new SVDSuperimposer(arr2, arr1);

        Matrix rotMatrix = svd.getRotation();
        Atom tranMatrix = svd.getTranslation();

        Iterator<Atom> it=src.iterator();

        while(it.hasNext()){
            Atom atomTmp=it.next();

            Calc.rotate(atomTmp,rotMatrix);
            Atom atomTmp2= Calc.add(atomTmp,tranMatrix);
            atomTmp.setCoords( atomTmp2.getCoords() );
        }

    }

    /**
     * gets the atom by name
     * 
     * @param name
     * @return
     */
    public Atom getAtom(String name) {

        if( atoms[0].getName().equals(name) )
            return atoms[0];

        if( atoms[1].getName().equals(name) )
            return atoms[1];

        if( atoms[2].getName().equals(name) )
            return atoms[2];

        if( atoms[3].getName().equals(name) )
            return atoms[3];

        if( atoms[4].getName().equals(name) )
            return atoms[4];

        return null;
    }

    /**
     * recalculates the coordinates of CoM, used after a rotamer index change
     */
    private void recalcCoM(){

        int idx=this.getRotIdx();

        Atom com =
                (Atom) getParameters().getRotLib().getSideChainCoM(getName(), idx).clone();

        Atom[] arr1 = new Atom[3];
        Atom[] arr2 = new Atom[3];

        arr2[0] = this.getN();
        arr2[1] = this.getCa();
        arr2[2] = this.getC();

        AminoAcid rotamer =
                (AminoAcid) getParameters().getRotLib().getRotamer( getName(), idx);
        try {
            arr1[0] = rotamer.getN();
            arr1[1] = rotamer.getCA();
            arr1[2] = rotamer.getC();

            // ok now we got the two arrays, do a SVD:
            SVDSuperimposer svd = new SVDSuperimposer(arr2, arr1);

            Calc.rotate( com, svd.getRotation());
            com = Calc.add( com, svd.getTranslation());


        } catch (StructureException ex) {
            Logger.getLogger(AABasicCaCbNC.class.getName()).log(Level.SEVERE, null, ex);
        }
        finally{
            this.getCoM().setX( com.getX() );
            this.getCoM().setY( com.getY() );
            this.getCoM().setZ( com.getZ() );
        }
        
    }

    /**
     * changes the rotamer index and recalculates the CoM
     * @param rotIdx
     */
    public void changeRotIdx(Integer rotIdx) {
        this.setRotIdx(rotIdx);
        recalcCoM();
    }

}
