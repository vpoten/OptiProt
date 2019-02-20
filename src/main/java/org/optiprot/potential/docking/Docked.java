/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking;

import org.optiprot.metaheuristic.de.IDESolutionEval;
import org.optiprot.potential.docking.element.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.jama.Matrix;
import org.optiprot.io.mol2.Mol2Atom;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcTransform;
import org.optiprot.potential.docking.element.DockingAtomClassify.simpletype;

/**
 *
 * @author victor
 */
public class Docked {


    private DockingProtein protein=null;
    private DockingLigand ligand = null;

    private Matrix matTransf=Matrix.identity(4, 4);
    private Matrix matAtomFix=Matrix.identity(4, 4);
    private Matrix matActPointFix=Matrix.identity(4, 4);

    private ArrayList<Atom> transLigAtoms=new ArrayList<Atom>();
    private ArrayList<Atom> transLigAtomsC1=new ArrayList<Atom>();
    private ArrayList<Atom> transLigAtomsC2=new ArrayList<Atom>();
    private ArrayList<Atom> origLigAtoms=new ArrayList<Atom>();
    
    private ArrayList<DockingActivePoint> ligActiveAtoms=new ArrayList<DockingActivePoint>();

    //lists of compatible active points for each ligand active atom
    private ArrayList<ArrayList<Integer>> ligCompActivePoints=new ArrayList<ArrayList<Integer>>();
    //lists of compatible active atoms for each protein active point
    private ArrayList<ArrayList<Integer>> protCompActiveAtom=new ArrayList<ArrayList<Integer>>();

    private IDESolutionEval evaluator=null;

    //optimization using h-bond favorable positions or nonpolar favorable positions
    private boolean modePolar=true;

    static private int MIN_ACT_ATOMS=3;//min active atoms in ligand to switch modePolar
    

    /**
     *
     * @param protein
     * @param ligand
     * @param center
     * @param grid used to calculate the active site points
     */
    public Docked(DockingProtein protein, DockingLigand ligand, Atom center, AtomGrid grid) {
        setProtein(protein);
        setLigand(ligand);

        getProtein().calcActivePoints( center, getLigand().calcRadii(), grid, isModePolar() );

        createLigandActiveAtoms();

        if( getLigActiveAtoms().size()<MIN_ACT_ATOMS || !isAbleAssignActivePoints() ){
            //if the ligand dont have active atoms (h-bond) use nonpolar (C)
            setModePolar(false);
            getProtein().calcActivePoints( center, getLigand().calcRadii(), grid, isModePolar() );
            createLigandActiveAtoms();
        }
        
    }

    
    
    /**
     * gets the ligand atoms (transformed by the matrix)
     *
     * @return
     */
    public List<Atom> getLigandAtoms() {
        return getTransAtoms();
    }

    public List<Atom> getLigandAtomsCopy1() {

        for( int i=0 ; i<getTransAtoms().size(); i++ ){
            Atom atC=transLigAtomsC1.get(i);
            Atom atL=getTransAtoms().get(i);

            atC.setX( atL.getX() );
            atC.setY( atL.getY() );
            atC.setZ( atL.getZ() );
        }

        return transLigAtomsC1;
    }

    public List<Atom> getLigandAtomsCopy2() {
        
        for( int i=0 ; i<getTransAtoms().size(); i++ ){
            Atom atC=transLigAtomsC2.get(i);
            Atom atL=getTransAtoms().get(i);

            atC.setX( atL.getX() );
            atC.setY( atL.getY() );
            atC.setZ( atL.getZ() );
        }

        return transLigAtomsC2;
    }

    
    /**
     * @return the protein
     */
    public DockingProtein getProtein() {
        return protein;
    }

    /**
     * gets the protein's atoms of active site
     * 
     * @return
     */
    public Collection<Atom> getProteinAtoms() {
        return getProtein().getActSiteProtAtoms();
    }

    /**
     * @param protein the protein to set
     */
    protected void setProtein(DockingProtein protein) {
        this.protein = protein;
    }

    /**
     * @return the ligand
     */
    public DockingLigand getLigand() {
        return ligand;
    }

    /**
     * @param ligand the ligand to set
     * @pre : the new ligand is the same molecule (another conformer)
     */
    public void setLigand(DockingLigand ligand) {

        boolean firstime=false;

        if( this.ligand==null )
            firstime=true;

        this.ligand = ligand;

        if( !getTransAtoms().isEmpty() )
            getTransAtoms().clear();

        if( !getOrigAtoms().isEmpty() )
            getOrigAtoms().clear();

        calcLigandTransAtoms();

        if( !firstime )
            createLigandActiveAtoms();
    }

    /**
     * sets the geometrical transformation and transforms the atoms of ligand
     *
     * @param axis
     * @param angle in radians
     * @param shift
     */
    public void setTransform(Atom axis, double angle, Atom shift){

        setTransform( getMatTransf(), axis, angle, shift);

        calcLigandTransAtoms();
    }

    public void setTransform(Matrix mat){

        getMatTransf().setMatrix(0, 3, 0, 3, mat);

        calcLigandTransAtoms();
    }

    /**
     *
     * @param mat : I/O
     * @param axis
     * @param angle in radians
     * @param shift
     */
    public static void setTransform( Matrix mat, Atom axis, double angle, Atom shift){

        double [][] m1=CalcTransform.matrixFromAxisAngle(axis, angle);

        mat.set(0, 0, m1[0][0]);
        mat.set(0, 1, m1[0][1]);
        mat.set(0, 2, m1[0][2]);
        mat.set(1, 0, m1[1][0]);
        mat.set(1, 1, m1[1][1]);
        mat.set(1, 2, m1[1][2]);
        mat.set(2, 0, m1[2][0]);
        mat.set(2, 1, m1[2][1]);
        mat.set(2, 2, m1[2][2]);

        mat.set(0, 3, shift.getX() );
        mat.set(1, 3, shift.getY() );
        mat.set(2, 3, shift.getZ() );
    }

    /**
     * @return the matTransf
     */
    protected Matrix getMatTransf() {
        return matTransf;
    }

    /**
     * this matrix translates the ligand's fixed atom to origin
     *
     * @return
     */
    protected Matrix getMatAtomFix() {
        return matAtomFix;
    }

    /**
     * this matrix translates the ligand's fixed atom on origin to
     * active point
     *
     * @return
     */
    protected Matrix getMatActPointFix() {
        return matActPointFix;
    }

    

    /**
     * creates ligand active atoms
     */
    private void createLigandActiveAtoms(){

        getLigActiveAtoms().clear();

        for( int i=0 ; i<getLigand().getAtoms().size(); i++ ){

            Atom atO=getOrigAtoms().get(i);
            Atom atL=getLigand().getAtoms().get(i);

            //checks if the atom is active
            DockingActivePoint actp=null;

            if( isModePolar() )
                actp=DockingMolecule.isActiveAtom((Mol2Atom)atL);
            else
                actp=DockingMolecule.isNonpolarAtom((Mol2Atom)atL);

            if( actp!=null ){
                //if is active store it in and keep the original coordinates
                actp.setPoint(atO);
                getLigActiveAtoms().add( actp );
            }

        }

        createLigandCompatibleList();
    }

    /**
     * transforms ligand atoms
     *
     * @return the transLigAtoms
     */
    private void calcLigandTransAtoms() {

        if( getTransAtoms().isEmpty() ){
            for( int i=0 ; i<getLigand().getAtoms().size(); i++ ){

                if( !(getLigand().getAtoms().get(i) instanceof Mol2Atom) )
                    throw new ClassCastException("Not Mol2Atom");

                Atom atO=new AtomImpl();
                Atom atL=getLigand().getAtoms().get(i);

                atO.setName( atL.getName() );
                atO.setX( atL.getX() );
                atO.setY( atL.getY() );
                atO.setZ( atL.getZ() );

                getOrigAtoms().add(atO);
                transLigAtomsC1.add( (Atom)atO.clone() );
                transLigAtomsC2.add( (Atom)atO.clone() );
                getTransAtoms().add( atL );

            }
        }

        Matrix matFinal = getFinalTransform();

        //copy original coordinates and transform by the matrix
        for( int i=0 ; i<getOrigAtoms().size(); i++ ){

            Atom atO=getOrigAtoms().get(i);
            Atom atL=getTransAtoms().get(i);

            atL.setX( atO.getX() );
            atL.setY( atO.getY() );
            atL.setZ( atO.getZ() );
            
            CalcTransform.applyTransform( atL, matFinal );
        }

    }

    /**
     * create list of compatible active points for ligand active atoms
     */
    private void createLigandCompatibleList() {

        //create list of ligand active compatible points
        getLigCompActivePoints().clear();

        for( DockingActivePoint dlig : getLigActiveAtoms() ){

            ArrayList<Integer> listIdx=new ArrayList<Integer>();
            getLigCompActivePoints().add( listIdx );

            simpletype ligtype=dlig.getLigandType();

            int index=0;

            for( DockingActivePoint dp : getProtein().getActivePoints() ){

                simpletype pointtype=dp.getLigandType();

                if( DockingActivePoint.isTypeComp(ligtype, pointtype) ){
                    listIdx.add(index);
                }

                index++;
            }
        }

        //create list of active points compatible active ligand atoms
        getProtCompActiveAtom().clear();

        for( DockingActivePoint dprot : getProtein().getActivePoints() ){

            ArrayList<Integer> listIdx=new ArrayList<Integer>();
            getProtCompActiveAtom().add( listIdx );

            simpletype pointtype=dprot.getLigandType();

            int index=0;

            for( DockingActivePoint dlig : getLigActiveAtoms() ){

                simpletype ligtype=dlig.getLigandType();

                if( DockingActivePoint.isTypeComp(ligtype, pointtype) ){
                    listIdx.add(index);
                }

                index++;
            }
        }
    }

    /**
     * gets the final transformation matrix
     * 
     * @return
     */
    private Matrix getFinalTransform() {
        return getMatActPointFix().times( getMatTransf().times( getMatAtomFix() ) );
    }

    /**
     * @return the transLigAtoms
     */
    private ArrayList<Atom> getTransAtoms() {
        return transLigAtoms;
    }

    /**
     * @return the origLigAtoms
     */
    private ArrayList<Atom> getOrigAtoms() {
        return origLigAtoms;
    }

    
    protected ArrayList<DockingActivePoint> getLigActiveAtoms() {
        return ligActiveAtoms;
    }


    

    /**
     * gets the number of active points in binding site
     * 
     * @return
     */
    public int getNumActPoints(){
        return this.getProtein().getActivePoints().size();
    }

    /**
     * gets the number of active atoms of ligand
     * 
     * @return
     */
    public int getNumActAtoms(){
        return this.getLigActiveAtoms().size();
    }

    /**
     * gets an index of an ligand active atom with active point
     * (id_point) compatible type
     *
     * @param id_point
     * @return
     */
    public Integer getCompActAtom( int id_point ){

        ArrayList<Integer> listIdx=getProtCompActiveAtom().get(id_point);

        if( listIdx.isEmpty() )
            return null;

        return listIdx.get(  (int)Math.floor(listIdx.size()*Math.random()) );

    }

    /**
     * gets an index of an active point with ligand active atom
     * (id_atom) compatible type
     *
     * @param id_atom
     * @return
     */
    public Integer getCompActPoint( int id_atom ){

        ArrayList<Integer> listIdx=getLigCompActivePoints().get(id_atom);

        if( listIdx.isEmpty() )
            return null;

        return listIdx.get( (int)Math.floor(listIdx.size()*Math.random()) );
    }


    /**
     *
     * @param id_atom : index of ligand active atom
     * @param id_point : index of binding site active point
     */
    public void assignAtomToActivePoint( int id_atom, int id_point ){

        Atom al=getLigActiveAtoms().get(id_atom).getPoint();
        //-L
        getMatAtomFix().set( 0, 3, -al.getX());
        getMatAtomFix().set( 1, 3, -al.getY());
        getMatAtomFix().set( 2, 3, -al.getZ());

        Atom ap=getProtein().getActivePoints().get(id_point).getPoint();
        //P
        getMatActPointFix().set( 0, 3, ap.getX());
        getMatActPointFix().set( 1, 3, ap.getY());
        getMatActPointFix().set( 2, 3, ap.getZ());

        ///calcLigandTransAtoms();
    }

    /**
     * calculates the overlap factor of the ligand atoms inside the pocket
     * 
     * @return
     */
    public double calcOverlapFactor(){

        return getProtein().calcOverlapFactor( this.getLigandAtoms() );
    }

    /**
     * @param evaluator the evaluator to set
     */
    public void setEvaluator(IDESolutionEval evaluator) {
        this.evaluator = evaluator;
    }

    public IDESolutionEval getEvaluator() {
        return this.evaluator;
    }

    /**
     * @return the modePolar
     */
    public boolean isModePolar() {
        return modePolar;
    }

    /**
     * @param modePolar the modePolar to set
     */
    protected void setModePolar(boolean modePolar) {
        this.modePolar = modePolar;
    }

    /**
     * 
     * @return : true if exists active points in binding site and 
     * active atoms in ligand and those points are compatible
     */
    public boolean isAbleAssignActivePoints(){

        if( getLigActiveAtoms().isEmpty() ||
                getProtein().getActivePoints().isEmpty() )
            return false;

        //checks if exists compatible active points is binding site
        for(int i=0;i<getNumActAtoms();i++){
            Integer val=getCompActPoint(i);

            if(val!=null)
                return true;
        }

        return false;
    }

    /**
     * gets, for current conformation, the ligand atom in favorable docking points
     *
     * @param list : I/O
     * @param dist : distance threshold
     */
    public void getLigandAtomsInActivePoints(List<Atom> list, double dist ){

        double sqdist=dist*dist;
        Matrix matFinal = getFinalTransform();
        Atom atransL=new AtomImpl();

        for( DockingActivePoint ligAt : this.getLigActiveAtoms() ){
            for( DockingActivePoint actP : this.getProtein().getActivePoints() ){

                if( DockingActivePoint.isTypeComp(ligAt.getLigandType(), actP.getLigandType()) ){

                    atransL.setName( ligAt.getPoint().getName() );
                    atransL.setX( ligAt.getPoint().getX() );
                    atransL.setY( ligAt.getPoint().getY() );
                    atransL.setZ( ligAt.getPoint().getZ() );
                    CalcTransform.applyTransform( atransL, matFinal );

                    if( CalcGeom.squareDistance( atransL, actP.getPoint()) < sqdist ){

                        list.add( (Atom)atransL.clone() );
                        break;
                    }
                }
            }
        }

    }

    /**
     * @return the ligCompActivePoints
     */
    private ArrayList<ArrayList<Integer>> getLigCompActivePoints() {
        return ligCompActivePoints;
    }

    /**
     * @param ligCompActivePoints the ligCompActivePoints to set
     */
    private void setLigCompActivePoints(ArrayList<ArrayList<Integer>> ligCompActivePoints) {
        this.ligCompActivePoints = ligCompActivePoints;
    }

    /**
     * @return the protCompActiveAtom
     */
    private ArrayList<ArrayList<Integer>> getProtCompActiveAtom() {
        return protCompActiveAtom;
    }

    /**
     * @param protCompActiveAtom the protCompActiveAtom to set
     */
    private void setProtCompActiveAtom(ArrayList<ArrayList<Integer>> protCompActiveAtom) {
        this.protCompActiveAtom = protCompActiveAtom;
    }
    
    
}
