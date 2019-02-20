/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.jama.Matrix;
import org.optiprot.io.mol2.Mol2Atom;
import org.optiprot.io.mol2.Mol2Bond;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcTransform;
import org.optiprot.pockets.CalcPockets;
import org.optiprot.pockets.GridProbeCluster;
import org.optiprot.potential.docking.element.DockingAtomClassify.simpletype;

/**
 *
 * @author victor
 */
public class DockingProtein  extends DockingMolecule {

    private HashSet<Atom> actSiteProtAtoms=new HashSet<Atom>();

    private Atom actSiteCenter=null;
    private BSPTree actSiteBtree=null;
    private List<Atom> actSite=null;
    private double actSiteRadii=0;

    private ArrayList<DockingActivePoint> activePoints=new ArrayList<DockingActivePoint>();
    
    static private double BOX_EDGE_LEN = 18.0;//in Angstroms

    static private double DIST_HBOND = 1.8;
    static private double DIST_METAL = 2.0;
    static private double DIST_NONPOLAR = 4.0;
    static private double DIST_POINT_PROBE = 1.5;
    static private double DIST_MAXRAD = 2.2;

    //vectors used for the calculation of active points
    private static Atom[] multActVectors=null;
    private static Atom[] multMetActVectors=null;

    static{
        multActVectors=calcMultActVectors();
        multMetActVectors=calcMultMetActVectors();
    }


    
    public DockingProtein(Chain chain) {
        super(chain);
    }

    /**
     *
     * @param center : center of active site
     * @param radii : radii of ligand
     * @param grid : auxiliar grid
     * @param polar : if calculates polar points or not
     */
    public void calcActivePoints(Atom center, double radii, AtomGrid grid, boolean polar ) {

        setActSiteCenter(center);
        setActSiteRadii(radii);
        setActSiteBtree( null );

        //calc the protein atoms in binding site
        calcActiveProtAtoms( BOX_EDGE_LEN, grid, false );
        
        try {
            //calc the active points induced by active atoms of the protein in binding site
            calcActivePoints(polar);

            if( getActivePoints().isEmpty() ){
                //if cannot calculates active points using pocket then recalculates
                // using the active site radii
                setActSiteBtree( null );
                calcActivePoints(polar);
            }
        } catch (StructureException ex) {
            Logger.getLogger(DockingProtein.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * calculates the overlap factor of the atoms inside the pocket
     *
     * @param listAtoms
     * @return
     */
    public double calcOverlapFactor(List<Atom> listAtoms) {

        int count=0;

        double sqlim=Math.pow( getActSiteRadii(),2 );

        for( Atom at : listAtoms){
            if( getActSiteBtree()!=null ){
                if( getActSiteBtree().areNeighbours(at, DIST_POINT_PROBE) )
                    count++;
            }
            else{
                //if not exists the bsptree uses active site radii
                if( CalcGeom.squareDistance(at, getActSiteCenter())<sqlim )
                    count++;
            }
        }

        return count/(double)listAtoms.size();
    }

    
    /**
     * gets the amino acids near to an atom
     * 
     * @param at
     * @param dist
     * @return
     */
    public Set<AminoAcid> getNeighResidues( Atom at, double dist ){

        Set<AminoAcid> set=new HashSet<AminoAcid>();

        getBtree().neighbours(at, dist, getListAtomAux() );

        for( Atom atom : getListAtomAux() ){
            if( atom.getParent() instanceof AminoAcid )
                set.add( (AminoAcid)atom.getParent() );
        }

        if( !getListAtomAux().isEmpty() )
            getListAtomAux().clear();

        return set;
    }


    /**
     * appends to list the atoms near to at
     *
     * @param at
     * @param dist max distance of neighbors
     * @param list : I/O
     */
    public void getNeighAtoms( Atom at, double dist, Collection<Atom> list ){


        getBtree().neighbours(at, dist, getListAtomAux() );

        for( Atom atom : getListAtomAux() ){
            list.add(atom);
        }

        if( !getListAtomAux().isEmpty() )
            getListAtomAux().clear();

    }

    

    /**
     * calculates the probes that aproximates the active site
     *
     * @param actSiteCenter
     * @param side_length
     * @param grid
     * @return
     */
    private GridProbeCluster calcActiveSite( Atom actSiteCenter, double side_length, AtomGrid grid ){
        try {
            return CalcPockets.calcActiveSite(getBtree(), actSiteCenter, side_length, grid);
        } catch (StructureException ex) {
            Logger.getLogger(DockingProtein.class.getName()).log(Level.SEVERE, null, ex);
        }

        return null;
    }

    /**
     * 
     * @param actSiteCenter
     * @param radii
     * @return
     */
    private GridProbeCluster generateActiveSite(Atom actSiteCenter, double radii, AtomGrid grid ) {
        try {
            return CalcPockets.generateActiveSite(getBtree(), actSiteCenter, radii, grid);
        } catch (StructureException ex) {
            Logger.getLogger(DockingProtein.class.getName()).log(Level.SEVERE, null, ex);
        }

        return null;
    }


    public Collection<Atom> getActSiteProtAtoms() {
        return actSiteProtAtoms;
    }

    
    /**
     * @return the actSiteCenter
     */
    public Atom getActSiteCenter() {
        return actSiteCenter;
    }

    /**
     * @param actSiteCenter the actSiteCenter to set
     */
    protected void setActSiteCenter(Atom activeCenter) {
        this.actSiteCenter = activeCenter;
    }
    
    /**
     * @return the activePoints
     */
    public ArrayList<DockingActivePoint> getActivePoints() {
        return activePoints;
    }


   /**
    * calculates the atoms of protein in active site using a box centered in
    * actSiteCenter, also calculates a btree with the probes of active site
    *
    * @param edge_len : length of the box that contains the active site
    * @param grid : auxiliary grid used for pocket calculation
    * @param calcCluster : if true calculates the active site, else use a sphere
    *   with ligand radii
    */
     private void calcActiveProtAtoms( double edge_len, AtomGrid grid, boolean calcCluster) {

         if( !getActSiteProtAtoms().isEmpty() )
            getActSiteProtAtoms().clear();

         Atom activeCenter=this.getActSiteCenter();
         double radii=this.getActSiteRadii();

         getAtomsInsideBox( activeCenter, edge_len, getActSiteProtAtoms() );

         GridProbeCluster cluster=null;

         if( calcCluster )
             cluster=calcActiveSite(getActSiteCenter(), edge_len, grid);
         else
             cluster=generateActiveSite( getActSiteCenter(), radii, grid);


         if( cluster!=null ){
             //constructs btree of active site points

             ArrayList<Atom> list=new ArrayList<Atom>();
             cluster.getAtoms(list);

             if( radii>0.0 ){
                 //keep the pocket points inside the sphere centered at activeCenter

                 double sqlim=radii*radii;
                 ArrayList<Atom> list2=new ArrayList<Atom>();

                 for( Atom at : list ){
                     if( CalcGeom.squareDistance(activeCenter, at)<=sqlim )
                         list2.add(at);
                 }

                 setActSite(list2);
             }
             else{
                 setActSite(list);
             }


             this.setActSiteBtree( new BSPTree(getActSite()) );
         }

    }


     /**
     * gets the active atoms of active site's protein atoms
     *
     * @param list
     * @param polar
     */
    private void calcActiveAtomsInActSite(Collection<Atom> list, boolean polar ) {

        for( Atom at : getActSiteProtAtoms() ){

            Mol2Atom at2=(Mol2Atom) at;
            simpletype type=DockingAtomClassify.getSimpletype(at2);

            if( polar ){
                if( type==simpletype.metal || type==simpletype.acceptor ||
                        type==simpletype.donor_accept ){
                    list.add(at2);
                }
                else if( type==simpletype.donorh ){
                    list.add(at2);
                }
            }
            else{
                if( type==simpletype.nonpolar ){
                    list.add(at2);
                }
            }

        }
    }



    /**
     * calculates the active points induced by active atoms of the protein
     * in binding site
     *
     * @param polar : if true calculates points for h-bond else calculates nonpolar
     * @throws org.biojava.bio.structure.StructureException
     */
    private void calcActivePoints( boolean polar ) throws StructureException{

        ArrayList<Atom> listActAtoms=new ArrayList<Atom>();
        calcActiveAtomsInActSite(listActAtoms, polar);
        getActivePoints().clear();

        ArrayList<Atom> listaux=new ArrayList<Atom>();

        for( Atom at : listActAtoms ){
            Mol2Atom at2=(Mol2Atom) at;
            simpletype type=DockingAtomClassify.getSimpletype(at2);

            if( !listaux.isEmpty() )
                listaux.clear();

            if( polar ){
                //polara cases, h-bond
                if( type==simpletype.donorh ){
                    //a H-bond H
                    if( at2.getValence()>1 )
                        continue;// H with more than 2 bonds??

                    Atom p=calcSimpleActivePoint( at2, DIST_HBOND );

                    if( isValidActivePoint( at2, p, DIST_HBOND+0.25) )
                        getActivePoints().add( new DockingActivePoint( p, simpletype.acceptor));

                }
                else if( type==simpletype.donor_accept || type==simpletype.acceptor ){
                    //a H-bond acceptor
                    if( at2.getValence()==1 ){
                        Atom p=calcSimpleActivePoint( at2, DIST_HBOND );

                        if( isValidActivePoint( at2, p, DIST_HBOND+0.25) )
                            getActivePoints().add( new DockingActivePoint( p, simpletype.donorh));
                    }
                    else{
                        calcMultipleActPoints( at2, DIST_HBOND, listaux );

                        for(Atom p: listaux )
                            getActivePoints().add( new DockingActivePoint( p, simpletype.donorh));
                    }
                }
                else{
                    //metal ion
                    calcMultipleMetalActPoints( at2, DIST_METAL, listaux );

                    for(Atom p: listaux )
                        getActivePoints().add( new DockingActivePoint( p, simpletype.acceptor));

                }
            }
            else{
                //nonpolar case
                if( type==simpletype.nonpolar ){

                    calcMultipleMetalActPoints( at2, DIST_NONPOLAR, listaux );

                    for(Atom p: listaux )
                        getActivePoints().add( new DockingActivePoint( p, simpletype.nonpolar));

                }
            }

        }
    }

    /**
     * calculates the place for a hydrogen bond atom connected to the given
     *
     * @param atom
     * @param opt_dist
     * @return : the active point (nots checked)
     * @throws org.biojava.bio.structure.StructureException
     */
    private Atom calcSimpleActivePoint( Mol2Atom atom, double opt_dist )
            throws StructureException {

        Atom dir=Calc.substract(atom, atom.getBonds().get(0).getNeighbor(atom));

        //normalize dir and set the opt_dist length
        CalcGeom.product2(dir, opt_dist/Calc.amount(dir));

        //translate the point relative to atom position
        Atom p=Calc.add(atom, dir);

        return p;
    }

    /**
     *
     * @param atom
     * @param opt_dist
     * @param list : valid active points
     *
     * @throws org.biojava.bio.structure.StructureException
     */
    private void calcMultipleActPoints( Mol2Atom atom, double opt_dist, List<Atom> list )
            throws StructureException{

        Atom vaux=new AtomImpl();
        vaux.setX(0); vaux.setY(0); vaux.setZ(0);

        for( Mol2Bond bond : atom.getBonds() ){
            CalcGeom.addEquals(vaux, Calc.substract( bond.getNeighbor(atom), atom) );
        }

        CalcGeom.product2(vaux, -1.0/atom.getValence());//calc the median vector
        CalcGeom.product2(vaux, 1.0/Calc.amount(vaux));//normalize vaux

        Atom axis=Calc.vectorProduct(multActVectors[0], vaux);
        double angle=Math.acos( Calc.skalarProduct(multActVectors[0], vaux) );

        double [][]m1=CalcTransform.matrixFromAxisAngle(axis, angle);
        Matrix mat=Matrix.identity(4, 4);
        mat.set(0, 0, m1[0][0]);
        mat.set(0, 1, m1[0][1]);
        mat.set(0, 2, m1[0][2]);
        mat.set(1, 0, m1[1][0]);
        mat.set(1, 1, m1[1][1]);
        mat.set(1, 2, m1[1][2]);
        mat.set(2, 0, m1[2][0]);
        mat.set(2, 1, m1[2][1]);
        mat.set(2, 2, m1[2][2]);

        for(int i=0;i<multActVectors.length;i++){
            vaux=CalcGeom.product( multActVectors[i], 1);
            //transforms the vector to meet the atom-bonds relative direction
            CalcTransform.applyTransform(vaux, mat);
            //set the opt_dist length
            CalcGeom.product2(vaux, opt_dist);
            //translate the point relative to atom position
            CalcGeom.addEquals(vaux, atom);

            if( isValidActivePoint( atom, vaux, opt_dist+0.25) )
                list.add(vaux);
        }
    }

    /**
     *
     * @param atom
     * @param opt_dist
     * @param list
     * @throws org.biojava.bio.structure.StructureException
     */
    private void calcMultipleMetalActPoints( Mol2Atom atom, double opt_dist, List<Atom> list )
            throws StructureException{

        Atom vaux=null;

        for(int i=0;i<multMetActVectors.length;i++){
            vaux=CalcGeom.product( multMetActVectors[i], 1);

            //set the opt_dist length
            CalcGeom.product2(vaux, opt_dist);
            //translate the point relative to atom position
            CalcGeom.addEquals(vaux, atom);

            if( isValidActivePoint( atom, vaux, opt_dist+0.25) )
                list.add(vaux);
        }

    }

    /**
     * precalculates the places for hydrogen bond H donors
     *
     * @return
     */
    private static Atom[] calcMultActVectors() {

        Atom [] vect=new Atom [9];
        Atom vaux=new AtomImpl();

        vect[0]=new AtomImpl();
        vect[0].setX(0); vect[0].setY(1); vect[0].setZ(0);


        vaux.setX(1); vaux.setY(0); vaux.setZ(0);
        vect[1]=Calc.unitVector( Calc.add(vaux, vect[0]) );

        vaux.setX(0); vaux.setY(0); vaux.setZ(1);
        vect[2]=Calc.unitVector( Calc.add(vaux, vect[0]) );

        vaux.setX(-1); vaux.setY(0); vaux.setZ(0);
        vect[3]=Calc.unitVector( Calc.add(vaux, vect[0]) );

        vaux.setX(0); vaux.setY(0); vaux.setZ(-1);
        vect[4]=Calc.unitVector( Calc.add(vaux, vect[0]) );

        vaux=Calc.add(vect[1], vect[0]);
        CalcGeom.product2(vaux, 0.5);
        vect[5]=Calc.unitVector( Calc.add(vect[2], vaux) );

        vaux=Calc.add(vect[2], vect[0]);
        CalcGeom.product2(vaux, 0.5);
        vect[6]=Calc.unitVector( Calc.add(vect[3], vaux) );

        vaux=Calc.add(vect[3], vect[0]);
        CalcGeom.product2(vaux, 0.5);
        vect[7]=Calc.unitVector( Calc.add(vect[4], vaux) );

        vaux=Calc.add(vect[4], vect[0]);
        CalcGeom.product2(vaux, 0.5);
        vect[8]=Calc.unitVector( Calc.add(vect[1], vaux) );

        return vect;
    }

    /**
     * precalculates the places for metal ion acceptors
     *
     * @return
     */
    private static Atom[] calcMultMetActVectors() {

        Atom [] l_rays=new Atom [14];

        //create the 6 vertices of the octahedron
        l_rays[0]=new AtomImpl();
        l_rays[0].setX(0); l_rays[0].setY(1); l_rays[0].setZ(0);

        l_rays[1]=new AtomImpl();
        l_rays[1].setX(1); l_rays[1].setY(0); l_rays[1].setZ(0);

        l_rays[2]=new AtomImpl();
        l_rays[2].setX(0); l_rays[2].setY(0); l_rays[2].setZ(-1);

        l_rays[3]=new AtomImpl();
        l_rays[3].setX(-1); l_rays[3].setY(0); l_rays[3].setZ(0);

        l_rays[4]=new AtomImpl();
        l_rays[4].setX(0); l_rays[4].setY(0); l_rays[4].setZ(1);

        l_rays[5]=new AtomImpl();
        l_rays[5].setX(0); l_rays[5].setY(-1); l_rays[5].setZ(0);

        Atom a=null;
        Atom b=null;
        Atom c=null;

        Atom vaux=new AtomImpl();

        //calc the centers of octahedron faces
        for(int i=0;i<4;i++){

            a=l_rays[0];
            b=l_rays[i+1];
            c=l_rays[((i+1)%4)+1];

            vaux.setX(0); vaux.setY(0); vaux.setZ(0);
            CalcGeom.addEquals(vaux, a);
            CalcGeom.addEquals(vaux, b);
            CalcGeom.addEquals(vaux, c);
            l_rays[5+i*2+1]=Calc.unitVector( CalcGeom.product(vaux, 1.0/3.0));

            a=l_rays[5];
            vaux.setX(0); vaux.setY(0); vaux.setZ(0);
            CalcGeom.addEquals(vaux, a);
            CalcGeom.addEquals(vaux, b);
            CalcGeom.addEquals(vaux, c);
            l_rays[5+i*2+2]=Calc.unitVector( CalcGeom.product(vaux, 1.0/3.0) );
        }

        return l_rays;
    }


    /**
     * checks if the point is in a protein atom different than at
     *
     * @param at
     * @param point
     * @param dist : a value greater than at-point distance
     * @return
     */
    private boolean isValidActivePoint(Mol2Atom at, Atom point, double dist){

        ArrayList<Atom> listaux=new ArrayList<Atom>();

       
        if( dist < (DIST_MAXRAD+0.5) ){
            getNeighAtoms(point, dist, listaux);

            if( listaux.size()>1 )
                return false;//point inside another protein atom
        }
        else{
            if( getBtree().areNeighbours(point, DIST_MAXRAD) )
                return false;//point inside another protein atom
           
        }

        // if the active point is to far from active site
        if( getActSiteBtree()!=null ){
            if( !getActSiteBtree().areNeighbours(point, DIST_POINT_PROBE) )
                return false;
        }
        else{
            //if not exists the bsptree uses active site radii
            double sqlim=Math.pow( getActSiteRadii()+DIST_POINT_PROBE,2);

            if( CalcGeom.squareDistance(point, getActSiteCenter())>sqlim )
                return false;
        }

        return true;
    }

    /**
     * @return the actSiteBtree
     */
    protected BSPTree getActSiteBtree() {
        return actSiteBtree;
    }

    /**
     * @param actSiteBtree the actSiteBtree to set
     */
    protected void setActSiteBtree(BSPTree actSiteBtree) {
        this.actSiteBtree = actSiteBtree;
    }

    /**
     * @return the actSite
     */
    public List<Atom> getActSite() {
        return actSite;
    }

    /**
     * @param actSite the actSite to set
     */
    protected void setActSite(List<Atom> actSite) {
        this.actSite = actSite;
    }

    /**
     * @return the actSiteRadii
     */
    public double getActSiteRadii() {
        return actSiteRadii;
    }

    /**
     * @param actSiteRadii the actSiteRadii to set
     */
    protected void setActSiteRadii(double actSiteRadii) {
        this.actSiteRadii = actSiteRadii;
    }

}
