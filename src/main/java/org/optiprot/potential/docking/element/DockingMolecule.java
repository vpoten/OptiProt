/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.optiprot.potential.docking.cdk.MMFF94potential;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.modeling.builder3d.ForceFieldConfigurator;
import org.optiprot.io.mol2.Mol2Atom;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcTransform;
import org.optiprot.potential.docking.element.DockingAtomClassify.simpletype;

/**
 *
 * @author victor
 */
public class DockingMolecule {

    private IMolecule molecule=null;//CDK
    private MMFF94potential fffunction=null;

    private BSPTree btree=null;
    private Chain chain = null;//biojava

    private Double internalEnergy=null;

    private int numRotableBonds = 0;

    //auxiliary list
    private ArrayList<Atom> listAtom=new ArrayList<Atom>();

    static private Double _J_2_cal = 1.0/4.184;//J to cal

    //OpenBabel tools tokens
    private static String OB_NUMROTBOND="Number of rotatable bonds:";
    private static String OB_TOTENERGY="TOTAL ENERGY =";


    ///////////////////////

    public DockingMolecule( Chain chain ) {
        this.setChain(chain);
        this.setBtree( new BSPTree(getChain()) );
    }

    
    /**
     * @return the btree
     */
    protected BSPTree getBtree() {
        return btree;
    }

    /**
     * @param btree the btree to set
     */
    protected void setBtree(BSPTree btree) {
        this.btree = btree;
    }

    /**
     * @return the chain
     */
    public Chain getChain() {
        return chain;
    }

    /**
     * @param chain the chain to set
     */
    protected void setChain(Chain chain) {
        this.chain = chain;
    }



    /**
     * @return the molecule
     */
    protected IMolecule getMolecule() {
        return molecule;
    }

    /**
     * @param molecule the molecule to set
     */
    protected void setMolecule(IMolecule molecule) {
        this.molecule = molecule;
    }


    /**
     * classify the atom by simpletype
     *
     * @param atP
     * @return
     */
    public simpletype getSimpletype(Atom atM) {
        return DockingAtomClassify.getSimpletype( (Mol2Atom)atM );
    }

    /**
     * @return the listAtom
     */
    protected ArrayList<Atom> getListAtomAux() {
        return listAtom;
    }

    /**
     * gets the internal potential energy using a ffield
     *
     * @return energy in kcal/mol
     */
    public Double getInternalEnergy(){

        if( internalEnergy==null && getMolecule()!=null && fffunction!=null ){
            setInternalEnergy( fffunction.energyFunctionOfAMolecule(getMolecule())*_J_2_cal );
        }

        return this.internalEnergy;
    }

    /**
     * @param internalEnergy the internalEnergy to set
     */
    private void setInternalEnergy(Double internalEnergy) {
        this.internalEnergy = internalEnergy;
    }

    public void setInternalEnergy(double internalEnergy) {
        this.internalEnergy = internalEnergy;
    }


    
    /**
     * read molecule from file
     *
     * @param filename : sybyl mol2 file
     * @return
     */
    public void readCDKMol2(String filename){

        IMolecule mol=null;

        try {
            Mol2Reader reader = new Mol2Reader(new BufferedReader(new FileReader(filename)));
            mol = (IMolecule) reader.read(new Molecule());
        } catch (FileNotFoundException ex) {
            setMolecule(null);
            return;
        } catch (CDKException ex) {
            setMolecule(null);
            return;
        } catch (Exception ex) {
            setMolecule(null);
            return;
        }

        setMolecule(mol);

        ForceFieldConfigurator ffconfig=new ForceFieldConfigurator();
        
        try {
            ffconfig.setForceFieldConfigurator("mmff94");
            ffconfig.assignAtomTyps(getMolecule());
            fffunction = new MMFF94potential(getMolecule(), ffconfig.getParameterSet());
        } catch (Exception ex) {
            fffunction = null;
        }

        
    }

    
    /**
     * gets the CoM of the molecule
     *
     * @return
     */
    public Atom getCoM(){
        return CalcGeom.getCentroid(getChain());
    }


    /**
     * shifts the atoms of the molecule
     *
     * @param shiftVector
     */
    public void shift( Atom shiftVector ){

        CalcTransform.translateSubChain(shiftVector, getChain(), 0, getChain().getAtomLength());
    }

    /**
     * fills the list with the atoms inside the box given by center and side length

     * @param center center of the box
     * @param edge_len length of the box edge
     * @param list : I/O
     */
    public void getAtomsInsideBox( Atom center, double edge_len, Collection<Atom> list ){

        double [] x=new double [] { center.getX()-edge_len*0.5, center.getX()+edge_len*0.5};
        double [] y=new double [] { center.getY()-edge_len*0.5, center.getY()+edge_len*0.5};
        double [] z=new double [] { center.getZ()-edge_len*0.5, center.getZ()+edge_len*0.5};

        for( Group grp : getChain().getAtomGroups() ){
            for( Atom at : grp.getAtoms() ){

                if( at.getX()>=x[0] && at.getX()<=x[1] &&
                    at.getY()>=y[0] && at.getY()<=y[1] &&
                    at.getZ()>=z[0] && at.getZ()<=z[1] ){
                    list.add(at);
                }
            }
        }
    }

    /**
     * fills the list with the high active atoms: H-bond donor(H), acceptor
     * donor_accept and metal ions
     *
     * @param list : I/O
     */
    public void getActiveAtoms( Collection<DockingActivePoint> list ){

        for( Group grp : getChain().getAtomGroups() ){
            for( Atom at : grp.getAtoms() ){

                DockingActivePoint ap=isActiveAtom( (Mol2Atom) at );

                if( ap!=null )
                    list.add( ap );
            }
        }
    }

    /**
     *
     * @param at2
     * @return a DockingAtivePoint if is an active atom, null otherwise
     */
    public static DockingActivePoint isActiveAtom( Mol2Atom at2 ){

        simpletype type=DockingAtomClassify.getSimpletype(at2);

        if( type==simpletype.metal || type==simpletype.acceptor ||
                type==simpletype.donor_accept || type==simpletype.donorh ){
            return new DockingActivePoint( at2, type);
        }

        return null;
    }

    /**
     *
     * @param at2
     * @return a DockingAtivePoint if is a nonpolar atom (not H), null otherwise
     */
    public static DockingActivePoint isNonpolarAtom( Mol2Atom at2 ){

        simpletype type=DockingAtomClassify.getSimpletype(at2);

        if( type==simpletype.nonpolar ){
            return new DockingActivePoint( at2, type);
        }

        return null;
    }

    /**
     * calculates the radii of the sphere centered at CoM that contains the atoms
     *
     * @return
     */
    public double calcRadii(){

        Atom com=getCoM();
        double sqmaxdist=-1;

        for( Group grp : getChain().getAtomGroups() ){
            for( Atom at : grp.getAtoms() ){

                double dist=CalcGeom.squareDistance(at, com);

                if( dist > sqmaxdist )
                    sqmaxdist=dist;
            }
        }

        return Math.sqrt(sqmaxdist);
    }

    /**
     * calculates some properties using OpenBabel commandline tools
     * 
     * - internal energy with obenergy -ff MMFF94 <file>
     * - num of rotable bonds with obrotamer <file>
     *
     * @param filepath
     */
    public void calcObProperties( String filepath ){

        Runtime run=Runtime.getRuntime();

        String command="obrotamer "+filepath;

        try {
            Process pro = run.exec(command);
            pro.waitFor();

            BufferedReader r = new BufferedReader( new InputStreamReader(pro.getErrorStream()) );

            String l_line=r.readLine();
            l_line=l_line.trim();

            if( l_line.startsWith(OB_NUMROTBOND) ){
                this.setNumRotableBonds(
                        Integer.parseInt( l_line.substring(OB_NUMROTBOND.length()).trim() )
                        );
            }

            r.close();

        } catch (IOException ex) {
            Logger.getLogger(DockingMolecule.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(DockingMolecule.class.getName()).log(Level.SEVERE, null, ex);
        }

        command="obenergy -ff MMFF94 "+filepath;

        try {
            Process pro = run.exec(command);
            pro.waitFor();

            BufferedReader r = new BufferedReader( new InputStreamReader(pro.getInputStream()) );

            String l_line=r.readLine();

            while( l_line!=null ){

                l_line=l_line.trim();

                if( l_line.startsWith(OB_TOTENERGY) ){

                    l_line=l_line.substring(OB_TOTENERGY.length()).trim();

                    this.setInternalEnergy(
                            Double.parseDouble( l_line.substring(0, l_line.indexOf(" ")) )
                            );

                    break;
                }

                l_line=r.readLine();
            }

            r.close();

        } catch (IOException ex) {
            Logger.getLogger(DockingMolecule.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(DockingMolecule.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     *
     * @return
     */
    public int getNumRotableBonds(){
        return numRotableBonds;
    }

    /**
     * @param numRotableBonds the numRotableBonds to set
     */
    public void setNumRotableBonds(int numRotableBonds) {
        this.numRotableBonds = numRotableBonds;
    }
}
