/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.*;
import org.optiprot.OptiProtParameters;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcIntegrals;
import org.optiprot.potential.element.AtomInfo;
import org.optiprot.potential.element.CharmmResidue;
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
public class MolecularElementsImpl implements IMolecularElements {

    private List<MolecularBond> m_bonds=new ArrayList<MolecularBond>();
    private List<MolecularAngle> m_angles=new ArrayList<MolecularAngle>();
    private List<MolecularDihedral> m_dihedrals=new ArrayList<MolecularDihedral>();
    private List<MolecularImproper> m_impropers=new ArrayList<MolecularImproper>();
    private List<MolecularNonbonded> m_nonbondeds=new ArrayList<MolecularNonbonded>();
    private List<MolecularPair> m_pairs=new ArrayList<MolecularPair>();

    private boolean m_useBonds=true;
    private boolean m_useAngles=true;
    private boolean m_useImpropers=true;
    private boolean m_calcBornRadii=true;
    private boolean m_aproxConnect=true;

    //stores CHARMM topology database
    private HashMap<String,CharmmResidue> m_topology=null;

    ////private HashMap<String,MolecularPair> m_pairsTable=new HashMap<String,MolecularPair>();

    //stores bonded atoms (adjacency matrix) (I-J) (J-I)
    private HashMap<String,List<String>> m_bondsTable=new HashMap<String,List<String>>();

    static final private String TABLE_SEP="_";
    static private double SQLENGTH_13_BOND=2.71*2.71;
    
    private Chain m_chain=null;
    
    //////////////////////////////


    /**
     * create the molecular elements for the potential energy calculation of the
     * given chain
     *
     * @param chain
     * @param grid
     * @param parameters
     */
    public MolecularElementsImpl( Chain chain, AtomGrid grid, OptiProtParameters parameters ) {

        this.setTopology( parameters.getTopology() );

        m_chain=chain;

        IForceField ffield=parameters.getForceField();
        BSPTree btree=new BSPTree( m_chain );
        CalcIntegrals.generateGrid(btree, ffield, grid);

        List<Group> groups = chain.getAtomGroups();

        AtomInfo [] listAtoms=new AtomInfo [btree.getNumAtoms()];
        int cont=0;

        //proccess chain's residues
        for(int i=0;i<groups.size();i++){

            //previous and next group
            Group prevGrp=null;
            Group nextGrp=null;

            if((i-1)>=0)
                prevGrp=groups.get(i-1);

            if((i+1)<groups.size())
                nextGrp=groups.get(i+1);

            createBonds( groups.get(i), i, prevGrp, nextGrp);
            createDihedrals( groups.get(i), i, prevGrp, nextGrp);
            createImpropers( groups.get(i), i, prevGrp, nextGrp);

            //get the amino acid's atoms
            for( Atom at : groups.get(i).getAtoms() ){

                AtomInfo atInfo=new AtomInfo( at, i);

                if( this.isCalcBornRadii() ){
                    atInfo.calcBornRadii( btree, grid, ffield );
                }

                listAtoms[cont++]=atInfo;
            }

        }

        int length=listAtoms.length;

        //proccess pairs
        for( int i=0; i<length; i++){
            for( int j=i;j<length;j++){

                AtomInfo at1=listAtoms[i];
                AtomInfo at2=listAtoms[j];

                String name1 = RotamerLibrary.residue2Charmm(
                        groups.get(at1.getResIdx()).getPDBName() );
                String name2 = RotamerLibrary.residue2Charmm(
                        groups.get(at2.getResIdx()).getPDBName() );

                CharmmResidue res1=this.getTopology().get(name1);
                CharmmResidue res2=this.getTopology().get(name2);

                MolecularNonbonded newpair=new MolecularNonbonded( at1.getAtom(), at2.getAtom());

                Integer typeA=res1.getAtomType(at1.getAtom().getName());
                Integer typeB=res2.getAtomType(at2.getAtom().getName());
                Double chargeA=res1.getAtomCharge(at1.getAtom().getName());
                Double chargeB=res2.getAtomCharge(at2.getAtom().getName());

                if( typeA==null || typeB==null ){
                    continue;
                }

                newpair.setAtomTypeA( typeA );
                newpair.setAtomTypeB( typeB );
                newpair.setChargeA(chargeA);
                newpair.setChargeB(chargeB);

                if( this.isCalcBornRadii() ){
                    newpair.setBornRadiusAB( at1.getBornRadii()*at2.getBornRadii() );
                }

                this.addPair(newpair,at1.getResIdx(), at2.getResIdx());

                if( i!=j )
                    this.addNonbonded(newpair,at1.getResIdx(), at2.getResIdx());
            }
        }

        btree=null;
        listAtoms=null;
    }


    /**
     * evaluates the potential energy of the given chain, without save
     * the molecular elements
     *
     * @param chain
     * @param grid
     * @param parameters
     * @param cMolMechanics : if calcs Molecular Mechanics terms
     * @param cGBorn : if calcs Generalized Born term (solvent polarization)
     * @param cSASA : if calcs SASA term (hydrophobic effect)
     * @return
     */
    static public Double calcEnergy( Chain chain, AtomGrid grid, OptiProtParameters parameters,
            boolean cMolMechanics, boolean cGBorn, boolean cSASA ){

        double energy=0.0;
        double energyBond=0.0;
        double energyDihe=0.0;
        double energyImpr=0.0;
        double energyCoulVDW=0.0;
        double energyGBorn=0.0;
        double energySASA=0.0;

        IForceField ffield=parameters.getForceField();


        if( cSASA && !cMolMechanics && !cGBorn ){
            return MMPotentialEnergy.calcSASA(chain, ffield, parameters);
        }

        BSPTree btree=new BSPTree( chain );

        if( cGBorn ){
            CalcIntegrals.generateGrid(btree, ffield, grid);
        }

        List<Group> groups = chain.getAtomGroups();

        AtomInfo [] listAtoms=new AtomInfo [btree.getNumAtoms()];
        int cont=0;

        //create topology
        HashMap<String,CharmmResidue> topology = new HashMap<String,CharmmResidue>();

        for(CharmmResidue residue : parameters.getTopology() ){
            topology.put( residue.getName().toUpperCase(), residue);
        }

        //proccess chain's residues
        for(int i=0;i<groups.size();i++){

            //previous and next group
            Group prevGrp=null;
            Group nextGrp=null;

            if((i-1)>=0)
                prevGrp=groups.get(i-1);

            if((i+1)<groups.size())
                nextGrp=groups.get(i+1);

            if( cMolMechanics ){
                energyBond += calcBonds( groups.get(i), prevGrp, nextGrp, ffield, topology);
                energyDihe += calcDihedrals( groups.get(i), prevGrp, nextGrp, ffield, topology);
                energyImpr += calcImpropers( groups.get(i), prevGrp, nextGrp, ffield, topology);
            }

            //get the amino acid's atoms
            for( Atom at : groups.get(i).getAtoms() ){

                AtomInfo atInfo=new AtomInfo( at, i);

                if( cGBorn ){
                    atInfo.calcBornRadii( btree, grid, ffield );
                }

                listAtoms[cont++] = atInfo;
            }

        }

        int length=listAtoms.length;
        MolecularNonbonded newpair=new MolecularNonbonded();

        //proccess pairs
        for( int i=0; i<length; i++){
            for( int j=i;j<length;j++){

                AtomInfo at1=listAtoms[i];
                AtomInfo at2=listAtoms[j];

                String name1 = RotamerLibrary.residue2Charmm(
                        groups.get(at1.getResIdx()).getPDBName() );
                String name2 = RotamerLibrary.residue2Charmm(
                        groups.get(at2.getResIdx()).getPDBName() );

                CharmmResidue res1=topology.get(name1);
                CharmmResidue res2=topology.get(name2);

                newpair.setAtoms( at1.getAtom().getCoords(),
                        at2.getAtom().getCoords() );

                Integer typeA=res1.getAtomType(at1.getAtom().getName());
                Integer typeB=res2.getAtomType(at2.getAtom().getName());
                Double chargeA=res1.getAtomCharge(at1.getAtom().getName());
                Double chargeB=res2.getAtomCharge(at2.getAtom().getName());

                if( typeA==null || typeB==null ){
                    continue;
                }

                newpair.setAtomTypeA( typeA );
                newpair.setAtomTypeB( typeB );
                newpair.setChargeA(chargeA);
                newpair.setChargeB(chargeB);

                if( cGBorn ){
                    newpair.setBornRadiusAB( at1.getBornRadii()*at2.getBornRadii() );
                    energyGBorn += calcGBorn(newpair, ffield);
                }


                if( i!=j && cMolMechanics )
                    energyCoulVDW += calcCoulombVDW(newpair, ffield);

            }
        }

        if( cSASA ){
            energySASA = MMPotentialEnergy.calcSASA(chain, ffield, parameters);
        }

        energy=energyBond + energyDihe + energyImpr + energyCoulVDW +
                energyGBorn + energySASA;

        btree=null;
        listAtoms=null;

        return energy;
    }


    public Chain getChain() {
        return m_chain;
    }

    public boolean hasAngles() {

        if( !isUseAngles() )
            return false;

        return !m_angles.isEmpty();
    }

    public boolean hasBonds() {

        if( !isUseBonds() )
            return false;

        return !m_bonds.isEmpty();
    }

    public boolean hasDihedrals() {
        return !m_dihedrals.isEmpty();
    }

    public boolean hasImpropers() {

        if( !isUseImpropers() )
            return false;

        return !m_impropers.isEmpty();
    }

    public boolean hasNonbondeds() {
        return !m_nonbondeds.isEmpty();
    }

    public boolean hasPairs() {
        return !m_pairs.isEmpty();
    }

    public List<MolecularAngle> getAngles() {
        return m_angles;
    }

    public List<MolecularBond> getBonds() {
        return m_bonds;
    }

    public List<MolecularDihedral> getDihedrals() {
        return m_dihedrals;
    }

    public List<MolecularImproper> getImpropers() {
        return m_impropers;
    }

    public List<MolecularNonbonded> getNonbondeds() {
        return m_nonbondeds;
    }

    public List<MolecularPair> getPairs() {
        return m_pairs;
    }

    /**
     * @return the m_topology
     */
    protected HashMap<String,CharmmResidue> getTopology() {
        return m_topology;
    }

    /**
     * @param m_topology the m_topology to set
     */
    public void setTopology(List<CharmmResidue> topology) {
        this.m_topology = new HashMap<String,CharmmResidue>();

        for(CharmmResidue residue : topology){
            this.m_topology.put( residue.getName().toUpperCase(), residue);
        }
    }

    /**
     * adds the pair to the pairs' list and hastable (dont check if exists)
     *
     * @param pair
     * @param idxA
     * @param idxB
     */
    private void addPair(MolecularPair pair, int idxA, int idxB) {

////        String key=makePairsKey(idxA, pair.getNameAtomA(), idxB, pair.getNameAtomB());
////        this.getPairsTable().put(key, pair);
        
        this.getPairs().add(pair);

    }

////    protected MolecularPair getPair(Atom atA, int idxA, Atom atB, int idxB) {
////
////        MolecularPair pair=null;
////
////        String key=makePairsKey( idxB, atB.getName(), idxA, atA.getName() );
////        pair=this.getPairsTable().get(key);
////
////        if( pair!=null )
////            return pair;
////
////        key=makePairsKey(idxA, atA.getName(), idxB, atB.getName());
////        pair=this.getPairsTable().get(key);
////
////        return pair;
////    }

    /**
     * adds the bond to the bonds' list and hastable (dont check if exists)
     * @param bond
     * @param idxA
     * @param idxB
     */
    private void addBond(MolecularBond bond, int idxA, int idxB) {

        if( !this.isAproxConnect() ){
            String key1=makePairsKey(idxA, bond.getNameAtomA());
            String key2=makePairsKey(idxB, bond.getNameAtomB());

            if( !this.getBondsTable().containsKey(key1) ){
                this.getBondsTable().put(key1, new ArrayList<String>());
            }

            this.getBondsTable().get(key1).add(key2);

            if( !this.getBondsTable().containsKey(key2) ){
                this.getBondsTable().put(key2, new ArrayList<String>());
            }

            this.getBondsTable().get(key2).add(key1);
        }
        
        this.getBonds().add(bond);
    }


    /**
     * adds the nonbonded to the nonbondeds' list
     * checks if the pair is a bond or a 1,3 atoms (separated by 2 bonds)
     *
     * @param pair
     * @param idxA
     * @param idxB
     */
    private void addNonbonded(MolecularNonbonded pair, int idxA, int idxB) {

        if( !this.isAproxConnect() ){
            String key1=makePairsKey(idxA, pair.getNameAtomA());
            String key2=makePairsKey(idxB, pair.getNameAtomB());

            // exclude bondeds
            if( existsBond( key1, key2 ) )
                return;

            // 1,3 bondeds exclude
            if( exists13Bond( key1, key2 ) )
                return;
        }
        else{
            if( is13Bond(pair) )
                return;
        }

        this.getNonbondeds().add(pair);

    }

    /**
     * aproximation for 13 bond detection based in distance
     *
     * @param pair
     * @return
     */
    private boolean is13Bond(MolecularNonbonded pair) {

        if( pair.getSqrDistance()<SQLENGTH_13_BOND )
            return true;

        return false;

    }


    /**
     * check in the bonds table if exists a bond with the given keys
     * 
     * @param key1
     * @param key2
     * @return
     */
    private boolean existsBond( String key1, String key2 ){

       if( !this.getBondsTable().containsKey(key1) )
            return false;

       for( String str : this.getBondsTable().get(key1) ){
           if( str.equals(key2) )
               return true;
       }

       return false;
    }

    /**
     * check if the atoms are separated by two bonds
     *
     * @param key1
     * @param key2
     * @return
     */
    private boolean exists13Bond( String key1, String key2 ){

        if( !this.getBondsTable().containsKey(key1) )
            return false;

       for( String str1 : this.getBondsTable().get(key1) ){

           if( !this.getBondsTable().containsKey(str1) )
                continue;

           for( String str2 : this.getBondsTable().get(str1) ){
               if( str2.equals(key2) )
                    return true;
           }
       }

       return false;
    }

    /**
     * create the bonds specified in the residue topology and adds them to
     * the bonds' list
     *
     * @param group : actual residue
     * @param idx : index of actual residue
     * @param prevGrp : previous residue
     * @param nextGrp : next residue
     */
    private void createBonds(Group group, int idx, Group prevGrp, Group nextGrp) {

        String name = RotamerLibrary.residue2Charmm(
                group.getPDBName().toUpperCase() );
        String prevName="";
        String nextName="";

        if( prevGrp!=null)
            prevName = RotamerLibrary.residue2Charmm(
                    prevGrp.getPDBName().toUpperCase() );

        if( nextGrp!=null)
            nextName = RotamerLibrary.residue2Charmm(
                    nextGrp.getPDBName().toUpperCase() );

        
        CharmmResidue res=this.getTopology().get(name);
        CharmmResidue prevRes=this.getTopology().get(prevName);
        CharmmResidue nextRes=this.getTopology().get(nextName);

        //proccess all bond's in residue topology
        for(MolecularBond bond : res.getBonds() ){

            String nameA=bond.getNameAtomA();
            String nameB=bond.getNameAtomB();

            if( nameA.charAt(0)=='-' && prevRes==null )
                continue;
            if( nameB.charAt(0)=='+' && nextRes==null )
                continue;

            Integer typeA=res.getAtomType(nameA);
            Integer typeB=res.getAtomType(nameB);

            Atom atA=null;
            Atom atB=null;
            int idxA=idx;
            int idxB=idx;

            try{
                if( nameA.charAt(0)=='-'  ){
                    atA=prevGrp.getAtom(nameA.substring(1));
                    idxA=idx-1;
                    typeA=prevRes.getAtomType(nameA.substring(1));
                }
                else{
                    atA=group.getAtom(nameA);
                }

                if( nameB.charAt(0)=='+' ){
                    atB=nextGrp.getAtom(nameB.substring(1));
                    idxB=idx+1;
                    typeB=nextRes.getAtomType(nameB.substring(1));
                }
                else{
                    atB=group.getAtom(nameB);
                }
            }
            catch( StructureException e ){
                atA=null; atB=null;
            }

            if( atA==null || atB==null )
                continue;

            MolecularBond newbond=new MolecularBond(atA,atB);
            newbond.setAtomTypeA(typeA);
            newbond.setAtomTypeB(typeB);

            this.addBond(newbond, idxA, idxB);
        }
    }

////    /**
////     * @return the m_pairsTable
////     */
////    protected HashMap<String, MolecularPair> getPairsTable() {
////        return m_pairsTable;
////    }
////
////
////    /**
////     * @param m_pairsTable the m_pairsTable to set
////     */
////    private void setPairsTable(HashMap<String, MolecularPair> m_pairsTable) {
////        this.m_pairsTable = m_pairsTable;
////    }

    /**
     * @return the m_bondsTable
     */
    protected HashMap<String,List<String>> getBondsTable() {
        return m_bondsTable;
    }


    /**
     * @param m_bondsTable the m_bondsTable to set
     */
    private void setBondsTable(HashMap<String,List<String>> m_bondsTable) {
        this.m_bondsTable = m_bondsTable;
    }

    /**
     * construct a key for the pairs table
     *
     * @param resA
     * @param nameA
     * @param resB
     * @param nameB
     * @return
     */
    static private String makePairsKey(int resA, String nameA, int resB, String nameB ){

        return resA+TABLE_SEP+nameA+TABLE_SEP+resB+TABLE_SEP+nameB;
    }

    static private String makePairsKey(int resA, String nameA ){

        return resA+TABLE_SEP+nameA;
    }

    /**
     * create the dihedrals and adds them to the dihedrals' list
     *
     * @param group
     * @param idx : index of actual residue
     * @param prevGrp
     * @param nextGrp
     */
    private void createDihedrals(Group group, int idx, Group prevGrp, Group nextGrp) {

        MolecularDihedral dihedral=null;

        String name = RotamerLibrary.residue2Charmm(
                group.getPDBName().toUpperCase() );
        CharmmResidue res=this.getTopology().get(name);

        // create phi/psi dihedral

        //phi :  -C N CA C
        if( prevGrp!=null ){
            try {
                dihedral=new MolecularDihedral();

                dihedral.setAtom(0, prevGrp.getAtom("C"), res.getAtomCharge("C"), res.getAtomType("C"));
                dihedral.setAtom(1, group.getAtom("N"), res.getAtomCharge("N"), res.getAtomType("N"));
                dihedral.setAtom(2, group.getAtom("CA"), res.getAtomCharge("CA"), res.getAtomType("CA"));
                dihedral.setAtom(3, group.getAtom("C"), res.getAtomCharge("C"), res.getAtomType("C"));
                dihedral.calcDihedralAng();

                this.getDihedrals().add(dihedral);
            } catch (StructureException ex) {
                Logger.getLogger(MolecularElementsImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        //psi : N CA C +N
        if( nextGrp!=null ){
            try {
                dihedral=new MolecularDihedral();

                dihedral.setAtom(0, group.getAtom("N"), res.getAtomCharge("N"), res.getAtomType("N"));
                dihedral.setAtom(1, group.getAtom("CA"), res.getAtomCharge("CA"), res.getAtomType("CA"));
                dihedral.setAtom(2, group.getAtom("C"), res.getAtomCharge("C"), res.getAtomType("C"));
                dihedral.setAtom(3, nextGrp.getAtom("N"), res.getAtomCharge("N"), res.getAtomType("N"));
                dihedral.calcDihedralAng();

                this.getDihedrals().add(dihedral);
            } catch (StructureException ex) {
                Logger.getLogger(MolecularElementsImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }

    /**
     * create the impropers and adds them to the impropers' list
     *
     * @param group
     * @param idx : index of actual residue
     * @param prevGrp
     * @param nextGrp
     */
    private void createImpropers(Group group, int idx, Group prevGrp, Group nextGrp) {

        String name = RotamerLibrary.residue2Charmm(
                group.getPDBName().toUpperCase() );
        String prevName="";
        String nextName="";

        if( prevGrp!=null)
            prevName = RotamerLibrary.residue2Charmm(
                    prevGrp.getPDBName().toUpperCase() );

        if( nextGrp!=null)
            nextName = RotamerLibrary.residue2Charmm(
                    nextGrp.getPDBName().toUpperCase() );


        CharmmResidue res=this.getTopology().get(name);
        CharmmResidue prevRes=this.getTopology().get(prevName);
        CharmmResidue nextRes=this.getTopology().get(nextName);

        for( MolecularImproper improper : res.getImpropers() ){

            try{

                Atom at1=group.getAtom(improper.getAtom(0).getName());
                Atom at2=null;
                Atom at3=null;
                Atom at4=group.getAtom(improper.getAtom(3).getName());

                String name2=improper.getAtom(1).getName();
                String name3=improper.getAtom(2).getName();

                if( name2.charAt(0)=='+'){
                    if( nextGrp==null )
                        continue;
                    at2=nextGrp.getAtom(name2.substring(1));
                }
                else if( name2.charAt(0)=='-' ){
                    if( prevGrp==null )
                        continue;
                    at2=prevGrp.getAtom(name2.substring(1));
                }
                else{
                    at2=group.getAtom(name2);
                }

                if( name3.charAt(0)=='+'){
                    if( nextGrp==null )
                        continue;
                    at3=nextGrp.getAtom(name3.substring(1));
                }
                else if( name3.charAt(0)=='-' ){
                    if( prevGrp==null )
                        continue;
                    at3=prevGrp.getAtom(name3.substring(1));
                }
                else{
                    at3=group.getAtom(name3);
                }

                MolecularImproper newimproper=new MolecularImproper();

                newimproper.setAtom(0, at1, res.getAtomCharge(at1.getName()), res.getAtomType(at1.getName()));
                newimproper.setAtom(1, at2, res.getAtomCharge(at2.getName()), res.getAtomType(at2.getName()));
                newimproper.setAtom(2, at3, res.getAtomCharge(at3.getName()), res.getAtomType(at3.getName()));
                newimproper.setAtom(3, at4, res.getAtomCharge(at4.getName()), res.getAtomType(at4.getName()));
                newimproper.calcDihedralAng();

                this.getImpropers().add(newimproper);

            } catch (StructureException ex) {
                ///Logger.getLogger(MolecularElementsImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    /**
     * @return the m_useBonds
     */
    public boolean isUseBonds() {
        return m_useBonds;
    }

    /**
     * @param m_useBonds the m_useBonds to set
     */
    public void setUseBonds(boolean m_useBonds) {
        this.m_useBonds = m_useBonds;
    }

    /**
     * @return the m_useAngles
     */
    public boolean isUseAngles() {
        return m_useAngles;
    }

    /**
     * @param m_useAngles the m_useAngles to set
     */
    public void setUseAngles(boolean m_useAngles) {
        this.m_useAngles = m_useAngles;
    }

    /**
     * @return the m_useImpropers
     */
    public boolean isUseImpropers() {
        return m_useImpropers;
    }

    /**
     * @param m_useImpropers the m_useImpropers to set
     */
    public void setUseImpropers(boolean m_useImpropers) {
        this.m_useImpropers = m_useImpropers;
    }

     /**
     * @return the m_calcBornRadii
     */
    public boolean isCalcBornRadii() {
        return m_calcBornRadii;
    }

    /**
     * @param m_calcBornRadii the m_calcBornRadii to set
     */
    public void setCalcBornRadii(boolean m_calcBornRadii) {
        this.m_calcBornRadii = m_calcBornRadii;
    }

    /**
     * @return the m_aproxConnect
     */
    public boolean isAproxConnect() {
        return m_aproxConnect;
    }

    /**
     * @param m_aproxConnect the m_aproxConnect to set
     */
    public void setAproxConnect(boolean m_aproxConnect) {
        this.m_aproxConnect = m_aproxConnect;
    }


    /**
     * calc. potential energy
     * @param group
     * @param prevGrp
     * @param nextGrp
     * @param ffield
     * @param topology
     * @return
     */
    static private double calcBonds(Group group, Group prevGrp, Group nextGrp, 
            IForceField ffield, HashMap<String,CharmmResidue> topology) {

        double energy=0;

        String name = RotamerLibrary.residue2Charmm(
                group.getPDBName().toUpperCase() );
        String prevName="";
        String nextName="";

        if( prevGrp!=null)
            prevName = RotamerLibrary.residue2Charmm(
                    prevGrp.getPDBName().toUpperCase() );

        if( nextGrp!=null)
            nextName = RotamerLibrary.residue2Charmm(
                    nextGrp.getPDBName().toUpperCase() );


        CharmmResidue res=topology.get(name);
        CharmmResidue prevRes=topology.get(prevName);
        CharmmResidue nextRes=topology.get(nextName);

        MolecularBond newbond=new MolecularBond();

        //proccess all bond's in residue topology
        for(MolecularBond bond : res.getBonds() ){

            String nameA=bond.getNameAtomA();
            String nameB=bond.getNameAtomB();

            if( nameA.charAt(0)=='-' && prevRes==null )
                continue;
            if( nameB.charAt(0)=='+' && nextRes==null )
                continue;

            Integer typeA=res.getAtomType(nameA);
            Integer typeB=res.getAtomType(nameB);

            Atom atA=null;
            Atom atB=null;
            try{
                if( nameA.charAt(0)=='-'  ){
                    atA=prevGrp.getAtom(nameA.substring(1));
                    typeA=prevRes.getAtomType(nameA.substring(1));
                }
                else{
                    atA=group.getAtom(nameA);
                }

                if( nameB.charAt(0)=='+' ){
                    atB=nextGrp.getAtom(nameB.substring(1));
                    typeB=nextRes.getAtomType(nameB.substring(1));
                }
                else{
                    atB=group.getAtom(nameB);
                }
            }
            catch( StructureException e ){
                atA=null; atB=null;
            }

            if( atA==null || atB==null )
                continue;

            newbond.setAtoms( atA.getCoords(), atB.getCoords() );
            newbond.setAtomTypeA(typeA);
            newbond.setAtomTypeB(typeB);

            energy+=MMPotentialEnergy.calcBondLengthEnergy(newbond, ffield);
        }

        return energy;
    }

    /**
     * calc. potential energy
     * @param group
     * @param prevGrp
     * @param nextGrp
     * @param ffield
     * @param topology
     * @return
     */
    static private double calcDihedrals(Group group, Group prevGrp, Group nextGrp, 
            IForceField ffield, HashMap<String,CharmmResidue> topology) {

        double energy=0;

        MolecularDihedral dihedral=new MolecularDihedral();

        String name = RotamerLibrary.residue2Charmm(
                group.getPDBName().toUpperCase() );
        CharmmResidue res=topology.get(name);

        // create phi/psi dihedral

        //phi :  -C N CA C
        if( prevGrp!=null ){
            try {

                dihedral.setAtom(0, prevGrp.getAtom("C"), res.getAtomCharge("C"), res.getAtomType("C"));
                dihedral.setAtom(1, group.getAtom("N"), res.getAtomCharge("N"), res.getAtomType("N"));
                dihedral.setAtom(2, group.getAtom("CA"), res.getAtomCharge("CA"), res.getAtomType("CA"));
                dihedral.setAtom(3, group.getAtom("C"), res.getAtomCharge("C"), res.getAtomType("C"));
                dihedral.calcDihedralAng();

                energy+=MMPotentialEnergy.calcTorsionEnergy(dihedral, ffield);
            } catch (StructureException ex) {
                Logger.getLogger(MolecularElementsImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        //psi : N CA C +N
        if( nextGrp!=null ){
            try {

                dihedral.setAtom(0, group.getAtom("N"), res.getAtomCharge("N"), res.getAtomType("N"));
                dihedral.setAtom(1, group.getAtom("CA"), res.getAtomCharge("CA"), res.getAtomType("CA"));
                dihedral.setAtom(2, group.getAtom("C"), res.getAtomCharge("C"), res.getAtomType("C"));
                dihedral.setAtom(3, nextGrp.getAtom("N"), res.getAtomCharge("N"), res.getAtomType("N"));
                dihedral.calcDihedralAng();

                energy+=MMPotentialEnergy.calcTorsionEnergy(dihedral, ffield);
            } catch (StructureException ex) {
                Logger.getLogger(MolecularElementsImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        return energy;
    }

    /**
     * calc. potential energy
     * @param group
     * @param prevGrp
     * @param nextGrp
     * @param ffield
     * @param topology
     * @return
     */
    static private double calcImpropers(Group group, Group prevGrp, Group nextGrp, 
            IForceField ffield, HashMap<String,CharmmResidue> topology ) {

        double energy=0;

        String name = RotamerLibrary.residue2Charmm(
                group.getPDBName().toUpperCase() );
        String prevName="";
        String nextName="";

        if( prevGrp!=null)
            prevName = RotamerLibrary.residue2Charmm(
                    prevGrp.getPDBName().toUpperCase() );

        if( nextGrp!=null)
            nextName = RotamerLibrary.residue2Charmm(
                    nextGrp.getPDBName().toUpperCase() );


        CharmmResidue res=topology.get(name);
        MolecularImproper newimproper=new MolecularImproper();

        for( MolecularImproper improper : res.getImpropers() ){

            try{

                Atom at1=group.getAtom(improper.getAtom(0).getName());
                Atom at2=null;
                Atom at3=null;
                Atom at4=group.getAtom(improper.getAtom(3).getName());

                String name2=improper.getAtom(1).getName();
                String name3=improper.getAtom(2).getName();

                if( name2.charAt(0)=='+'){
                    if( nextGrp==null )
                        continue;
                    at2=nextGrp.getAtom(name2.substring(1));
                }
                else if( name2.charAt(0)=='-' ){
                    if( prevGrp==null )
                        continue;
                    at2=prevGrp.getAtom(name2.substring(1));
                }
                else{
                    at2=group.getAtom(name2);
                }

                if( name3.charAt(0)=='+'){
                    if( nextGrp==null )
                        continue;
                    at3=nextGrp.getAtom(name3.substring(1));
                }
                else if( name3.charAt(0)=='-' ){
                    if( prevGrp==null )
                        continue;
                    at3=prevGrp.getAtom(name3.substring(1));
                }
                else{
                    at3=group.getAtom(name3);
                }

                

                newimproper.setAtom(0, at1, res.getAtomCharge(at1.getName()), res.getAtomType(at1.getName()));
                newimproper.setAtom(1, at2, res.getAtomCharge(at2.getName()), res.getAtomType(at2.getName()));
                newimproper.setAtom(2, at3, res.getAtomCharge(at3.getName()), res.getAtomType(at3.getName()));
                newimproper.setAtom(3, at4, res.getAtomCharge(at4.getName()), res.getAtomType(at4.getName()));
                newimproper.calcDihedralAng();

                energy+=MMPotentialEnergy.calcTorsionEnergy(newimproper, ffield);


            } catch (StructureException ex) {
                ///Logger.getLogger(MolecularElementsImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        return energy;
    }

    /**
     * calc. potential energy
     * @param pair
     * @param ffield
     * @return
     */
    static private double calcGBorn( MolecularPair pair, IForceField ffield ){
        return MMPotentialEnergy.calcGBornSolvationEnergy( pair, ffield );
    }

    /**
     * calc. potential energy
     * @param nonbonded
     * @param ffield
     * @return
     */
    static private double calcCoulombVDW( MolecularNonbonded nonbonded, IForceField ffield ){

        if( nonbonded.getSqrDistance()<SQLENGTH_13_BOND )
            return 0.0;

        return MMPotentialEnergy.calcCoulombVDWEnergy(nonbonded, ffield);
    }
}
