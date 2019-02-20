/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.rotamer;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.optiprot.io.BioJavaStructureReader;
import org.biojava.bio.structure.*;
import org.optiprot.maths.CalcSurfVol;
import org.optiprot.potential.IForceField;
import org.optiprot.potential.TestForceField;

/**
 *
 * @author victor
 */
public class RotamerLibrary {

    //table of rotamers indexed by AA name
    private HashMap<String, Group []> AATable=new HashMap<String, Group []>();

    //table of rotamer angles indexed by AA name
    private HashMap<String, double [][]> AnglesTable=new HashMap<String, double [][]>();

    //table of SA and VDW surface area of rotamers indexed by AA name
    private HashMap<String, double []> SAAreaTable=new HashMap<String, double []>();
    private HashMap<String, double []> VDWAreaTable=new HashMap<String, double []>();
    private HashMap<String, double []> SideVDWAreaTable=new HashMap<String, double []>();
    //Average SA and VdW
    private HashMap<String, Double> AvgSAAreaTable=new HashMap<String, Double>();
    private HashMap<String, Double> AvgVDWAreaTable=new HashMap<String, Double>();
    private HashMap<String, Double> AvgSideVDWAreaTable=new HashMap<String, Double>();

    //table of sidechain's CoM of rotamers indexed by AA name
    private HashMap<String, Atom []> SidechainCoMTable=new HashMap<String, Atom []>();

    //table of sidechain atoms
    private HashMap<String, Atom [][]> SidechainAtomTable=new HashMap<String, Atom [][]>();

    //table of number of sidechain angles indexed by AA name
    private HashMap<String,Integer> NumAnglesTable=new HashMap<String,Integer>();
    
   
    private static String [] AminoAcids={"ala",
                            "arg",
                            "asn",
                            "asp",
                            "cys",
                            "gln",
                            "glu",
                            "gly",
                            "his",
                            "ile",
                            "leu",
                            "lys",
                            "met",
                            "phe",
                            "pro",
                            "ser",
                            "thr",
                            "trp",
                            "tyr",
                            "val"};

    public RotamerLibrary()
    {
        this.fillNumAngles();

    }

    
    /**
     * fills the table of # angles per aminoacid
     */
    private void fillNumAngles() {
        
        this.getNumAnglesTable().put( AminoAcids[0], 0);
        this.getNumAnglesTable().put( AminoAcids[1], 4);
        this.getNumAnglesTable().put( AminoAcids[2], 2);
        this.getNumAnglesTable().put( AminoAcids[3], 2);
        this.getNumAnglesTable().put( AminoAcids[4], 1);
        this.getNumAnglesTable().put( AminoAcids[5], 3);
        this.getNumAnglesTable().put( AminoAcids[6], 3);
        this.getNumAnglesTable().put( AminoAcids[7], 0);
        this.getNumAnglesTable().put( AminoAcids[8], 2);
        this.getNumAnglesTable().put( AminoAcids[9], 2);
        this.getNumAnglesTable().put( AminoAcids[10], 2);
        this.getNumAnglesTable().put( AminoAcids[11], 4);
        this.getNumAnglesTable().put( AminoAcids[12], 3);
        this.getNumAnglesTable().put( AminoAcids[13], 2);
        this.getNumAnglesTable().put( AminoAcids[14], 1);
        this.getNumAnglesTable().put( AminoAcids[15], 1);
        this.getNumAnglesTable().put( AminoAcids[16], 1);
        this.getNumAnglesTable().put( AminoAcids[17], 2);
        this.getNumAnglesTable().put( AminoAcids[18], 2);
        this.getNumAnglesTable().put( AminoAcids[19], 1);
    }

    private Group get( String cad ){

//        Alanine 	Ala 	A 	nonpolar 	neutral
//        Arginine 	Arg 	R 	polar 	positive
//        Asparagine 	Asn 	N 	polar 	neutral
//        Aspartic acid 	Asp 	D 	polar 	negative
//        Cysteine 	Cys 	C 	nonpolar 	neutral
//        Glutamic acid 	Glu 	E 	polar 	negative
//        Glutamine 	Gln 	Q 	polar 	neutral
//        Glycine 	Gly 	G 	nonpolar 	neutral
//        Histidine 	His 	H 	polar 	positive
//        Isoleucine 	Ile 	I 	nonpolar 	neutral
//        Leucine 	Leu 	L 	nonpolar 	neutral
//        Lysine 	Lys 	K 	polar 	positive
//        Methionine 	Met 	M 	nonpolar 	neutral
//        Phenylalanine 	Phe 	F 	nonpolar 	neutral
//        Proline 	Pro 	P 	nonpolar 	neutral
//        Serine 	Ser 	S 	polar 	neutral
//        Threonine 	Thr 	T 	polar 	neutral
//        Tryptophan 	Trp 	W 	nonpolar 	neutral
//        Tyrosine 	Tyr 	Y 	polar 	neutral
//        Valine 	Val 	V 	nonpolar 	neutral

//        cad=cad.toUpperCase();
//
//        if( cad.equals( "ALA") ){
//        }
//        else if( cad.equals( "ARG") ){
//        }
//        else if( cad.equals( "ASN") ){
//        }
//        else if( cad.equals( "ASP") ){
//        }
//        else if( cad.equals( "CYS") ){
//        }
//        else if( cad.equals( "GLN") ){
//        }
//        else if( cad.equals( "GLU") ){
//        }
//        else if( cad.equals( "GLY") ){
//        }
//        else if( cad.equals( "HIS") ){
//        }
//        else if( cad.equals( "ILE") ){
//        }
//        else if( cad.equals( "LEU") ){
//        }
//        else if( cad.equals( "LYS") ){
//        }
//        else if( cad.equals( "MET") ){
//        }
//        else if( cad.equals( "PHE") ){
//        }
//        else if( cad.equals( "PRO") ){
//        }
//        else if( cad.equals( "SER") ){
//        }
//        else if( cad.equals( "THR") ){
//        }
//        else if( cad.equals( "TRP") ){
//        }
//        else if( cad.equals( "TYR") ){
//        }
//        else if( cad.equals( "VAL") ){
//        }


        return null;
    }

    /**
     * translates the name to the internal amino acid name
     * (3-letter code lower case)
     * 
     * @param aaname
     * @return
     */
    protected static String internalName( String aaname ){

        if( aaname.length()==3 )
            return aaname.toLowerCase();

        return StandardAminoAcid.getAminoAcid(aaname).getPDBName().toLowerCase();
    }

    /**
     * returns a rotamer chose randomly
     *
     * @param aa
     * @return
     */
    public Group getRotamer( String aa ){

        double val=Math.random();

        Group [] array=this.getAATable().get( internalName(aa) );

        return array[ (int)Math.floor(val*(double)array.length) ];
    }


    /**
     * gets the idx of the closer rotamer to aa
     *
     * @param aa
     * @return
     */
    public Integer getRotamerIdx( Group aa ){


        double [] chi=null;

        chi=this.calcSidechainAngles(aa);

        if(chi.length==0)
            return 0;

        return getRotamerIdx( internalName(aa.getPDBName()), chi );
    }


    /**
     * gets the closer rotamer to aa
     *
     * @param aa
     * @return
     */
    public Group getRotamer( AminoAcid aa ){

        return this.getAATable().get(internalName(aa.getPDBName()))[this.getRotamerIdx(aa)];

    }

    /**
     * returns the rotamer index whose angles are closer to angles
     * in array chi
     *
     * @param aa
     * @param chi rotamer angles in degrees
     * @return
     */
    private int getRotamerIdx( String aa, double [] chi ){

        double [][] array_angles=this.getAnglesTable().get(aa);

        double min=java.lang.Double.MAX_VALUE;
        int idx_min=-1;

        for(int i=0;i<array_angles.length;i++){

            double distance=0.0;

            for(int j=0;j<chi.length;j++){

                double ang1=Math.toRadians(array_angles[i][j]);
                double ang2=Math.toRadians(chi[j]);

                double x1=Math.cos(ang1);
                double x2=Math.cos(ang2);
                double y1=Math.sin(ang1);
                double y2=Math.sin(ang2);


                distance+=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
            }

            if( distance<min ){
                idx_min=i;
                min=distance;
            }
        }


        return idx_min;

    }

    /**
     * returns the rotamer whose angles are closer to angles
     * in array chi
     *
     * @param aa
     * @param chi chi angles in degrees
     * @return
     */
    public Group getRotamer( String aa, double [] chi ){

        return this.getAATable().get(internalName(aa))[this.getRotamerIdx( internalName(aa), chi )];

    }

    /**
     * returns the rotamer with a given index
     *
     * @param aa
     * @param idx
     * @return
     */
    public Group getRotamer( String aa, int idx ){

        return this.getAATable().get(internalName(aa))[idx];

    }

    /**
     * returns a rotamer distinct (in angles) than the stored in aa_struct
     *
     * @param aa
     * @param aa_struct
     * @return
     */
    public Group getRotamerDistinct( String aa, Group aa_struct ){

        aa=internalName(aa);
        
        int nangles=this.getNumAngles(aa);

        if( nangles==0 )
            return null;

        List<Atom> rot_atoms=this.getRotamerAtoms(aa_struct, nangles);

        double [] angles=new double[nangles];

        //calculate the angles
        for(int i=0;i<nangles;i++){

            try{
                angles[i]=Calc.torsionAngle( rot_atoms.get(i) ,
                    rot_atoms.get(i+1),
                    rot_atoms.get(i+2),
                    rot_atoms.get(i+3) );

            }
            catch(Exception e){

            }
        }

        int idx_closer=this.getRotamerIdx( aa, angles );

        int length=this.getAATable().get(aa).length;

        //returns a rotamer that is not the closer (chose randomly)
        int idx_distinct=(int)Math.floor(Math.random()*(double)length);

        while( idx_distinct==idx_closer ){

            idx_distinct=(int)Math.floor(Math.random()*(double)length);
        }

        return this.getAATable().get(aa)[idx_distinct];
    }

    /**
     * gets the SASA of the closer rotamer to group
     *
     * @param group
     * @return
     */
    public double getSASA(Group group) {
        
        int idx=this.getRotamerIdx(group);

        return this.getSAAreaTable().get(internalName(group.getPDBName()))[idx];
    }

    public double getAvgSASA(Group group) {
        return this.getAvgSAAreaTable().get(internalName(group.getPDBName()));
    }

    public double getAvgVdWSA(Group group) {
        return this.getAvgVDWAreaTable().get(internalName(group.getPDBName()));
    }

    public double getAvgSideVdWSA(Group group) {
        return this.getAvgSideVDWAreaTable().get(internalName(group.getPDBName()));
    }

    public double getAvgSideVdWSA(String aaname) {
        return this.getAvgSideVDWAreaTable().get( internalName(aaname) );
    }

    /**
     * loads the library of rotamers stored in a directory
     *
     * @param path
     * @param sufix
     * @throws java.lang.Exception
     */
    public void loadLibrary( String path, String sufix ) throws Exception
    {

        for( String aaname : AminoAcids ){

            String path_file=path+File.separator+aaname+sufix;

            Structure struc=BioJavaStructureReader.readStructure( path_file );

            this.indexStructure( aaname, struc );
            
        }

        calcAvgSurfAreas();
    }

    /**
     * @return the AATable
     */
    private HashMap<String, Group []> getAATable() {
        return AATable;
    }


    /**
     * @param AATable the AATable to set
     */
    private void setAATable(HashMap<String,Group []> AATable) {
        this.AATable = AATable;
    }

    /**
     * gets the atoms of the aminoacid used to calculate the angles of
     * the rotamers
     * 
     * @param group
     * @param nangles
     * @return
     */
    private List<Atom> getRotamerAtoms(Group group, int nangles) {

        ArrayList<Atom> result=new ArrayList<Atom>();

        Iterator<Atom> it=group.getAtoms().iterator();

        Atom [] side_chain={null,null,null,null};

        try{

            result.add( group.getAtom("N") );
            result.add( group.getAtom("CA") );
            result.add( group.getAtom("CB") );

            while(it.hasNext()){

                Atom at=it.next();

                if( at.getName().length()==2 && at.getName().charAt(0)!='H' )
                {
                    switch( at.getName().charAt(1) )
                    {
                        case 'G': side_chain[0]=at; break;
                        case 'D':
                            if(nangles>1)
                                side_chain[1]=at;
                            break;
                        case 'E':
                            if(nangles>2)
                                side_chain[2]=at;
                            break;
                        case 'Z':
                            if(nangles>3)
                                side_chain[3]=at;
                            break;

                    }

                }
            }
        }
        catch(Exception e)
        {
            return null;
        }

        for(int i=0;i<4;i++){
            if( side_chain[i]!=null ){
                result.add( side_chain[i] );
            }
        }

        if( (result.size()-3)==nangles )
        {
            return result;
        }

        char last=' ';

        switch(nangles){
            case 1: last='G';break;
            case 2: last='D';break;
            case 3: last='E';break;
            case 4: last='Z';break;
        }

        try{

            if( group.hasAtom( "C"+last+"1") ){
                result.add( group.getAtom("C"+last+"1") );
            }
            else if( group.hasAtom( "C"+last+"2") ){
                result.add( group.getAtom("C"+last+"2") );
            }
            else if( group.hasAtom( "N"+last+"1") ){
                result.add( group.getAtom("N"+last+"1") );
            }
            else if( group.hasAtom( "N"+last+"2") ){
                result.add( group.getAtom("N"+last+"2") );
            }
            else if( group.hasAtom( "O"+last+"1") ){
                result.add( group.getAtom("O"+last+"1") );
            }
            else if( group.hasAtom( "O"+last+"2") ){
                result.add( group.getAtom("O"+last+"2") );
            }

        }
        catch(Exception e)
        {

        }

        return result;

        
    }

    

    /**
     * fills the HashMaps used to index the rotamers by AA name and calculates
     * the angles and some properties of the different rotamers
     * 
     * @param aa
     * @param struct
     */
    private void indexStructure( String aa, Structure struct )
    {
        Group [] array=null;

        if( struct.size()==0){
            //case of GLY
            AminoAcid gly=StandardAminoAcid.getAminoAcid("GLY");
            array = new Group [] {gly};
        }
        else{
             List<Group> groups=struct.getChain(0).getAtomGroups();
             array=groups.toArray(new Group [groups.size()]);
        }

        //fill the HashMap of aminoacids (structures) indexed by AA names
        

        //normalize the atoms' names
        for(int i=0;i<array.length;i++){
            normalizeNames( array[i] );
        }

        this.getAATable().put(aa, array);

        IForceField ff=new TestForceField();
        Atom [] atomArr=null;
        
        if( array.length<=1 ){
            //case of ALA and GLY

            this.getAnglesTable().put( aa, new double[0][0]);

            try {

                atomArr=getSideChainAtoms((AminoAcid)array[0]);

                this.getSidechainCoMTable().put( aa,
                    new Atom [] { Calc.getCentroid( atomArr ) } );

                this.getSidechainAtomTable().put(aa,
                        new Atom [][] {atomArr} );

                this.getSideVDWAreaTable().put(aa,
                        new double[] {CalcSurfVol.VDWsurface(atomArr, ff, 4)} );
                
                atomArr = array[0].getAtoms().toArray(new Atom [array[0].size()]);

                this.getSAAreaTable().put(aa,
                        new double [] {CalcSurfVol.SAsurface(atomArr, ff, 4)});

                this.getVDWAreaTable().put(aa,
                        new double[] {CalcSurfVol.VDWsurface(atomArr, ff, 4)} );
            } catch (StructureException ex) {
                Logger.getLogger(RotamerLibrary.class.getName()).log(Level.SEVERE, null, ex);
            }

            return;
        }

        //calculate the angles and other properties of the rotamers

        int nangles=this.getNumAngles(aa);

        List<Atom> rot_atoms=this.getRotamerAtoms(array[0], nangles);


        double [][] angles=new double[array.length][nangles];
        Atom [] centerMass=new Atom [array.length];
        double [] SAareas=new double [array.length];
        double [] VDWareas=new double [array.length];
        Atom [][] sidechains=null;
        double [] sideVDWareas=new double [array.length];
        try {
            sidechains = new Atom[array.length][getSideChainAtoms((AminoAcid) array[0]).length];
        } catch (StructureException ex) {
            Logger.getLogger(RotamerLibrary.class.getName()).log(Level.SEVERE, null, ex);
        }

       
        for(int i=0;i<array.length;i++){
            for(int j=0;j<nangles;j++){

                try{
                    //calcs the torsion angles of the rotamer
                    angles[i][j]=Calc.torsionAngle( array[i].getAtom( rot_atoms.get(j).getName() ),
                        array[i].getAtom( rot_atoms.get(j+1).getName() ),
                        array[i].getAtom( rot_atoms.get(j+2).getName() ),
                        array[i].getAtom( rot_atoms.get(j+3).getName() ) );

                }
                catch(Exception e){

                }
            }


            try {
                //calcs the center of mass of the rotamer
                sidechains[i]=getSideChainAtoms((AminoAcid) array[i]);
                centerMass[i] = Calc.getCentroid( sidechains[i] );

                if( aa.equals("lys") ){
                    // TODO fix this CGAL bug
                    ///sideVDWareas[i]=CalcSurfVol.VDWsurface( sidechains[i], ff, 4);
                    sideVDWareas[i]=115.0;
                }
                else{
                    sideVDWareas[i]= CalcSurfVol.VDWsurface( sidechains[i], ff, 4);
                }
            } catch (StructureException ex) {
                centerMass[i] =null;
                Logger.getLogger(RotamerLibrary.class.getName()).log(Level.SEVERE, null, ex);
            }

            

            //cals the SA and VdW surface area of the rotamer
            atomArr = array[i].getAtoms().toArray(new Atom [array[i].size()]);
            try {
                SAareas[i] = CalcSurfVol.SAsurface(atomArr, ff, 4);
                VDWareas[i] = CalcSurfVol.VDWsurface(atomArr, ff, 4);
            } catch (StructureException ex) {
                Logger.getLogger(RotamerLibrary.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        //fill the HashMap of angles indexed by AA name
        this.getAnglesTable().put( aa, angles);
        //fill the HashMap of CoM indexed by AA name
        this.getSidechainCoMTable().put(aa, centerMass );
        //fill the HashMap of sidechain's atoms
        this.getSidechainAtomTable().put(aa, sidechains);
        //fill the HashMap of SA and VDW surface area indexed by AA name
        this.getSAAreaTable().put(aa, SAareas);
        this.getVDWAreaTable().put(aa, VDWareas);
        this.getSideVDWAreaTable().put(aa, sideVDWareas);

    }

    /**
     * returns an array with the sidechains chi angles
     *
     * @param aa
     * @return
     */
    private double [] calcSidechainAngles( Group aa ){

        int nangles=this.getNumAngles(aa.getPDBName());

        List<Atom> rot_atoms=this.getRotamerAtoms(aa, nangles);

        double [] angles = new double [nangles];

        if( nangles==0 ){
            return angles;
        }

        for(int j=0;j<nangles;j++){

            try{
                angles[j]=Calc.torsionAngle( aa.getAtom( rot_atoms.get(j).getName() ),
                    aa.getAtom( rot_atoms.get(j+1).getName() ),
                    aa.getAtom( rot_atoms.get(j+2).getName() ),
                    aa.getAtom( rot_atoms.get(j+3).getName() ) );

            }
            catch(Exception e){

            }
        }

        return angles;
    }

    /**
     * @return the AnglesTable
     */
    private HashMap<String,double[][]> getAnglesTable() {
        return AnglesTable;
    }

    /**
     * @param AnglesTable the AnglesTable to set
     */
    private void setAnglesTable(HashMap<String,double[][]> AnglesTable) {
        this.AnglesTable = AnglesTable;
    }

    /**
     * return the sidechain's  number of angles (chi)
     *
     * @param aa
     * @return
     */
    public Integer getNumAngles(String aa){

        return this.getNumAnglesTable().get(internalName(aa));
    }

    /**
     * @return the NumAnglesTable
     */
    private HashMap<String,Integer> getNumAnglesTable() {
        return NumAnglesTable;
    }

    /**
     * @param NumAnglesTable the NumAnglesTable to set
     */
    private void setNumAnglesTable(HashMap<String,Integer> NumAnglesTable) {
        this.NumAnglesTable = NumAnglesTable;
    }

    /**
     * Utility metod, get an array with the sidechain atoms (CB excluded)
     *
     * @param aa
     * @return
     */
     public static Atom[] getSideChainAtoms(AminoAcid aa) throws StructureException {

         ArrayList<Atom> array=new ArrayList<Atom>();

         Iterator<Atom> it=aa.getAtoms().iterator();

         while(it.hasNext()){

             Atom atom=it.next();

             String name=atom.getName().trim();

             if( name.equals("N") ){
                 continue;
             }
             else if( name.equals("C") ){
                 continue;
             }
             else if( name.equals("O") ){
                 continue;
             }
             else if( name.equals("CA") ){
                 continue;
             }
             else if( name.equals("CB") ){
                 continue;
             }
             else if( name.equals("HA") ){
                 continue;
             }
             else if( name.equals("HN") ){
                 continue;
             }
             else if( name.equals("H") ){
                 continue;
             }
             else if( name.equals("H1") ){
                 continue;
             }
             else if( name.equals("H2") ){
                 continue;
             }
             else if( name.equals("H3") ){
                 continue;
             }
             
             array.add(atom);
         }

         if( array.isEmpty() ){
             array.add( Calc.createVirtualCBAtom(aa) );
         }


         Atom [] a=new Atom [array.size()];

        return array.toArray( a );
    }

     /**
      * get the sidechain's atoms of a rotamer
      * @param aa
      * @param idx
      * @return
      */
     public Atom[] getSideChainAtoms(String aa, int idx){
         return this.getSidechainAtomTable().get(internalName(aa))[idx];
     }

     public Atom getSideChainCoM(String aa, int idx){
         return this.getSidechainCoMTable().get(internalName(aa))[idx];
     }

    /**
     * @return the SidechainCoMTable
     */
    private HashMap<String, Atom[]> getSidechainCoMTable() {
        return SidechainCoMTable;
    }

    /**
     * @param SidechainCoMTable the SidechainCoMTable to set
     */
    private void setSidechainCoMTable(HashMap<String, Atom[]> SidechainCoMTable) {
        this.SidechainCoMTable = SidechainCoMTable;
    }

    /**
     * returns true if the atom is a Hidrogen
     * 
     * @param atom
     * @return
     */
    public static boolean isHidrogen(Atom atom) {

        String name=atom.getName();

        if( name.contains("H")){
            if( name.startsWith("N") || name.startsWith("O") ||
                   name.startsWith("C") ){
                return false;
            }
            else{
                return true;
            }
        }

        return false;
    }

    /**
     * @return the AreaTable
     */
    private HashMap<String, double[]> getSAAreaTable() {
        return SAAreaTable;
    }

    /**
     * @param AreaTable the AreaTable to set
     */
    private void setSAAreaTable(HashMap<String, double[]> AreaTable) {
        this.SAAreaTable = AreaTable;
    }

    /**
     * @return the VDWAreaTable
     */
    private HashMap<String, double[]> getVDWAreaTable() {
        return VDWAreaTable;
    }

    /**
     * @param VDWAreaTable the VDWAreaTable to set
     */
    private void setVDWAreaTable(HashMap<String, double[]> VDWAreaTable) {
        this.VDWAreaTable = VDWAreaTable;
    }

    /**
     * convert rotamer library's names to CHARMM
     *
     * @param charmm_name
     * @return
     */
    private static String toCharmmNames( String name ){

        try{
            String number=name.substring(0,1);
            Integer.parseInt( number );
            name=name.substring(1)+number;
        }
        catch( NumberFormatException ex){

        }

        if( name.equals("H") || name.equals("H3") ){
            return "HN";
        }

        return name;
    }

    /**
     * normalize the atom names to CHARMM
     * 
     * @param group
     */
    private static void normalizeNames(Group group) {

        for( Atom at : group.getAtoms() ){
            at.setName( toCharmmNames(at.getName()) );
        }
    }

    /**
     * converts the PDB residue name to charmm
     * @param pdb_name
     * @return
     */
    static public String residue2Charmm(String pdb_name){

        if( pdb_name.equals("HIS") ){
            // HSD neu. prot ND1
            // HSE neu. prot NE2
            // HSP protonated HIS (NE2+)
            return "HSE";
        }

        return pdb_name;
    }

    /**
     * @return the AvgSAAreaTable
     */
    private HashMap<String, Double> getAvgSAAreaTable() {
        return AvgSAAreaTable;
    }

    /**
     * @return the AvgVDWAreaTable
     */
    private HashMap<String, Double> getAvgVDWAreaTable() {
        return AvgVDWAreaTable;
    }
    
    /**
     * calcs avg SA and VdW surf. area
     */
    private void calcAvgSurfAreas() {

        for( String aa : AminoAcids ) {

            //calcs avg SA surf. area
            double [] areas=this.getSAAreaTable().get(aa);

            double avg=0;

            for(int i=0;i<areas.length;i++){
                avg+=areas[i];
            }
            avg/=(double)areas.length;

            this.getAvgSAAreaTable().put(aa, avg);

            //calcs avg VDW surf. Area
            areas=this.getVDWAreaTable().get(aa);

            avg=0;

            for(int i=0;i<areas.length;i++){
                avg+=areas[i];
            }
            avg/=(double)areas.length;

            this.getAvgVDWAreaTable().put(aa, avg);

            //calcs avg sidechain VDW surf. Area
            areas=this.getSideVDWAreaTable().get(aa);

            avg=0;

            for(int i=0;i<areas.length;i++){
                avg+=areas[i];
            }
            avg/=(double)areas.length;

            this.getAvgSideVDWAreaTable().put(aa, avg);
        }
    }

    /**
     * @return the SidechainAtomTable
     */
    private HashMap<String, Atom[][]> getSidechainAtomTable() {
        return SidechainAtomTable;
    }

    /**
     * @param SidechainAtomTable the SidechainAtomTable to set
     */
    private void setSidechainAtomTable(HashMap<String, Atom[][]> SidechainAtomTable) {
        this.SidechainAtomTable = SidechainAtomTable;
    }

    /**
     * @return the AvgSideVDWAreaTable
     */
    private HashMap<String, Double> getAvgSideVDWAreaTable() {
        return AvgSideVDWAreaTable;
    }

    /**
     * @param AvgSideVDWAreaTable the AvgSideVDWAreaTable to set
     */
    private void setAvgSideVDWAreaTable(HashMap<String, Double> AvgSideVDWAreaTable) {
        this.AvgSideVDWAreaTable = AvgSideVDWAreaTable;
    }

    /**
     * @return the SideVDWAreaTable
     */
    private HashMap<String, double[]> getSideVDWAreaTable() {
        return SideVDWAreaTable;
    }

    /**
     * @param SideVDWAreaTable the SideVDWAreaTable to set
     */
    private void setSideVDWAreaTable(HashMap<String, double[]> SideVDWAreaTable) {
        this.SideVDWAreaTable = SideVDWAreaTable;
    }
    
}
