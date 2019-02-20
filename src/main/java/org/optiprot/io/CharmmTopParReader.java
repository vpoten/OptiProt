/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.optiprot.potential.IForceFieldTable;
import org.optiprot.potential.element.CharmmElement;
import org.optiprot.potential.element.CharmmResidue;
import org.optiprot.potential.element.MolecularBond;
import org.optiprot.potential.element.MolecularGroup;
import org.optiprot.potential.element.MolecularImproper;

/**
 * Class for reading of CHARMM parameter and topology files
 * @author victor
 */
public class CharmmTopParReader {

    //topology file fields
    final static private String MASS="MASS";
    final static private String RESI="RESI";
    final static private String GROUP="GROUP";
    final static private String ATOM="ATOM";
    final static private String BOND="BOND";
    final static private String DOUBLE="DOUBLE";
    final static private String TRIPLE="TRIPLE";
    final static private String AROMATIC="AROMATIC";
    final static private String IMPR="IMPR";
    final static private String ANGLE="ANGLE";
    final static private String CMAP="CMAP";
    final static private String DONOR="DONOR";
    final static private String ACCEPTOR="ACCEPTOR";
    final static private String IC="IC";
    final static private String PRES="PRES";
    final static private String COMMENT="!";
    final static private String HEADER="*";
    final static private String END="END";

    //parameters file fields
    final static private String BONDS="BONDS";
    final static private String ANGLES="ANGLES";
    final static private String DIHEDRALS="DIHEDRALS";
    final static private String IMPROPER="IMPROPER";
    final static private String NONBONDED="NONBONDED";
    final static private String HBOND="HBOND";

    ////////////////////////////////////////////

    /**
     *
     * @param path
     * @param tableParams : I/O
     * @param elemTable : charmm elements (readed in topology file)
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     */
    static public void parseParFile( String path, IForceFieldTable tableParams,
            ArrayList<CharmmElement> elemList )
            throws FileNotFoundException, IOException, URISyntaxException {

        BufferedReader r = null;

        if( path.startsWith("/") ){
            //load from filesystem
            r = new BufferedReader( new FileReader(path) );
        }
        else{
            //load as resource
            URL url=CharmmTopParReader.class.getClassLoader().getResource(path);
            r = new BufferedReader( new FileReader(new File(url.toURI())) );
        }

        //fill the table indexed by elem type
        HashMap<String,CharmmElement> elemTable = new HashMap<String,CharmmElement>();

        for( CharmmElement ele : elemList ){
            if(ele!=null)
                elemTable.put( ele.getType(), ele);
        }

        String l_line=r.readLine();

        String [] tokens = new String [16];
        int ntokens=0;
        String section="";

        while( l_line!=null )
        {
           l_line=l_line.trim();

           if( l_line.isEmpty() ){
                //nothing to do
           }
           else if( l_line.startsWith(COMMENT) ){
               // a comment
           }
           else if( l_line.startsWith(HEADER) ){
               // a header
           }
           else if( l_line.startsWith(BONDS) ){
               // bonds section
               section=BONDS;
           }
           else if( l_line.startsWith(ANGLES) ){
               // angles section
               section=ANGLES;
           }
           else if( l_line.startsWith(DIHEDRALS) ){
               // dihedrals section
               section=DIHEDRALS;
           }
           else if( l_line.startsWith(IMPROPER) ){
               // improper section
               section=IMPROPER;
           }
           else if( l_line.startsWith(NONBONDED) ){
               // nonbonded section
               section=NONBONDED;
           }
           else if( l_line.startsWith(CMAP) ){
               // cmap section
               section=CMAP;
           }
           else if( l_line.startsWith(HBOND) ){
               // hbond section
               section=HBOND;
           }
           else if( l_line.startsWith(END) ){
               // end
               section=END;
           }
           else{
               // parse line
               ntokens=getTokens( l_line, " \t", 12, tokens);

               if( section.equals(BONDS) ){
                   // at1 at2 kb b0
                    tableParams.addBond( elemTable.get(tokens[0]).getIndex(),
                            elemTable.get(tokens[1]).getIndex(),
                            Double.parseDouble(tokens[2]), Double.parseDouble(tokens[3]));
               }
               else if( section.equals(ANGLES) ){
                   // at1 at2 at3 ktheta theta0 kub s0
                   double kub=0;
                   double s0=0;

                   if(ntokens>=7){
                       try{ kub=Double.parseDouble(tokens[5]); }
                       catch(Exception e){ kub=0; }

                       try{ s0=Double.parseDouble(tokens[6]); }
                       catch(Exception e){ s0=0; }
                   }

                   tableParams.addAngle( elemTable.get(tokens[0]).getIndex(),
                           elemTable.get(tokens[1]).getIndex(),
                           elemTable.get(tokens[2]).getIndex(),
                           Double.parseDouble(tokens[3]), Double.parseDouble(tokens[4]),
                           kub, s0);

               }
               else if( section.equals(DIHEDRALS) ){
                   // at1 at2 at3 at4 kchi n delta
                   
                   tableParams.addDihedral( 
                           (elemTable.get(tokens[0])==null) ? CharmmElement.WILDCARD : elemTable.get(tokens[0]).getIndex(),
                           elemTable.get(tokens[1]).getIndex(),
                           elemTable.get(tokens[2]).getIndex(),
                           (elemTable.get(tokens[3])==null) ? CharmmElement.WILDCARD : elemTable.get(tokens[3]).getIndex(),
                           Double.parseDouble(tokens[4]), Double.parseDouble(tokens[5]),
                           Double.parseDouble(tokens[6]));
                   
               }
               else if( section.equals(IMPROPER) ){
                   // at1 at2 at3 at4 kpsi (val) psi0

                   tableParams.addImproper( 
                           elemTable.get(tokens[0]).getIndex(),
                           (elemTable.get(tokens[1])==null) ? CharmmElement.WILDCARD : elemTable.get(tokens[1]).getIndex(),
                           (elemTable.get(tokens[2])==null) ? CharmmElement.WILDCARD : elemTable.get(tokens[2]).getIndex(),
                           elemTable.get(tokens[3]).getIndex(),
                           Double.parseDouble(tokens[4]), Double.parseDouble(tokens[6]));

               }
               else if( section.equals(NONBONDED) ){
                   // at (val) epsilon rmin/2 (val) eps14 rmin/2_14
                   double eps14=0;
                   double rmin2_14=0;

                   try{
                       Double.parseDouble(tokens[1]);
                       Double.parseDouble(tokens[2]);
                       Double.parseDouble(tokens[3]);
                   }
                   catch(Exception e){
                       l_line=r.readLine();
                       continue;
                   }

                   if(ntokens>=7){
                       try{ eps14=Double.parseDouble(tokens[5]); }
                       catch(Exception e){ eps14=0; }
                       
                       try{ rmin2_14=Double.parseDouble(tokens[6]); }
                       catch(Exception e){ rmin2_14=0; }
                   }

                   tableParams.addNonbonded(
                           elemTable.get(tokens[0]).getIndex(),
                           Double.parseDouble(tokens[2]), Double.parseDouble(tokens[3]),
                           eps14, rmin2_14);

               }
           }


           l_line=r.readLine();
        }

        r.close();
    }

    /**
     *
     * @param path
     * @param elemList : I/O
     * @return
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     */
    static public ArrayList<CharmmResidue> parseTopFile( String path,
             ArrayList<CharmmElement> elemList )
            throws FileNotFoundException, IOException, URISyntaxException {

        ArrayList<CharmmResidue> residues = new ArrayList<CharmmResidue>();

        int lastIndex=-1;//index of the last charmm element
        HashMap<String,CharmmElement> elemTable = new HashMap<String,CharmmElement>();
        HashMap<Integer,CharmmElement> elemTable2 = new HashMap<Integer,CharmmElement>();

        BufferedReader r = null;

        if( path.startsWith("/") ){
            //load from filesystem
            r = new BufferedReader( new FileReader(path) );
        }
        else{
            //load as resource
            URL url=CharmmTopParReader.class.getClassLoader().getResource(path);
            r = new BufferedReader( new FileReader(new File(url.toURI())) );
        }

        String l_line=r.readLine();

        CharmmResidue residue = null;
        MolecularGroup group = null;
        Atom atom=null;
        String [] tokens = new String [16];
        int ntokens=0;

        while( l_line!=null )
        {
           l_line=l_line.trim();

           if( l_line.isEmpty() ){
                //nothing to do
           }
           else if( l_line.startsWith(COMMENT) ){
               // a comment
           }
           else if( l_line.startsWith(HEADER) ){
               // a header
           }
           else if( l_line.startsWith(RESI) ){
               // residue start (name , charge)

               if( residue!=null ){
                   if( group!=null ){
                        residue.addGroup(group);
                        group=null;
                   }
                   residues.add(residue);
               }

               residue=new CharmmResidue();
               //getTokens( l_line, " \t", 3, tokens);
               splitString( l_line, "\\s", 3, tokens);

               residue.setName(tokens[1]);
               residue.setCharge( Double.parseDouble(tokens[2]) );

           }
           else if( l_line.startsWith(MASS) ){
               // atom mass (index, type, mass, element ! description )

               getTokens( l_line, " \t", 5, tokens);

               CharmmElement ele = new CharmmElement();
               ele.setIndex(Integer.parseInt(tokens[1]));
               ele.setType(tokens[2]);
               ele.setMass(Double.parseDouble(tokens[3]));
               ele.setElement(tokens[4]);

              
               ele.setDescription( l_line.substring(l_line.indexOf(COMMENT)+1) );

               //index the elements by type
               if( ele.getIndex() > lastIndex )
                   lastIndex=ele.getIndex();

               elemTable.put(ele.getType(),ele);
               elemTable2.put(ele.getIndex(),ele);
           }
           else if( l_line.startsWith(GROUP) ){
               // group of atoms

               if( residue==null ){
                   l_line=r.readLine();
                   continue;
               }

               if( group!=null )
                   residue.addGroup(group);

               group=new MolecularGroup();
           }
           else if( l_line.startsWith(ATOM) ){
               // atom (name , type , charge)

               if( residue==null ){
                   l_line=r.readLine();
                   continue;
               }

               //getTokens( l_line, " \t", 4, tokens);
               splitString( l_line, "\\s", 4, tokens);
               atom=new AtomImpl();
               atom.setName( tokens[1] );
               group.addAtom( atom, Double.parseDouble(tokens[3]),
                       elemTable.get(tokens[2]).getIndex() );
           }
           else if( l_line.startsWith(BOND) ){
               // simple bond

               if( residue==null ){
                   l_line=r.readLine();
                   continue;
               }

               if( group!=null ){
                   residue.addGroup(group);
                   group=null;
               }

               addBonds( residue, l_line, MolecularBond.BondType.SIMPLE );
           }
           else if( l_line.startsWith(DOUBLE) ){
               // double bond

               if( residue==null ){
                   l_line=r.readLine();
                   continue;
               }

               addBonds( residue, l_line, MolecularBond.BondType.DOUBLE );
           }
           else if( l_line.startsWith(TRIPLE) ){
               // triple bond

               if( residue==null ){
                   l_line=r.readLine();
                   continue;
               }

               addBonds( residue, l_line, MolecularBond.BondType.TRIPLE );
           }
           else if( l_line.startsWith(AROMATIC) ){
               // aromatic bond

               if( residue==null ){
                   l_line=r.readLine();
                   continue;
               }

               addBonds( residue, l_line, MolecularBond.BondType.AROMATIC );
           }
           else if( l_line.startsWith(IMPR) ){
               // improper

               if( residue==null ){
                   l_line=r.readLine();
                   continue;
               }

               addImpropers( residue, l_line);
           }
           else if( l_line.startsWith(ANGLE) ){
               // angle
           }
           else if( l_line.startsWith(DONOR) ){
               // donor
           }
           else if( l_line.startsWith(ACCEPTOR) ){
               // acceptor
           }
           else if( l_line.startsWith(CMAP) ){
               // cmap correction
           }
           else if( l_line.startsWith(IC) ){
               // internal coordinates
           }
           else if( l_line.startsWith(PRES) ){
               // patch residues

               if( residue!=null ){
                   if( group!=null ){
                        residue.addGroup(group);
                        group=null;
                   }
                   residues.add(residue);
                   residue=null;
               }
           }
           else if( l_line.startsWith(END) ){
               // end
               
               if( residue!=null ){
                   if( group!=null ){
                        residue.addGroup(group);
                        group=null;
                   }
                   residues.add(residue);
                   residue=null;
               }
           }
           else{
               // nothing to do
           }

       
           l_line=r.readLine();
        }

        r.close();

        //prepare elemList (insert null element in missing index)
        for( int i=0; i<=lastIndex; i++){
            if( elemTable2.containsKey(i) ){
                elemList.add( elemTable2.get(i) );
            }
            else{
                elemList.add(null);
            }
        }

        return residues;
    }


    /**
     * add molecular impropers to the residue
     * @param residue
     * @param line
     */
    private static void addImpropers(CharmmResidue residue, String line) {

        String [] tokens=new String [16];

        if( line.indexOf( COMMENT )>0 )
            line=line.substring(0,line.indexOf(COMMENT));

        int ntokens=getTokens( line, " \t", 16, tokens);

        for( int i=1; i<ntokens; i+=4){

            Atom at1=new AtomImpl();
            at1.setName( tokens[i]);
            Atom at2=new AtomImpl();
            at2.setName( tokens[i+1]);
            Atom at3=new AtomImpl();
            at3.setName( tokens[i+2]);
            Atom at4=new AtomImpl();
            at4.setName( tokens[i+3]);

            MolecularImproper mi=new MolecularImproper();
            mi.setAtom(0, at1, residue.getAtomCharge(at1.getName()), residue.getAtomType(at1.getName()) );
            mi.setAtom(1, at2, residue.getAtomCharge(at2.getName()), residue.getAtomType(at2.getName()) );
            mi.setAtom(2, at3, residue.getAtomCharge(at3.getName()), residue.getAtomType(at3.getName()) );
            mi.setAtom(3, at4, residue.getAtomCharge(at4.getName()), residue.getAtomType(at4.getName()) );
            residue.getImpropers().add(mi);
        }
    }


    /**
     *
     * @param line : line to parse
     * @param delimit : delimiters
     * @param num : number of tokens to obtain
     * @param tokens : array of tokens I/O
     * @return number of readed tokens
     */
    static private int getTokens( String line, String delimit, int num, String [] tokens ){

        StringTokenizer l_tokenizer = new StringTokenizer(line, delimit);

        int i=0;

        while (l_tokenizer.hasMoreTokens())
        {
            tokens[i++] = l_tokenizer.nextToken();

            if( i==(num-1) )
                break;
        }

        return i;
    }

    static private int splitString( String line, String delimit, int num, String [] tokens ){

        String [] arr=line.split(delimit);

        int j=0;
        for( int i=0;i<arr.length;i++){

            if( !arr[i].isEmpty() ){
                tokens[j++]=arr[i];
            }
            
            if( j==num )
                break;
        }

        return j;
    }

    /**
     * add molecular bonds to the residue
     * @param residue
     * @param line
     * @param type
     */
    static private void addBonds( CharmmResidue residue, String line, MolecularBond.BondType type ){

        String [] tokens=new String [16];

        if( line.indexOf( COMMENT )>0 )
            line=line.substring(0,line.indexOf(COMMENT));


        int ntokens=getTokens( line, " \t", 16, tokens);

        for( int i=1; i<ntokens; i+=2){

            Atom at1=new AtomImpl();
            at1.setName( tokens[i]);
            Atom at2=new AtomImpl();
            at2.setName( tokens[i+1]);

            MolecularBond mb=new MolecularBond(at1,at2);
            mb.setBondType(type);
            mb.setAtomTypeA( residue.getAtomType(at1.getName()) );
            mb.setAtomTypeB( residue.getAtomType(at2.getName()) );
            mb.setChargeA( residue.getAtomCharge(at1.getName()) );
            mb.setChargeB( residue.getAtomCharge(at2.getName()) );
            residue.getBonds().add(mb);
        }
    }
    
}
