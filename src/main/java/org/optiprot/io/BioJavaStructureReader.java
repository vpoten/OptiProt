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
import java.util.ArrayList;
import org.biojava.bio.BioException;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 *
 * @author victor
 */
public class BioJavaStructureReader {

    public BioJavaStructureReader()
    {
        
    }

    public static ArrayList<Structure> readStructuresFromDir( String path )
    {
        ArrayList<Structure> array=null;

        File dir=new File(path);

        File[] list_files = dir.listFiles();

        if( list_files==null )
            return null;

        array=new ArrayList<Structure>();

        try{

            for( int i=0; i<list_files.length; i++)
            {
                if( list_files[i].isDirectory() )
                    continue;

                array.add( readStructure( list_files[i].getCanonicalPath() ) );
            }
        }
        catch( Exception e )
        {
            return null;
        }

        return array;
    }

    public static Structure readStructure( String path, String id )
    {

         PDBFileReader pdbreader = new PDBFileReader();

         pdbreader.setPath(path);

         BioJavaStructureReader.setParametersDefault( pdbreader );

         try{

             return pdbreader.getStructureById( id );

         } catch (Exception e) {

             return null;
         }
    }

    public static Structure readStructure( String path_file) throws IOException
    {

         PDBFileReader pdbreader = new PDBFileReader();

         BioJavaStructureReader.setParametersDefault( pdbreader );

         return pdbreader.getStructure( path_file );

    }

    private static  void setParametersDefault( PDBFileReader pdbreader )
    {
        // the following parameters are optional:

         //the parser can read the secondary structure
         // assignment from the PDB file header and add it to the amino acids
         pdbreader.setParseSecStruc(true);

         // align the SEQRES and ATOM records, default = true
         // slows the parsing speed slightly down, so if speed matters turn it off.
         pdbreader.setAlignSeqRes(true);

         // parse the C-alpha atoms only, default = false
         pdbreader.setParseCAOnly(false);

         // download missing PDB files automatically from EBI ftp server, default = false
         pdbreader.setAutoFetch(false);
    }

    /**
     * read a sequence from fasta file
     *
     * @param filename
     * @return
     */
    public static String readSequence(String filename)
            throws FileNotFoundException, BioException {

        String seq="";

        //setup file input
        BufferedReader in = new BufferedReader(new FileReader(filename));

        RichSequenceIterator it=RichSequence.IOTools.readFastaProtein(in, null);

        //gets the first sequence
        if( it.hasNext() ){
            seq=it.nextRichSequence().seqString();
        }

        return seq;
    }

    /**
     * 
     * @param pdbdir
     * @param pdbname
     * @return
     */
    public static Structure readStructurePDB( String pdbdir, String pdbname )
            throws IOException{

        return readStructure( pdbdir+File.separator+
                pdbname.substring(1,3).toLowerCase()+
                File.separator+"pdb"+pdbname.toLowerCase()+".ent.gz" );
    }

   

    /**
     * get the chains (only amino acids) whose indices are in numchains (joined)
     *
     * @param struct
     * @param numchains
     * @return
     */
    public static Chain getChains( Structure struct, int [] numchains ){
        return getChains( struct, numchains, true );
    }

    /**
     * get the chains whose indices are in numchains (joined)
     * @param struct
     * @param numchains
     * @param onlyAA : if true gets only AminoAcids
     * @return
     */
    public static Chain getChains( Structure struct, int [] numchains, boolean onlyAA ){

        Chain ch2=new ChainImpl();

        for( int i=0;i<numchains.length;i++){
            Chain ch=struct.getChain(numchains[i]);

            for( Group group : ch.getAtomGroups() ){
                if( onlyAA ){
                    if( group instanceof AminoAcid ){
                        ch2.addGroup(group);
                    }
                }
                else{
                    ch2.addGroup(group);
                }
            }
        }

        return ch2;
    }

    /**
     * get all chains in structure (joined)
     *
     * @param struct
     * @param remove : chains to remove ( X,Y,... )
     * @param onlyAA : if true gets only AminoAcids
     * @return
     */
    public static Chain getAllChains( Structure struct, String remove, boolean onlyAA ){

        Chain ch2=new ChainImpl();

        for( Chain chain : struct.getChains() ){

           if( isChainToRemove( chain ,remove ) )
                continue;

            for( Group group : chain.getAtomGroups() ){
                if( onlyAA ){
                    if( group instanceof AminoAcid ){
                        ch2.addGroup(group);
                    }
                }
                else{
                    ch2.addGroup(group);
                }
            }
        }

        return ch2;
    }


    /**
     * 
     * @param ch
     * @return
     */
    public static Chain getAminoAcids( Chain ch ){

        Chain ch2=new ChainImpl();
        ch2.setName( ch.getName() );

        for( Group group : ch.getAtomGroups() ){

            if( group instanceof AminoAcid ){
                ch2.addGroup(group);
            }
        }

        return ch2;
    }


    /**
     * get all chains in structure (joined)
     *
     * @param struct
     * @param remove : chains to remove ( X,Y,... )
     * @return
     */
    public static Chain getAllChains( Structure struct, String remove ){

        Chain ch2=new ChainImpl();

        for( Chain chain : struct.getChains() ){

           if( isChainToRemove( chain ,remove ) )
                continue;

            for( Group group : chain.getAtomGroups() ){

                if( group instanceof AminoAcid ){
                    ch2.addGroup(group);
                }
            }
        }

        return ch2;
    }


    /**
     * returns the heteroatom with given name name (if exists)
     *
     * @param ch
     * @param hetname : name of the required heteroatom
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static HetatomImpl getHetatom( Chain ch, String hetname)
            throws StructureException{

        for( Group group : ch.getAtomGroups() ){

            if( group instanceof HetatomImpl ){
                if( group.getPDBName().equals(hetname) )
                    return (HetatomImpl) group;
            }
        }

        throw new StructureException("Hetatom not found");
    }

    /**
     * returns the chain index wich contains the heteroatom
     *
     * @param struct
     * @param hetname
     * @return
     */
    public static int findHetatom( Structure struct, String hetname )
            throws StructureException{

        int idx=0;

        for( Chain chain : struct.getChains() ){

            for( Group group : chain.getAtomGroups() ){
                if( group instanceof HetatomImpl ){
                    if( group.getPDBName().equals(hetname) )
                        return idx;
                }
            }

            idx++;
        }

        throw new StructureException("Hetatom not found");
    }

    /**
     * 
     * @param struct
     * @param hetname
     * @param remove : chains to remove ( X,Y,... )
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static HetatomImpl getHetatom( Structure struct, String hetname,
            String remove ) throws StructureException{

        for( Chain chain : struct.getChains() ){

            if( isChainToRemove( chain ,remove ) )
                continue;

            for( Group group : chain.getAtomGroups() ){
                if( group instanceof HetatomImpl ){
                    if( group.getPDBName().equals(hetname) )
                        return (HetatomImpl) group;
                }
            }
        }

        throw new StructureException("Hetatom not found");
    }

    /**
     * returns the chain that contains the heteroatom
     *
     * @param struct
     * @param hetname
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static Chain findHetatomChain( Structure struct, String hetname )
            throws StructureException{

        int index=findHetatom(struct, hetname);
        return struct.getChain( index );
    }
    
    /**
     * returns true if the chain's name is in the remove string
     * 
     * @param chain
     * @param remove
     * @return
     */
    private static boolean isChainToRemove( Chain chain , String remove ){

        if( remove!=null ){
            for( int i=0;i<remove.length();i++){
                if( remove.charAt(i)==chain.getName().charAt(0) ){
                     return true;
                }
            }
        }

        return false;
    }

}
