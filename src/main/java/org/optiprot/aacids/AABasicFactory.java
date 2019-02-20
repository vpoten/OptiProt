/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.aacids;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.optiprot.OptiProtParameters;
import org.optiprot.maths.CalcClashes;
import org.optiprot.maths.CalcTransform;
import org.optiprot.rotamer.RotamerLibrary;

/**
 *
 * @author victor
 */
public class AABasicFactory {

    /**
     * create an array of IAABasic for a ProtFANSSolution
     *
     * @param a_protChain
     * @param parameters
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static IAABasic[] create( Chain a_protChain, OptiProtParameters parameters )
            throws  StructureException {

        Iterator<Group> aacids=a_protChain.getAtomGroups().iterator();

        IAABasic [] residues=new IAABasic [a_protChain.getAtomGroups().size()];

        int i=0;

        while( aacids.hasNext() ){

            residues[i++] = new AABasicCaCbNC( aacids.next(), parameters );
        }

        return residues;
    }

    /**
     * converts the given solution to structure chain
     *
     * @param chrom
     * @return
     */
    public static Chain toChain( IAABasic[] solution, boolean generatesH ){

        Chain chain=new ChainImpl();

        List<Group> list = new ArrayList<Group>();

        int contAtoms=1;

        for(int i=0;i<solution.length;i++){

            AminoAcid aa=solution[i].toAminoAcid(generatesH);

            aa.setPDBCode( Integer.toString(i+1) );

            for( Atom atom : aa.getAtoms() ){
                atom.setPDBserial(contAtoms++);
            }

            list.add( aa  );
        }

        chain.setAtomGroups(list);
        return chain;
    }
    
    /**
     * clone a IAABasic []
     * @param solution
     * @return
     */
    public static IAABasic[] clone( IAABasic[] solution ){

        IAABasic[] cloned=new IAABasic[solution.length];

        for(int i=0;i<solution.length;i++){
            cloned[i]=(IAABasic) solution[i].clone();
        }

        return cloned;
    }

    /**
     * converts the given solution to Atom array
     *
     * @param solution
     * @return
     */
     public static Atom [] toAtoms( IAABasic[] solution ){

         Atom [] arr=new Atom [ solution.length*solution[0].getNumAtoms()];

         int cont=0;

         for( IAABasic aabasic : solution ){
            for( Atom atom : aabasic.getAtoms() ){
                arr[cont++]=atom;
            }
         }

         return arr;
     }


     /**
      * creates a random structure given a sequence of amino acids
      *
      * @param sequence
      * @param parameters
      * @return
      * @throws org.biojava.bio.structure.StructureException
      */
      public static IAABasic[] create( String sequence, OptiProtParameters parameters )
              throws StructureException {

          IAABasic [] residues=new IAABasic [sequence.length()];

          RotamerLibrary rotLib=parameters.getRotLib();

          Group aacid=rotLib.getRotamer( String.valueOf(sequence.charAt(0)) );
          residues[0]=new AABasicCaCbNC( aacid, parameters );

          for( int i=1; i<sequence.length(); i++){
              
              aacid=rotLib.getRotamer( String.valueOf(sequence.charAt(i)) );
              IAABasic newAA = new AABasicCaCbNC( aacid, parameters );

              CalcTransform.addAminoAcid(newAA, residues, i-1);

              residues[i]=newAA;
          }

          CalcClashes.fixClashes( 0.6, residues, rotLib, 100);
          
          //center the chain at origin
          CalcTransform.center(residues);
          
          return residues;
      }
      
}
