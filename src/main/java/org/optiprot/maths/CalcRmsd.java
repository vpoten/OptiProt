/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.Iterator;
import java.util.List;
import org.biojava.bio.structure.*;
import org.biojava.bio.structure.jama.EigenvalueDecomposition;
import org.biojava.bio.structure.jama.Matrix;

/**
 *
 * @author victor
 */
public class CalcRmsd {

    /**
     * calcs the RMSD between the given chains using C,N,CA,CB,O
     *
     * @param chain1
     * @param chain2
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static double rmsd(Chain chain1, Chain chain2)
            throws StructureException{

        if( !chain1.getSeqResSequence().equals(chain2.getSeqResSequence()) ){
            throw new StructureException("Chains not equals.");
        }

       
        Iterator<Group> it1=chain1.getAtomGroups().iterator();
        Iterator<Group> it2=chain2.getAtomGroups().iterator();

        Atom [] atoms1=new Atom [chain1.getAtomGroups().size()*5];
        Atom [] atoms2=new Atom [chain1.getAtomGroups().size()*5];

        int i=0;

        while(it1.hasNext()){

            AminoAcid aa1=(AminoAcid) it1.next();
            AminoAcid aa2=(AminoAcid) it2.next();

            atoms1[i]=  (Atom) aa1.getCA().clone();
            atoms2[i]=  (Atom) aa2.getCA().clone();
            atoms1[i+1]=  (Atom) aa1.getC().clone();
            atoms2[i+1]=  (Atom) aa2.getC().clone();
            
            try{
                atoms1[i+2]=  (Atom) aa1.getCB().clone();
                atoms2[i+2]=  (Atom) aa2.getCB().clone();
            }
            catch( StructureException ex){
                atoms1[i+2]= Calc.createVirtualCBAtom(aa1);
                atoms2[i+2]= Calc.createVirtualCBAtom(aa2);
            }

            atoms1[i+3]=  (Atom) aa1.getN().clone();
            atoms2[i+3]=  (Atom) aa2.getN().clone();
            atoms1[i+4]=  (Atom) aa1.getO().clone();
            atoms2[i+4]=  (Atom) aa2.getO().clone();

            i+=5;
        }

        minimizeRmsd(atoms1,atoms2,5);

        return rmsd(atoms1,atoms2);

    }

    /**
     * calcs the RMSD between the given chains using CA
     *
     * @param chain1
     * @param chain2
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static double rmsdCA(Chain chain1, Chain chain2)
            throws StructureException{

        if( !chain1.getSeqResSequence().equals(chain2.getSeqResSequence()) ){
            throw new StructureException("Chains not equals.");
        }

        Iterator<Group> it1=chain1.getAtomGroups().iterator();
        Iterator<Group> it2=chain2.getAtomGroups().iterator();

        Atom [] atoms1=new Atom [chain1.getAtomGroups().size()];
        Atom [] atoms2=new Atom [chain1.getAtomGroups().size()];

        int i=0;

        while(it1.hasNext()){

            AminoAcid aa1=(AminoAcid) it1.next();
            AminoAcid aa2=(AminoAcid) it2.next();

            atoms1[i]=  (Atom) aa1.getCA().clone();
            atoms2[i]=  (Atom) aa2.getCA().clone();
            
            i++;
        }

        minimizeRmsd(atoms1,atoms2,1);

        return rmsd(atoms1,atoms2);
    }


    /**
     * Translates and rotates chain1 to get the min rmsd between the chains.
     * chain1 is modified
     *
     * @param chain1
     * @param chain2
     * @param step  number of atoms per aminoacid (the group's first one is CA)
     */
    private static void minimizeRmsd(Atom [] chain1, Atom [] chain2, int step)
            throws StructureException{

        Atom center1=CalcGeom.getCentroid(chain1, step);
        Atom center2=CalcGeom.getCentroid(chain2, step);

        // translate chain1 & chain2 to the origin
        Atom shiftVect1=new AtomImpl();
        shiftVect1.setCoords( new double [] {-center1.getX(),-center1.getY(),-center1.getZ()});
        Atom shiftVect2=new AtomImpl();
        shiftVect2.setCoords( new double [] {-center2.getX(),-center2.getY(),-center2.getZ()});

        for(int i=0;i<chain1.length;i++){
            Calc.shift(chain1[i], shiftVect1);
            Calc.shift(chain2[i], shiftVect2);
        }

        // Do the alignment using quaternions, paper link;
        // http://cnx.org/content/m11608/latest/

        Matrix A=new Matrix(4,4);
        Matrix B=new Matrix(4,4);
        Matrix N=new Matrix(4,4);

        for(int i=0;i<chain1.length;i+=step){

            A.set(1, 0, -chain1[i].getX() );
            A.set(2, 0, -chain1[i].getY() );
            A.set(3, 0, -chain1[i].getZ() );
            A.set(0, 1, chain1[i].getX() );
            A.set(2, 1, chain1[i].getZ() );
            A.set(3, 1, -chain1[i].getY() );
            A.set(0, 2, chain1[i].getY() );
            A.set(1, 2, -chain1[i].getZ() );
            A.set(3, 2, chain1[i].getX() );
            A.set(0, 3, chain1[i].getZ() );
            A.set(1, 3, chain1[i].getY() );
            A.set(2, 3, -chain1[i].getX() );

            B.set(0, 1, -chain2[i].getX() );
            B.set(0, 2, -chain2[i].getY() );
            B.set(0, 3, -chain2[i].getZ() );
            B.set(1, 0, chain2[i].getX() );
            B.set(1, 2, -chain2[i].getZ() );
            B.set(1, 3, chain2[i].getY() );
            B.set(2, 0, chain2[i].getY() );
            B.set(2, 1, chain2[i].getZ() );
            B.set(2, 3, -chain2[i].getX() );
            B.set(3, 0, chain2[i].getZ() );
            B.set(3, 1, -chain2[i].getY() );
            B.set(3, 2, chain2[i].getX() );

            N.plusEquals( A.times(B) );
        }

        EigenvalueDecomposition eigen=N.eig();

        double [] eigValues=eigen.getRealEigenvalues();

        // get the most positive eigenvalue
        int max=0;
        double maxval=eigValues[0];

        for(int i=1;i<4;i++){

            if( eigValues[i]>maxval ){
                maxval=eigValues[i];
                max=i;
            }
        }

        Matrix eigVect=eigen.getV();

        // the eigenvector associated to the most eigenvalue is the quaternion
        Quaternion quat=new Quaternion( eigVect.get(0,max), eigVect.get(1,max),
                eigVect.get(2,max), eigVect.get(3,max) );

        quat=quat.normalize();

        for(int i=0;i<chain1.length;i++){

            chain1[i]=quat.rotateVector(chain1[i]);
        }
        
    }

    /**
     * 
     * @param chain1
     * @param chain2
     * @return
     */
    private static double rmsd(Atom [] chain1, Atom [] chain2) {

        double rmsd=0;

        for(int i=0;i<chain1.length;i++){

            rmsd+=CalcGeom.squareDistance(chain1[i],chain2[i]);
        }

        return Math.sqrt(rmsd/chain1.length);
    }


    /**
     * calculates the rmsd without preprocessing
     * @pre : chain1 and chain2 have the same number of atoms.
     *
     * @param chain1
     * @param chain2
     * @return
     */
    public static double bruteRmsd( List<Atom> chain1, List<Atom> chain2) {

        double rmsd=0;

        for(int i=0;i<chain1.size();i++){

            rmsd+=CalcGeom.squareDistance(chain1.get(i),chain2.get(i));
        }

        return Math.sqrt(rmsd/chain1.size());
    }


    /**
     * calculates the rmsd without preprocessing
     * 
     * @param chain1
     * @param chain2
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static double bruteRmsd(Chain chain1, Chain chain2)
            throws StructureException {

        double rmsd=0;

        List<Group> l1=chain1.getAtomGroups();
        List<Group> l2=chain2.getAtomGroups();

        for(int i=0;i<l1.size();i++){

            rmsd+=CalcGeom.squareDistance( l1.get(i).getAtom("CA"),
                    l2.get(i).getAtom("CA") );
        }

        return Math.sqrt(rmsd/l1.size());
    }


}
