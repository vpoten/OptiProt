/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.align.ClusterAltAligs;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.jama.Matrix;

/**
 *
 * @author victor
 */
public class CalcAlignment {

    /**
     * calculates the alignment and transform chain2 to meet chain1
     *
     * @param chain1
     * @param chain2
     * @throws org.biojava.bio.structure.StructureException
     */
    public static void align( Chain chain1, Chain chain2 )
            throws StructureException {

        StructurePairAligner sc = new StructurePairAligner();

        Structure s1=new StructureImpl();
        s1.addChain(chain1);

        Structure s2=new StructureImpl();
        s2.addChain(chain2);

        // do the calculations
        sc.align(s1, s2);

        AlternativeAlignment[] aligs = sc.getAlignments();

        // cluster similar results together
        ClusterAltAligs.cluster(aligs);

        // transform chain 2 to meet chain1
        if (aligs.length > 0){

            Matrix mtrans=Matrix.identity(4, 4);
            Matrix mrot=Matrix.identity(4, 4);

            mrot.setMatrix(0, 2, 0, 2, aligs[0].getRotationMatrix().transpose() );

            mtrans.set(0, 3, aligs[0].getShift().getX() );
            mtrans.set(1, 3, aligs[0].getShift().getY() );
            mtrans.set(2, 3, aligs[0].getShift().getZ() );

            mtrans = mtrans.times(mrot);

           
            CalcTransform.applyTransform( 
                    chain2, mtrans, 0, chain2.getAtomLength() );

        }


    }

}
