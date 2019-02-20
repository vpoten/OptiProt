/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.*;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcGeom;
import org.optiprot.rotamer.RotamerLibrary;

/**
 * Neighbor Vector SASA algorithm calculation
 *
 * @author victor
 */
class NVSASACalc {

    private static double LBOUND=3.3;//in angmstroms
    private static double UBOUND=11.1;
    private static double KNWEIGHT=Math.PI/(UBOUND-LBOUND);


    /**
     *
     * @param molElements
     * @param rotLib
     * @return
     */
    public static double calc(Chain chain, RotamerLibrary rotLib) {

        double sasa=0;

        List<Group> groups = chain.getAtomGroups();

        double ncount=0;
        Atom nvector=new AtomImpl();
        Atom CBi=new AtomImpl();
        Atom vector =new AtomImpl();
        double nweight=0;
        
        for( int i=0;i<groups.size();i++ ){

            ncount=0;
            nvector.setX(0);
            nvector.setY(0);
            nvector.setZ(0);
            String cbname="CB";

            if( !groups.get(i).hasAtom("CB"))
                cbname="CA";

            try {
                Atom tmp = groups.get(i).getAtom(cbname);
                CBi.setX( tmp.getX()*-1.0 );
                CBi.setY( tmp.getY()*-1.0 );
                CBi.setZ( tmp.getZ()*-1.0 );
            } catch (StructureException ex) {
                Logger.getLogger(NVSASACalc.class.getName()).log(Level.SEVERE, null, ex);
                continue;
            }

            //calculate neighbor weights, neighbor counts and neighbor vector
            for( int j=0;j<groups.size();j++ ){

                if(i==j)
                    continue;

                cbname="CB";
                if( !groups.get(j).hasAtom("CB"))
                    cbname="CA";

                try {
                    Atom CBj=groups.get(j).getAtom(cbname);
                    vector.setX( CBj.getX() );
                    vector.setY( CBj.getY() );
                    vector.setZ( CBj.getZ() );
                    CalcGeom.addEquals(vector, CBi);
                    double distance=Calc.amount(vector);

                    if( distance>UBOUND )
                        continue;

                    nweight=neighborWeight2( distance );
                    ncount+=nweight;
                    CalcGeom.product2( vector, nweight/distance );
                    CalcGeom.addEquals( nvector, vector);
                } 
                catch (StructureException ex) {
                    Logger.getLogger(NVSASACalc.class.getName()).log(Level.SEVERE, null, ex);
                }
            }//

            // finish the neighbor vector
            CalcGeom.product2( nvector, 1.0/ncount );

            sasa+=Calc.amount(nvector)*rotLib.getAvgSASA(groups.get(i));
        }

        
        return sasa;
    }

    /**
     * 
     * @param chain
     * @param btree
     * @param rotLib
     * @return
     */
    public static double calc(Chain chain, BSPTree btree, RotamerLibrary rotLib) {

        double sasa=0;

        List<Group> groups = chain.getAtomGroups();

        double ncount=0;
        Atom nvector=new AtomImpl();
        Atom CBi=new AtomImpl();
        Atom vector =new AtomImpl();
        double nweight=0;
        ArrayList<Atom> listNeighb=new ArrayList<Atom>();

        for( int i=0;i<groups.size();i++ ){

            ncount=0;
            nvector.setX(0);
            nvector.setY(0);
            nvector.setZ(0);
            String cbname="CB";
            Atom tmp = null;

            if( !groups.get(i).hasAtom("CB"))
                cbname="CA";

            try {
                tmp = groups.get(i).getAtom(cbname);
            } catch (StructureException ex) {
                Logger.getLogger(NVSASACalc.class.getName()).log(Level.SEVERE, null, ex);
                continue;
            }

            CBi.setX( tmp.getX()*-1.0 );
            CBi.setY( tmp.getY()*-1.0 );
            CBi.setZ( tmp.getZ()*-1.0 );

            btree.neighbours(tmp, UBOUND, listNeighb);

            //calculate neighbor weights, neighbor counts and neighbor vector
            for( Atom neighb : listNeighb ){

                if( !neighb.getName().equals("CB") )
                    continue;

                vector.setX( neighb.getX() );
                vector.setY( neighb.getY() );
                vector.setZ( neighb.getZ() );
                CalcGeom.addEquals(vector, CBi);
                double distance=Calc.amount(vector);

                if( distance<1e-9 )
                    continue;

                nweight=neighborWeight2( distance );
                ncount+=nweight;
                CalcGeom.product2( vector, nweight/distance );
                CalcGeom.addEquals( nvector, vector);
                
            }//

            listNeighb.clear();

            // finish the neighbor vector
            CalcGeom.product2( nvector, 1.0/ncount );

            sasa+=Calc.amount(nvector)*rotLib.getAvgSASA(groups.get(i));
        }


        return sasa;
    }

    /**
     * calcs neighbor weight function
     *
     * @param dist
     * @return
     */
    static private double neighborWeight( double dist ){

        if( dist<=LBOUND )
            return 1.0;

        if( dist>LBOUND && dist<UBOUND )
            return 0.5*( Math.cos( (dist-LBOUND)*KNWEIGHT )+1.0 );

        return 0.0;
    }

    /**
     * calcs neighbor weight function
     * pre : dist < UBOUND
     *
     * @param dist
     * @return
     */
    static private double neighborWeight2( double dist ){

        if( dist<=LBOUND )
            return 1.0;

        return 0.5*( Math.cos( (dist-LBOUND)*KNWEIGHT )+1.0 );
    }
    
}
