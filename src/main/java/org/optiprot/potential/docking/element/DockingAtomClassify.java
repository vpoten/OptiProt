/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.element;

import java.util.HashMap;
import org.optiprot.io.mol2.Mol2Atom;

/**
 *
 * @author victor
 */
public class DockingAtomClassify {

    //simpletype classification
    public enum simpletype {donorh,donor,donor_accept,acceptor,metal,nonpolar,nonpolarh};

    static public HashMap<Integer,String> metals=new HashMap<Integer,String>();

    static{
        //fill metals table (atomic number)
        metals.put( 3, "");//Li
        metals.put( 4, "");
        metals.put( 5, "");

        metals.put( 11, "");//Na
        metals.put( 12, "");//Mg
        metals.put( 13, "");
        metals.put( 14, "");//Si

        metals.put( 19, "");//K
        metals.put( 20, "");//Ca
        metals.put( 21, "");
        metals.put( 22, "");
        metals.put( 23, "");
        metals.put( 24, "");
        metals.put( 25, "");
        metals.put( 26, "");//Fe
        metals.put( 27, "");
        metals.put( 28, "");
        metals.put( 29, "");
        metals.put( 30, "");
        metals.put( 31, "");
        metals.put( 32, "");
        metals.put( 33, "");

        metals.put( 42, "");//Mo
        metals.put( 50, "");//Sn
    }


//    static simpletype getSimpletype( OBAtom at ){
//
//        if( at.IsHbondDonorH() )
//            return simpletype.donor;
//        else if( at.IsHydrogen() )
//            return simpletype.nonpolar;
//
//        if( at.IsHbondAcceptor() )
//            return simpletype.acceptor;
//
//        if( at.IsHbondDonor() )
//            return simpletype.donor;
//
//        if( at.IsCarbon() )
//            return simpletype.nonpolar;
//
//        if( metals.containsKey( at.GetAtomicNum() ) )
//            return simpletype.metal;
//
//        return simpletype.nonpolar;
//
//    }

    static public simpletype getSimpletype( Mol2Atom at ){

        if( at.isHbondDonorH() )
            return simpletype.donorh;
        else if( at.isHydrogen() )
            return simpletype.nonpolarh;

        if( at.isHbondAcceptor() && at.isHbondDonor() )
            return simpletype.donor_accept;
        
        if( at.isHbondAcceptor() )
            return simpletype.acceptor;

        if( at.isHbondDonor() )
            return simpletype.donor;

        if( at.isCarbon() )
            return simpletype.nonpolar;

        if( metals.containsKey( at.getAtomicNum() ) )
            return simpletype.metal;

        return simpletype.nonpolar;

    }

    
}
