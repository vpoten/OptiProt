/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.optiprot.aacids.AABasicCaCbNC;
import org.optiprot.aacids.IAABasic;
import org.optiprot.aacids.AABasicFactory;
import org.optiprot.rotamer.RotamerLibrary;

/**
 *
 * @author victor
 */
public class CalcClashes {

    static double _1_4_PI=1.0/(Math.PI*4.0);
    static double radCA=3;
    static double SQFAR_DIST=17*17;
    static double radAtom=2;
    static double acceptRate=0.25;
    static int mutRotamerLim=3;


    //////////////////////////////////////

    /**
     * crisp clash detector
     *
     * @param genes
     * @param start
     * @param end
     * @param rot
     * @return : true if clashed
     */
    public static boolean isSelfClash( IAABasic[] chain, int start, int end, RotamerLibrary rot ){

        ArrayList<IAABasic> list=new ArrayList<IAABasic>();

        for(int i=0;i<start;i++){
            list.add( chain[i] );
        }

        for(int i=end;i<chain.length;i++){
            list.add( chain[i] );
        }

        double sqradCA=radCA*radCA;
        double sqradCM=0.0;

        for( int i=start;i<end;i++){

            IAABasic aacid=chain[i];
            Atom CA=aacid.getCa();
            Atom CoM=((AABasicCaCbNC)aacid).getCoM();

            sqradCM= _1_4_PI*rot.getAvgSideVdWSA( aacid.getName() );

            for( IAABasic aatest : list ){

                double sqdistInterCA=CalcGeom.squareDistance(CA, aatest.getCa());

                if( sqdistInterCA > SQFAR_DIST )
                    continue;

                if( sqdistInterCA < sqradCA )
                    return true;

                double sqradCMtest=_1_4_PI*rot.getAvgSideVdWSA( aatest.getName() );

                if( CalcGeom.squareDistance(CA, ((AABasicCaCbNC)aatest).getCoM() ) <
                        sqradCMtest )
                    return true;

                if( CalcGeom.squareDistance(CoM, aatest.getCa()) < sqradCM )
                    return true;

                if( CalcGeom.squareDistance(CoM, ((AABasicCaCbNC)aatest).getCoM()) <
                        Math.min(sqradCM,sqradCMtest) )
                    return true;

            }
        }//


        return false;
    }

    /**
     * evaluate the membership function for fuzzy set "Atom Clash"
     *
     * @param dist : distance
     * @param distLim : limit distance
     * @param acceptRate : distance acceptability factor
     * @param lambda : acceptability param
     * @return
     */
    static protected boolean isClashed( double dist, double distLim,
            double acceptRate, double lambda){

        if( dist >= distLim )
            return false;

        double alpha = distLim*(1.0-acceptRate);

        if( (dist < distLim) && (dist > alpha) ){
            double val = (distLim-dist)/(distLim-alpha);

            if( val<lambda )
                return false;
        }

        return true;
    }


    /**
     * 
     * @param lambda
     * @param chain
     * @param rot
     * @param maxTrials
     */
    static public void fixClashes( double lambda, IAABasic[] chain, 
            RotamerLibrary rot, int maxTrials ){

        ///double lambdaDec=(lambda-0.25)/(double)maxTrials;

        while( maxTrials>0 ){

            //detect clashes
            List<AAClash> clashes = clashDetection( lambda, chain, rot );

            if( clashes.size()==0 )
                return;

            //sort in reverse order
            Collections.sort(clashes, new AAClashComparator() );

            while( !clashes.isEmpty() ){

                //choose a clash to fix (the closest)
                AAClash clashToFix=clashes.get( clashes.size()-1 );
                
                boolean fixed=modifyFix( chain, clashToFix, lambda );

                clashes.remove( clashes.size()-1 );
            }

            maxTrials--;
            ///lambda-=lambdaDec;
        }
    }

    /**
     * fuzzy detection of the chain's clashes
     *
     * @param lambda
     * @param chain
     * @param rot
     * @return : the list of clashes
     */
    static public List<AAClash> clashDetection(double lambda, IAABasic[] chain, RotamerLibrary rot) {

        List<AAClash> clashes=new ArrayList<AAClash>();

        double sqradCA=radCA*radCA;
        double sqradCM=0.0;

        for( int i=0;i<chain.length-1;i++){

            IAABasic aacid=chain[i];
            Atom CA=aacid.getCa();
            Atom CoM=((AABasicCaCbNC)aacid).getCoM();

            sqradCM= _1_4_PI*rot.getAvgSideVdWSA( aacid.getName() );

            for( int j=i+1; j<chain.length; j++ ){

                IAABasic aatest=chain[j];

                double sqdistInterCA=CalcGeom.squareDistance(CA, aatest.getCa());

                if( sqdistInterCA > SQFAR_DIST )
                    continue;

                if( isClashed( sqdistInterCA, sqradCA, acceptRate, lambda) ){
                    clashes.add( new AAClash( i, CA, j, aatest.getCa(), sqradCA) );
                    continue;
                }

                double sqradCMtest=_1_4_PI*rot.getAvgSideVdWSA( aatest.getName() );

                Atom CoMtest=((AABasicCaCbNC)aatest).getCoM();

                double sqdist=CalcGeom.squareDistance(CA, CoMtest );

                if( isClashed( sqdist, sqradCMtest, acceptRate, lambda) ){
                    clashes.add( new AAClash( i, CA, j, CoMtest, sqradCMtest) );
                    continue;
                }

                sqdist=CalcGeom.squareDistance(CoM, aatest.getCa());

                if( isClashed( sqdist, sqradCM, acceptRate, lambda) ){
                    clashes.add( new AAClash( i, CoM, j, aatest.getCa(), sqradCM) );
                    continue;
                }

                sqdist=CalcGeom.squareDistance(CoM, CoMtest);

                if( isClashed( sqdist, Math.min(sqradCM,sqradCMtest), acceptRate, lambda) ){
                    clashes.add( new AAClash( i, CoM, j, CoMtest, Math.min(sqradCM,sqradCMtest)) );
                    continue;
                }
               
            }
        }//

        return clashes;
    }

    
    /**
     * modify the chain until the clash is fixed, if the clash couldnt be fixed
     * the chain is unmodified.
     *
     * @param chain
     * @param clash
     * @return : true if the clash is fixed
     */
    static protected boolean modifyFix( IAABasic[] chain, AAClash clash, double lambda ){

        boolean clashed = clash.isClashed(lambda);

        if( !clashed )
            return true;

        IAABasic[] copyChain =  AABasicFactory.clone(chain);

        int aatarget=clash.idx1+(clash.distance()/2);


        AAClash clash2=new AAClash( clash.idx1,
                copyChain[clash.idx1].getAtom( clash.atom1.getName() ),
                clash.idx2,
                copyChain[clash.idx2].getAtom( clash.atom2.getName() ),
                clash.sqDistLim);

        //change the rotamer to fix the clash
        Integer rotIdx1=null;
        Integer rotIdx2=null;
        
        if( clash2.atom1.getName().toUpperCase().equals("COM") ){
            rotIdx1=copyChain[clash2.idx1].getRotIdx();
        }

        if( clash2.atom2.getName().toUpperCase().equals("COM") ){
            rotIdx2=copyChain[clash2.idx2].getRotIdx();
        }

        if( rotIdx1!=null ){

            for(int i=0;i<mutRotamerLim;i++){
                copyChain[clash2.idx1].applyMutationChangeRotamer();
                clashed=clash2.isClashed(lambda);
                
                if( !clashed )
                    break;
            }

            if( !clashed ){
                chain[clash.idx1].changeRotIdx(
                        copyChain[clash2.idx1].getRotIdx() );
                return !clashed;
            }
        }

        if( rotIdx2!=null ){

            for(int i=0;i<mutRotamerLim;i++){
                copyChain[clash2.idx2].applyMutationChangeRotamer();
                clashed=clash2.isClashed(lambda);

                if( !clashed )
                    break;
            }

            if( !clashed ){
                chain[clash.idx2].changeRotIdx(
                        copyChain[clash2.idx2].getRotIdx() );
                return !clashed;
            }
        }

        //change the angle to fix the clash
        int iter=1;
        double step=Math.toRadians(1.0);//angle step
        double angle=0;
        boolean psi=false;

        //choose the backbone angle to modify
        if( Math.random()<0.5 )
            psi=true;

        while( clashed ){

            double sign=1.0;

            if( (iter%2)==1 )
                sign=-1.0;

            angle=iter*sign*step;

            if( psi )
                CalcTransform.changePsi( copyChain, aatarget, angle);
            else
                CalcTransform.changePhi( copyChain, aatarget, angle);

            iter++;

            clashed=clash2.isClashed(lambda);

            if( iter>20 )
                break;
        }

        if( !clashed ){
            //if the clash is fixed aply the modification to the original chain
            if( psi )
                CalcTransform.changePsi( chain, aatarget, angle);
            else
                CalcTransform.changePhi( chain, aatarget, angle);
        }

        copyChain=null;
        clash2=null;

        return !clashed;
    }
    
}
