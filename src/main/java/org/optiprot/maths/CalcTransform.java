/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.ArrayList;
import java.util.List;
import org.optiprot.aacids.IAABasic;
import org.biojava.bio.structure.*;
import org.biojava.bio.structure.jama.Matrix;
import org.optiprot.aacids.AABasicFactory;

/**
 *
 * @author victor
 */
public class CalcTransform {
    private static double ANG_ROTATION_MIN=1.0;//in degrees
    private static double ANG_LIMIT=11.0;

   
    /**
     * change phi angle of the given chain at given position
     *
     * @param chain
     * @param aaindex
     * @param angle in radians
     */
    public static void changePhi( IAABasic[] chain, int aaindex, double angle ){

        if( aaindex==0 )
            return;
        
        IAABasic aacid=chain[aaindex];

        CalcTransform.rotateSubChain( aacid.getN(), aacid.getVectorNCa(),
                angle, chain, 0, aaindex);
    }

    /**
     * change psi angle of the given chain at given position
     *
     * @param chain
     * @param aaindex
     * @param angle in radians
     */
    public static void changePsi( IAABasic[] chain, int aaindex, double angle ){

        if( aaindex==chain.length-1 )
            return;

        IAABasic aacid=chain[aaindex];

        CalcTransform.rotateSubChain( aacid.getC(), aacid.getVectorCaC(),
                angle, chain, aaindex+1, chain.length);
    }

    /**
     * Rotates a subchain given a center of rotation and an axis
     *
     * @param center
     * @param vaxis
     * @param angle in radians
     * @param genes
     * @param start
     * @param end
     */
    public static void rotateSubChain(Atom center, Atom vaxis, double angle,
            IAABasic[] genes, int start, int end) {

        vaxis=Calc.unitVector(vaxis);

        Matrix mtrans=createPointRotatMatrix( center, vaxis, angle );

        for( int i=start;i<end;i++ ){

            IAABasic aacid = genes[i];
            applyTransform( aacid.getAtoms(), mtrans);
        }

    }

    /**
     * 
     *
     * @param center
     * @param vaxis
     * @param angle in radians
     * @param chain
     * @param start
     * @param end
     */
    public static void rotateSubChain(Atom center, Atom vaxis, double angle,
            Chain chain, int start, int end) {

        vaxis=Calc.unitVector(vaxis);

        Matrix mtrans=createPointRotatMatrix( center, vaxis, angle );

        List<Group> list = chain.getAtomGroups();

        Atom[] array=null;

        for( int i=start;i<end;i++ ){

            Group aacid=list.get(i);
            array=new Atom [aacid.size()];
            applyTransform( aacid.getAtoms().toArray(array), mtrans);
        }

    }

    /**
     * translate a subchain by a shift vector
     *
     * @param shiftVector
     * @param genes
     * @param start
     * @param end
     */
    public static void translateSubChain(Atom shiftVector,
            IAABasic[] genes, int start, int end) {

        for( int i=start;i<end;i++ ){

            IAABasic aacid = genes[i];

            for( Atom atom : aacid.getAtoms() ){
               CalcGeom.addEquals(atom, shiftVector);
            }
        }

    }

    /**
     *
     * @param shiftVector
     * @param chain
     * @param start
     * @param end
     */
    public static void translateSubChain(Atom shiftVector,
            Chain chain, int start, int end) {

        List<Group> list = chain.getAtomGroups();

        for( int i=start;i<end;i++ ){

            Group aacid = list.get(i);

            for( Atom atom : aacid.getAtoms() ){
               CalcGeom.addEquals(atom, shiftVector);
            }
        }

    }

    /**
     * Transforms a subchain to meet the constraints of the peptidic bond,
     * used to cross two proteins by one cut point.
     *
     * @param genesSrc
     * @param genesDst
     * @param cutPoint
     * @param start
     * @param end
     */
    public static void placeSubChain( IAABasic[] genesSrc, IAABasic[] genesDst ,
            int cutPoint, int start, int end) throws Exception {

        IAABasic aaSrc1 =  genesSrc[cutPoint+1];
        IAABasic aaSrc2 =  genesSrc[cutPoint];
        IAABasic aaDst1 =  genesDst[cutPoint];
        IAABasic aaDst2 =  genesDst[cutPoint+1];

        Atom center1=aaSrc1.getN();//origin point
        Atom center2=aaDst2.getN();//destination point

        //calcs omega angles (around 180 degrees)
        double omega1=Calc.torsionAngle(aaDst1.getCa(), aaDst1.getC(),
                aaDst2.getN(), aaDst2.getCa());
        double omega2=Calc.torsionAngle(aaSrc2.getCa(), aaSrc2.getC(),
                aaSrc1.getN(), aaSrc1.getCa());

        omega1=(omega1+omega2)*0.5;

        // Calc the parameters for R1
        // R1 makes the torsion angle Ca-C-N-Ca joining point equals to 180

        Atom ab = Calc.substract(aaDst1.getCa(),aaDst1.getC());
        Atom cb = Calc.substract(aaDst2.getN(),aaDst1.getC());
        Atom bc = Calc.substract(aaSrc2.getC(),aaSrc1.getN());
        Atom dc = Calc.substract(aaSrc1.getCa(),aaSrc1.getN());

        Atom abc = Calc.vectorProduct(ab,cb);
        Atom bcd = Calc.vectorProduct(bc,dc);

        double angle1 =Calc.angle(abc,bcd);

        // calc the sign:
        Atom axis1 = Calc.vectorProduct(abc,bcd);
        double val = Calc.skalarProduct(cb,axis1);
        if (val<0.0) angle1 = -angle1 ;

        angle1 = angle1>0 ?  omega1-angle1 : -omega1-angle1 ;

        axis1=Calc.unitVector(axis1);
        
        Matrix t1 = Matrix.identity(4,4);
        Matrix r1=null;
        
        if( Math.abs(angle1)>CalcTransform.ANG_ROTATION_MIN ){
            angle1 = Math.toRadians(angle1);
            r1 = new Matrix( matrixFromAxisAngle(axis1, angle1) );
        }
        else{
            r1=Matrix.identity(4,4);
        }

        Matrix t2 = Matrix.identity(4,4);

        t1.set(0, 3, -center1.getX() );
        t1.set(1, 3, -center1.getY() );
        t1.set(2, 3, -center1.getZ() );

        t2.set(0, 3, center2.getX() );
        t2.set(1, 3, center2.getY() );
        t2.set(2, 3, center2.getZ() );

        Matrix mtrans=t2.times( r1.times(t1) );

        // apply the transformation to the first aminoacid of the subchain so
        // with the new coordinates we can calculate R2 parameterers

        Atom [] array= new Atom [3];

        array[0]=(Atom)aaSrc1.getN().clone();
        array[1]=(Atom)aaSrc1.getCa().clone();
        array[2]=(Atom)aaSrc1.getC().clone();

        applyTransform( array, mtrans);

        ab = Calc.substract(aaDst1.getC(),array[0]);
        cb = Calc.substract(array[1],array[0]);

        Atom axis2 = Calc.vectorProduct(ab,cb);

        axis2=Calc.unitVector(axis2);

        double angle2 = Calc.angle(cb,ab);

        //calc backbone C-N-Ca angle ( around 120 degrees)
        ab = Calc.substract(aaSrc1.getCa(),aaSrc1.getN());
        cb = Calc.substract(aaSrc2.getC(),aaSrc1.getN());
        double backCNCaangle1=Calc.angle(ab, cb);
        ab = Calc.substract(aaDst2.getCa(),aaDst2.getN());
        cb = Calc.substract(aaDst1.getC(),aaDst2.getN());
        double backCNCaangle2=Calc.angle(ab, cb);

        backCNCaangle1=(backCNCaangle1+backCNCaangle2)*0.5;

        angle2 = (angle2 > backCNCaangle1) ? -angle2+backCNCaangle1 : backCNCaangle1-angle2;


        // we have R2 parameters, build the final transformation and apply to
        // the entire subchain

        Matrix r2 =null;

        if( Math.abs(angle2)>CalcTransform.ANG_ROTATION_MIN ){
            angle2=Math.toRadians(angle2);
            r2 = new Matrix( matrixFromAxisAngle(axis2, angle2) );
        }
        else{
            r2 = Matrix.identity(4,4);
        }

        mtrans=t2.times( r2.times(r1.times(t1)) );

        for( int i=start;i<end;i++ ){

            IAABasic aacid = genesSrc[i];
            applyTransform( aacid.getAtoms(), mtrans);
        }
        


    }

    /**
     * transforms the Amino Acid to bond next to position given by target,
     * meeting the constraints of the peptidic bond
     *
     * @param aacid
     * @param genesDst
     * @param target
     * @throws java.lang.Exception
     */
    public static void addAminoAcid( IAABasic aacid, IAABasic[] genesDst ,
            int target ) throws StructureException {

        IAABasic aaDst =  genesDst[target];

        Atom center1=aacid.getN();//origin point
        Atom center2=generRandNPlace(aaDst);//destination point

        Atom shift=Calc.substract(center2, center1);

        //move aacid to the new location
        for( Atom atom : aacid.getAtoms() ){
               CalcGeom.addEquals(atom, shift);
        }

        Matrix mtrans=null;

        //normal to plane Ca-N-C'
        Atom norm = CalcGeom.normalPlane( aaDst.getCa(), aacid.getN(), aaDst.getC());
        Atom nc = Calc.substract( aaDst.getC(), aacid.getN() );
        Atom nca1 = Calc.substract( aacid.getCa(), aacid.getN() );
     
        double dist2=Calc.amount(nca1);//N-Ca bond length

        //generates C-N-Ca angle
        double backCNCaangle=119+(Math.random()-0.5)*4.0;

        //calculates the new location for Ca
        Atom nca2=Calc.unitVector( Calc.vectorProduct( norm, nc) );
        
        mtrans = createPointRotatMatrix( null, norm, Math.toRadians(backCNCaangle-90) );
        applyTransform( nca2, mtrans);
        CalcGeom.product2(nca2, dist2);

        //rotate Ca to get the new location
        Atom axis1 = Calc.vectorProduct(nca1, nca2);
        axis1=Calc.unitVector(axis1);

        double angle1=Calc.angle(nca1, nca2);

        if( Math.abs(angle1)>CalcTransform.ANG_ROTATION_MIN ){
            angle1 = Math.toRadians(angle1);
            mtrans = createPointRotatMatrix( center2, axis1, angle1 );
        }
        else{
            mtrans=Matrix.identity(4,4);
        }

        // rotate the amino acid to perform the omega angle constraint and
        //  the C-N-Ca angle
        applyTransform( aacid.getAtoms(), mtrans);

    }

    /**
     *
     * @param atoms
     * @param m : a 3*3 matrix
     */
    public static void applyTransform(Atom[] atoms, double[][] m) {

        Matrix mt=Matrix.identity(4,4);
        mt.set(0, 0, m[0][0]);
        mt.set(0, 1, m[0][1]);
        mt.set(0, 2, m[0][2]);

        mt.set(1, 0, m[1][0]);
        mt.set(1, 1, m[1][1]);
        mt.set(1, 2, m[1][2]);

        mt.set(2, 0, m[2][0]);
        mt.set(2, 1, m[2][1]);
        mt.set(2, 2, m[2][2]);

        applyTransform( atoms,  mt);
    }

    /**
     *
     * @param chain
     * @param mtrans
     * @param start
     * @param end
     */
    public static void applyTransform( Chain chain, Matrix mtrans, int start, int end) {

        Matrix t=null;

        if( mtrans.getRowDimension()==3 ){
            t = Matrix.identity(4,4);
            t.setMatrix(0, 2, 0, 2, mtrans);
        }
        else{
            t=mtrans;
        }

        List<Group> list = chain.getAtomGroups();

        for( int i=start;i<end;i++ ){

            Group aacid = list.get(i);

            for( Atom atom : aacid.getAtoms() ){
               applyTransform( atom, t);
            }
        }
    }

    /**
     *
     * @param atoms
     * @param mtrans
     */
    private static void applyTransform(Atom[] atoms, Matrix mtrans) {

        double [] column ={0,0,0,1};
        double [] result ={0,0,0,1};

        for(int i=0;i<atoms.length;i++){

            column[0]=atoms[i].getX();
            column[1]=atoms[i].getY();
            column[2]=atoms[i].getZ();

            result[0]=mtrans.get(0,0)*column[0] +
                    mtrans.get(0,1)*column[1] +
                    mtrans.get(0,2)*column[2] +
                    mtrans.get(0,3)*column[3];

            result[1]=mtrans.get(1,0)*column[0] +
                    mtrans.get(1,1)*column[1] +
                    mtrans.get(1,2)*column[2] +
                    mtrans.get(1,3)*column[3];

            result[2]=mtrans.get(2,0)*column[0] +
                    mtrans.get(2,1)*column[1] +
                    mtrans.get(2,2)*column[2] +
                    mtrans.get(2,3)*column[3];
            
            atoms[i].setX(result[0]);
            atoms[i].setY(result[1]);
            atoms[i].setZ(result[2]);
        }

    }

    public static void applyTransform(Atom atom, Matrix mtrans) {

        double [] column ={0,0,0,1};
        double [] result ={0,0,0,1};

        column[0]=atom.getX();
        column[1]=atom.getY();
        column[2]=atom.getZ();

        result[0]=mtrans.get(0,0)*column[0] +
                mtrans.get(0,1)*column[1] +
                mtrans.get(0,2)*column[2] +
                mtrans.get(0,3)*column[3];

        result[1]=mtrans.get(1,0)*column[0] +
                mtrans.get(1,1)*column[1] +
                mtrans.get(1,2)*column[2] +
                mtrans.get(1,3)*column[3];

        result[2]=mtrans.get(2,0)*column[0] +
                mtrans.get(2,1)*column[1] +
                mtrans.get(2,2)*column[2] +
                mtrans.get(2,3)*column[3];

        atom.setX(result[0]);
        atom.setY(result[1]);
        atom.setZ(result[2]);
    }


    /**
     *
     * @param center : center of the rotation (null if origin)
     * @param axis
     * @param angle in radians
     * @return
     */
    private static Matrix createPointRotatMatrix( Atom center, Atom axis, double angle ){

        Matrix t1 = Matrix.identity(4,4);
        Matrix r = new Matrix( matrixFromAxisAngle(axis, angle) );
        Matrix t2 = Matrix.identity(4,4);

        if( center!=null ){
            t1.set(0, 3, -center.getX() );
            t1.set(1, 3, -center.getY() );
            t1.set(2, 3, -center.getZ() );

            t2.set(0, 3, center.getX() );
            t2.set(1, 3, center.getY() );
            t2.set(2, 3, center.getZ() );
        }

        return t2.times( r.times(t1) );

    }


    /**
     * 
     * @param axis
     * @param angle in radians
     * @return
     */
    public static double[][] matrixFromAxisAngle(Atom axis, double angle) {

        double c = Math.cos(angle);
        double s = Math.sin(angle);
        double t = 1.0 - c;
        //  if axis is not already normalised then uncomment this
        // double magnitude = Math.sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
        // if (magnitude==0) throw error;
        // a1.x /= magnitude;
        // a1.y /= magnitude;
        // a1.z /= magnitude;

        double [][]m=new double[4][4];

        m[0][3]=m[1][3]=m[2][3]=0.0;
        m[3][0]=m[3][1]=m[3][2]=0.0;
        m[3][3]=1.0;

        m[0][0] = c + axis.getX()*axis.getX()*t;
        m[1][1] = c + axis.getY()*axis.getY()*t;
        m[2][2] = c + axis.getZ()*axis.getZ()*t;


        double tmp1 = axis.getX()*axis.getY()*t;
        double tmp2 = axis.getZ()*s;
        m[1][0] = tmp1 + tmp2;
        m[0][1] = tmp1 - tmp2;
        tmp1 = axis.getX()*axis.getZ()*t;
        tmp2 = axis.getY()*s;
        m[2][0] = tmp1 - tmp2;
        m[0][2] = tmp1 + tmp2;
        tmp1 = axis.getY()*axis.getZ()*t;
        tmp2 = axis.getX()*s;
        m[2][1] = tmp1 + tmp2;
        m[1][2] = tmp1 - tmp2;

        return m;
    }

    /**
     * copy the chain's atoms to a linear array
     *
     * @param chain
     * @return
     */
    public static Atom [] toAtomArray( Chain chain ){

        ArrayList<Atom> list=new ArrayList<Atom>();

        for( Group group : chain.getAtomGroups() ){
            for(Atom atom : group.getAtoms() ){
                list.add(atom);
            }
        }

        Atom [] molecule=list.toArray(new Atom [list.size()]);

        return molecule;
    }

    /**
     * generates a random place for a N atom bonded to
     * AminoAcid's C
     *
     *  angle choose randomly around 117 degrees
     *  bond length choose around 1.33 A
     *
     *  the place is in a circle of a plane with normal C-Ca
     *
     * @param aaDst : amino acid
     * @return
     */
    private static Atom generRandNPlace(final IAABasic aaDst) throws StructureException {

        double angle=117+(Math.random()-0.5)*8.0;
        double blength=1.33;//bond length
        double l2=Math.cos( Math.toRadians(180-angle) )*blength;
        double r2=Math.sin( Math.toRadians(180-angle) )*blength;

        //constructs the normal to the plane
        Atom vCaC=Calc.substract(aaDst.getC(),aaDst.getCa());
        double lCaC=Calc.amount(vCaC);
        vCaC=CalcGeom.product(vCaC, 1.0/lCaC);
        
        //the center of the circle on plane
        Atom center=Calc.add( aaDst.getCa(), CalcGeom.product(vCaC, lCaC+l2) );

        //pick a point of the plane (distinct than center)
        
        Atom vort=CalcGeom.genOrtVector( vCaC );

        Atom place=Calc.add(center, CalcGeom.product(vort, r2));

        
        //checks Ca-C-N
        Atom ab = Calc.substract( aaDst.getCa(), aaDst.getC() );
        Atom cb = Calc.substract( place, aaDst.getC() );
        angle=Calc.angle(ab, cb);
        
        return place;

    }

    /**
     * checks the correctness of backbone angles
     * 
     * @param chain
     * @param start
     * @param end
     * @return
     */
    public static boolean checksBackbone( IAABasic [] chain, int start, int end)
            throws StructureException {

        for( int i=start;i<end;i++ ){

            //checks N-Ca-C
            Atom ab = Calc.substract( chain[i].getN(), chain[i].getCa() );
            Atom cb = Calc.substract( chain[i].getC(), chain[i].getCa() );
            double angle=Calc.angle(ab, cb);

            if( Math.abs(angle-115) > ANG_LIMIT )
                return false;

            if( i==(end-1) )
                continue;

            //checks omega
            angle=Calc.torsionAngle( chain[i].getCa(), chain[i].getC(),
                    chain[i+1].getN(), chain[i+1].getCa() );

            if( 180-Math.abs(angle) > ANG_LIMIT*0.5 )
                return false;

            //checks Ca-C-N'
            ab = Calc.substract( chain[i].getCa(), chain[i].getC() );
            cb = Calc.substract( chain[i+1].getN(), chain[i].getC() );
            angle=Calc.angle(ab, cb);

            if( Math.abs(angle-115) > ANG_LIMIT )
                return false;

            //checks C-N'-Ca
            ab = Calc.substract( chain[i].getC(), chain[i+1].getN() );
            cb = Calc.substract( chain[i+1].getCa(), chain[i+1].getN() );
            angle=Calc.angle(ab, cb);

            if( Math.abs(angle-115) > ANG_LIMIT )
                return false;

        }//

        return true;
    }

    /**
     * center the chain at origin
     *
     * @param chain
     */
    public static void center( IAABasic [] chain ){

        Atom shiftVector=CalcGeom.getCentroid( AABasicFactory.toAtoms(chain), 1);
        CalcGeom.product2(shiftVector, -1);
        CalcTransform.translateSubChain(shiftVector, chain, 0, chain.length);
    }

    /**
     * center the chain at origin
     *
     * @param chain
     */
    public static void center( Chain chain ){

        Atom shiftVector=CalcGeom.getCentroid( chain );
        CalcGeom.product2(shiftVector, -1);
        CalcTransform.translateSubChain(shiftVector, chain, 0, chain.getAtomLength());
    }
    
}
