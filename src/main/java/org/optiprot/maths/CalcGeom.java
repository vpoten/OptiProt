/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.*;

/**
 *
 * @author victor
 */
public class CalcGeom {

    private static final double C_1_6=1.0/6.0;
    public static final double MIN_DISTANCE=1e-20;



    /**
     *
     * @param at1
     * @param at2
     * @return
     */
    public static double squareDistance( Atom at1, Atom at2 ){

        
        double a=at1.getX()-at2.getX();
        double b=at1.getY()-at2.getY();
        double c=at1.getZ()-at2.getZ();

        return a*a+b*b+c*c;
    }

    /**
     * Returns the center  of mass of the set of atoms.
     *
     * @param atomSet
     * @param step  number of atoms to skip
     * @return
     */
    public static Atom getCentroid(Atom[] atomSet, int step ){

        double[] coords = new double[3];

        coords[0] = 0;
        coords[1] = 0;
        coords[2] = 0 ;


        for (int i =0 ; i < atomSet.length; i+=step){
            Atom a = atomSet[i];
            coords[0] += a.getX();
            coords[1] += a.getY();
            coords[2] += a.getZ();
        }

        double n = 1.0/(atomSet.length/step);
        coords[0] = coords[0] * n;
        coords[1] = coords[1] * n;
        coords[2] = coords[2] * n;

        Atom vec = new AtomImpl();
        vec.setCoords(coords);
        return vec;

    }

    /**
     * Returns the center  of mass of a chain.
     *
     * @param chain
     * @return
     */
    public static Atom getCentroid(Chain chain){

        double[] coords = new double[3];

        coords[0] = 0;
        coords[1] = 0;
        coords[2] = 0 ;

        int length=0;

        for ( Group group : chain.getAtomGroups() ){
            for( Atom atom : group.getAtoms() ){
                coords[0] += atom.getX();
                coords[1] += atom.getY();
                coords[2] += atom.getZ();
                length++;
            }
        }

        double n = 1.0/length;

        coords[0] = coords[0] * n;
        coords[1] = coords[1] * n;
        coords[2] = coords[2] * n;

        Atom vec = new AtomImpl();
        vec.setCoords(coords);
        return vec;

    }

    /**
     * multiplies atom by an skalar
     *
     * @param a
     * @param skalar
     * @return a new atom
     */
    public static Atom product(Atom a, double skalar){

        Atom b=(Atom) a.clone();

        b.getCoords()[0]*=skalar;
        b.getCoords()[1]*=skalar;
        b.getCoords()[2]*=skalar;

        return b;
    }

    /**
     * multiplies atom by an skalar
     *
     * @param a
     * @param skalar
     */
    public static void product2(Atom a, double skalar){

        a.getCoords()[0]*=skalar;
        a.getCoords()[1]*=skalar;
        a.getCoords()[2]*=skalar;

    }

    /**
     * vector product res = a*b
     *
     * @param res
     * @param a
     * @param b
     */
    public static void vectorProduct(Atom res, Atom a, Atom b){

        res.setX( a.getY() * b.getZ() - a.getZ() * b.getY() );
        res.setY( a.getZ() * b.getX() - a.getX() * b.getZ() );
        res.setZ( a.getX() * b.getY() - a.getY() * b.getX() );

    }

    /**
     *
     * @param org  ray origin
     * @param vector  ray direction
     * @param normal  plane's normal
     * @param x0  a point of the plane
     * @return
     */
    public static Atom rayPlaneIntersec( Atom org, Atom vector, Atom normal, Atom x0)
            throws StructureException {

        double cos=Calc.skalarProduct(vector,normal);

        if( Math.abs(cos)<MIN_DISTANCE)
            throw new StructureException("Ray and plane are paralell");
        
        double d=-Calc.skalarProduct(normal, x0);

        double t =
                -(Calc.skalarProduct(org,normal)+d)/cos;

        return Calc.add(org, CalcGeom.product(vector, t) );
    }

    /**
     * return the normal of the plane given three points
     * (x3-x1)*(x2-x1)
     * 
     * @param x1
     * @param x2
     * @param x3
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static Atom normalPlane( Atom x1, Atom x2, Atom x3)
            throws StructureException{

        Atom b1=Calc.substract(x2,x1);
        Atom b2=Calc.substract(x3,x1);

        Atom b3=Calc.vectorProduct(b2, b1);

        if( Calc.amount(b3)<MIN_DISTANCE ){
            //if the points are collinear
            Atom a=new AtomImpl();
            
            double lb1=Calc.amount(b1);
            double lb2=Calc.amount(b2);
            
            if( lb1<MIN_DISTANCE && lb2<MIN_DISTANCE ){
                a.setCoords(new double [] {1,0,0});
            }
            else if( lb2<MIN_DISTANCE ){
                return randomPlane(x1,x2);
            }
            else{
                return randomPlane(x1,x3);
            }
        }

        return Calc.unitVector(b3);
    }

    /**
     * 
     * @param vector
     * @return
     */
    public static double squareLength(Atom vector) {
        return Calc.skalarProduct(vector, vector);
    }

    /**
     * 
     * @param z
     * @param x1
     * @param x2
     * @param x3
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    static public double tetrahedronVol(Atom z, Atom x1, Atom x2, Atom x3) 
            throws StructureException {

       return C_1_6*Math.abs( Calc.skalarProduct(
               Calc.substract(x1, z),
               Calc.vectorProduct(Calc.substract(x2, z), Calc.substract(x3, z)))
               );
    }
    
    /**
     * calcs a 3*3 determinant
     * 
     * @param a
     * @param b
     * @param c
     * @return
     */
    static public double determinant(Atom va, Atom vb, Atom vc){

        double [] a=va.getCoords();
        double [] b=vb.getCoords();
        double [] c=vc.getCoords();

        return a[0]*b[1]*c[2] - a[0]*b[2]*c[1] - a[1]*b[0]*c[2] +
                a[1]*b[2]*c[0] + a[2]*b[0]*c[1] - a[2]*b[1]*c[0];
    }

    /**
     * returns the normal to a plane that contains x1 & x2
     * 
     * @param x1
     * @param x2
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    private static Atom randomPlane(Atom x1, Atom x2) throws StructureException {

        Atom x3=new AtomImpl();
        x3.setX( Math.random() );
        x3.setY( Math.random() );
        x3.setZ( Math.random() );

        Atom b1=Calc.substract(x2,x1);
        Atom b2=Calc.substract(x3,x1);
        Atom b3=Calc.vectorProduct(b2, b1);

        while( Calc.amount(b3)<MIN_DISTANCE ){

            x3.setX( Math.random() );
            x3.setY( Math.random() );
            x3.setZ( Math.random() );

            b1=Calc.substract(x2,x1);
            b2=Calc.substract(x3,x1);
            b3=Calc.vectorProduct(b2, b1);
        }

        return Calc.unitVector(b3);
    }


    /**
     * generates points on the surface of a sphere uniformly distributed
     *
     * @param radii
     * @return
     */
    public static Atom marsagliaUnitVectGen(double radii){

        while(true){

            double x1=Math.random()*2.0-1.0;
            double x2=Math.random()*2.0-1.0;

            double sx1=x1*x1;
            double sx2=x2*x2;

            if( sx1+sx2 >=1.0 )
                continue;

            double sq=Math.sqrt(1.0-sx1-sx2);

            Atom at=new AtomImpl();
            at.setX( 2.0*x1*sq*radii );
            at.setY( 2.0*x2*sq*radii );
            at.setZ( (1.0-2.0*(sx1+sx2))*radii );

            return at;
        }
    }
    
    /**
     *  v1 += v2
     * @param atom
     * @param vector
     */
    public static void addEquals(Atom v1, Atom v2) {

        v1.setX( v1.getX()+v2.getX() );
        v1.setY( v1.getY()+v2.getY() );
        v1.setZ( v1.getZ()+v2.getZ() );
    }

    /**
     * generates a unit vector ortogonal to the given
     *
     * @param vector
     * @return
     */
    public static Atom genOrtVector( final Atom vector ){

        //choose a coordinate with value not zero
        int cnotz=0;

        if( Math.abs(vector.getX())>1e-6 ){
            cnotz=0;
        }
        else if( Math.abs(vector.getY())>1e-6 ){
            cnotz=1;
        }
        else{
            cnotz=2;
        }

        //generates a vector (x,y,z) ortogonal to vector (a,b,c)
        // ax+by+cz = 0

        double [] sol=new double [3];

        sol[(cnotz+1)%3] = (Math.random()-0.5)*2.0;
        sol[(cnotz+2)%3] = (Math.random()-0.5)*2.0;

        while( Math.abs(sol[(cnotz+2)%3])<1e-6 ){
            sol[(cnotz+2)%3] = (Math.random()-0.5)*2.0;
        }

        double [] coords= vector.getCoords();

        sol[cnotz]= ( -sol[(cnotz+1)%3]*coords[(cnotz+1)%3] -
                sol[(cnotz+2)%3]*coords[(cnotz+2)%3] ) / coords[cnotz];

        Atom vort=new AtomImpl();
        vort.setCoords(sol);

        return Calc.unitVector(vort);
    }

    /**
     * Faster pow function a^b
     * 
     * @param a
     * @param b
     * @return
     */
    public static double pow(final double a, final double b) {
        final int x = (int) (Double.doubleToLongBits(a) >> 32);
        final int y = (int) (b * (x - 1072632447) + 1072632447);
        return Double.longBitsToDouble(((long) y) << 32);
    }
    
}
