/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.*;
import org.biojava.bio.structure.jama.Matrix;

/**
 * formula from mathworld.wolfram.com
 *
 * @author victor
 */
public class CalcSphere {

    private static final double C_4_3_PI=Math.PI*(4.0/3.0);
    private static final double C_4_PI=Math.PI*4.0;
    private static final double C_1_3_PI=Math.PI/3.0;
    private static final double C_1_3=1.0/3.0;
    private static final double C_2_PI=2.0*Math.PI;


    static public double volume(double radius){

        return C_4_3_PI*radius*radius*radius;
    }


    static public double area(double radius){

        return C_4_PI*radius*radius;
    }

    /**
     * 
     * @param c1
     * @param r1
     * @param c2
     * @param r2
     * @return : true if exists intersection between two spheres
     */
    static public boolean isIntersection( Atom c1, double r1, Atom c2, double r2){

        double sqdist=CalcGeom.squareDistance(c1, c2);

        if (sqdist >= (r1 + r2)*(r1 + r2)) {
            return false;
        }

        return true;
    }

    static public double intersecVol( Atom c1, double r1, Atom c2, double r2){
        
        double dist=0.0;
        
        try {
            dist = Calc.getDistance(c1, c2);
        } catch (StructureException ex) {
            Logger.getLogger(CalcSphere.class.getName()).log(Level.SEVERE, null, ex);
            return 0;
        }
        
        if (dist >= r1 + r2) {
            return 0.0;
        }


        return (Math.PI*(r1+r2-dist)*(r1+r2-dist)*
                (dist*dist + 2.0*dist*r2 - 3.0*r2*r2 + 2.0*dist*r1 + 6.0*r1*r2 -3.0*r1*r1)) /
                (12.0*dist);

    }


    static public double intersecArea( Atom c1, double r1, Atom c2, double r2){

        double dist=0.0;

        try {
            dist = Calc.getDistance(c1, c2);
        } catch (StructureException ex) {
            Logger.getLogger(CalcSphere.class.getName()).log(Level.SEVERE, null, ex);
            return 0;
        }

        if (dist >= r1 + r2) {
            return 0.0;
        }

        double distInter = distInterSpheres(dist,r1,r2);

        return sphericCapArea(r1,r1-distInter)+sphericCapArea(r2,r2-(dist-distInter));
    }

    /**
     * return the dist of the intersection between 2 spheres (from (c1,r1))
     * 
     * @param dist
     * @param r1 radii of sphere at origin
     * @param r2
     * @return
     */
    private static double distInterSpheres(double dist, double r1, double r2) {

        return (dist*dist-r2*r2+r1*r1)/(2.0*dist);
    }



    static public double sphericCapArea( double radius, double h){

        return C_2_PI*radius*h;
    }

    
    static public double sphericCapVol( double radius, double h){

        return C_1_3_PI*h*h*(3.0*radius-h);
    }


    /**
     * 
     * @param c1
     * @param r1
     * @param c2
     * @param r2
     * @param c3
     * @param r3
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    static public double intersecArea( Atom c1, double r1, Atom c2, double r2,
            Atom c3, double r3, int maxDepth) throws StructureException{

        double dist12=0.0;
        double dist31=0.0;
        double dist23=0.0;

        try {
            dist12 = Calc.getDistance(c1, c2);
            dist31 = Calc.getDistance(c1, c3);
            dist23 = Calc.getDistance(c2, c3);
        } catch (StructureException ex) {
            Logger.getLogger(CalcSphere.class.getName()).log(Level.SEVERE, null, ex);
            return 0;
        }

        //--------------------------------------------
        // if any intersection of order 2 are empty
        if (dist12 >= r1 + r2) {
            return 0.0;
        }
        else if (dist31 >= r1 + r3) {
            return 0.0;
        }
        else if (dist23 >= r2 + r3) {
            return 0.0;
        }

        //calculates the contact poins between spheres (in c1-c2-c3 plane)
      
        boolean contact12=contact(r1, r2, dist12);
        boolean contact23=contact(r2, r3, dist23);
        boolean contact31=contact(r3, r1, dist31);


        //------------------------------------------------------------------
        //topological situation 1 - if theres no contact between the spheres
        if( !contact12 && !contact23 && !contact31 ){

            double rmin=0;
            rmin=Math.min(r1, r2);
            rmin=Math.min(r3, rmin);
            return area(rmin);
        }

        Atom normal=CalcGeom.normalPlane(c1, c2, c3);

        //points of contact between the i and j sphere
        Atom [] pcontact12=new Atom [2];
        Atom [] pcontact23=new Atom [2];
        Atom [] pcontact31=new Atom [2];

        Boolean [] pcontact12_i3=new Boolean [] {false,false};
        Boolean [] pcontact23_i1=new Boolean [] {false,false};
        Boolean [] pcontact31_i2=new Boolean [] {false,false};

        //number of points of contact between i-j inside the k-sphere
        int npcontact12_i3=0;
        int npcontact23_i1=0;
        int npcontact31_i2=0;

        if( contact12 ){

            npcontact12_i3 = pointsOfContact(c1, r1, c2, r2, c3, r3,
                 dist12, normal, pcontact12, pcontact12_i3 );
        }
        if( contact23 ){

            npcontact23_i1 = pointsOfContact(c2, r2, c3, r3, c1, r1,
                 dist23, normal, pcontact23, pcontact23_i1 );
        }
        if( contact31 ){

            npcontact31_i2 = pointsOfContact(c3, r3, c1, r1, c2, r2,
                 dist31, normal, pcontact31, pcontact31_i2 );
        }

        
        //-------------------------------------------------------------------
        //general situation
        if( contact12 && contact23 && contact31 ){
            if( npcontact12_i3==1 && npcontact23_i1==1 && npcontact31_i2==1 ){

                Atom [] T = trilateration( c1, r1, c2, r2, c3, r3,
                        new double [] {dist12,dist23,dist31} );

                double area=0;
                double [] results=null;

                //Xi are the contact points inside the i sphere

                Atom x1=pcontact23[0];
                if( pcontact23_i1[1] ){
                    x1=pcontact23[1];
                }

                Atom x2=pcontact31[0];
                if( pcontact31_i2[1] ){
                    x2=pcontact31[1];
                }

                Atom x3=pcontact12[0];
                if( pcontact12_i3[1] ){
                    x3=pcontact12[1];
                }
               
                results=trihedronSphereIntersect( c1, r1, T[0], x2, x3, T[1],maxDepth);
                area+=2.0*results[0];
                
                results=trihedronSphereIntersect( c2, r2, T[0], x3, x1, T[1],maxDepth);
                area+=2.0*results[0];
                
                results=trihedronSphereIntersect( c3, r3, T[0], x1, x2, T[1],maxDepth);
                area+=2.0*results[0];
                
                return area;
            }
        }

        //--------------------------------------------------------------------
        //if empty simplicia1 topology : 3 spheres empty intersection and all
        // (n-1) intersection non empty
        if( contact12 && contact23 && contact31 ){
            if( npcontact12_i3==0 && npcontact23_i1==0 && npcontact31_i2==0 ){
                return 0;
            }
        }
        
        //------------------------------------------------------------------
        // if an intersection of order 2 (i,j) is inside the remaining sphere (k)
        // and the union of i,j is inside k
        if( npcontact12_i3==2 && !contact23 && !contact31 ){

            return intersecArea(c1,r1,c2,r2);
        }
        if( npcontact23_i1==2 && !contact12 && !contact31 ){

            return intersecArea(c2,r2,c3,r3);
        }
        if( npcontact31_i2==2 && !contact23 && !contact12 ){

            return intersecArea(c3,r3,c1,r1);
        }


        //-------------------------------------------------------------
        // if a sphere i is inside the intersection between j,k
        if( contact12 && npcontact12_i3==0 && !contact23 && !contact31 ){

            return area(r3);
        }
        if( contact23 && npcontact23_i1==0 && !contact12 && !contact31 ){

            return area(r1);
        }
        if( contact31 && npcontact31_i2==0 && !contact23 && !contact12 ){

            return area(r2);
        }

        //-----------------------------------------------------------
        // if a sphere i is inside a sphere k and the intersection i,j is inside k
        if( (npcontact12_i3==2 && contact23 && !contact31)||
                (npcontact12_i3==2 && !contact23 && contact31) ){

            return intersecArea(c1,r1,c2,r2);
        }
        if( (npcontact23_i1==2 && contact12 && !contact31)||
                (npcontact23_i1==2 && !contact12 && contact31) ){

            return intersecArea(c2,r2,c3,r3);
        }
        if( (npcontact31_i2==2 && contact23 && !contact12)||
                (npcontact31_i2==2 && !contact23 && contact12) ){

            return intersecArea(c3,r3,c1,r1);
        }

        //-----------------------------------------------------
        // if the intersection i,j is inside k and k intersects with i and j
        if( npcontact12_i3==2 && contact23 && contact31 &&
                npcontact23_i1==0 && npcontact31_i2==0 ){

            return intersecArea(c1,r1,c2,r2);
        }
        if( npcontact23_i1==2 && contact31 && contact12 &&
                npcontact12_i3==0 && npcontact31_i2==0 ){

            return intersecArea(c2,r2,c3,r3);
        }
        if( npcontact31_i2==2 && contact12 && contact23 &&
                npcontact23_i1==0 && npcontact12_i3==0 ){

            return intersecArea(c3,r3,c1,r1);
        }


        //------------------------------------------------------
        // if a sphere k intersects the i,j intersection where i,j intersection
        // is not inside k
        if( npcontact12_i3==0 && npcontact23_i1==2 && npcontact31_i2==2 &&
                contact23 && contact31 ){

            return intersecArea(c3,r3,c1,r1)+intersecArea(c2,r2,c3,r3)-area(r3);
        }
        if( npcontact23_i1==0 && npcontact12_i3==2 && npcontact31_i2==2 &&
                contact12 && contact31 ){

            return intersecArea(c1,r1,c2,r2)+intersecArea(c3,r3,c1,r1)-area(r1);
        }
        if( npcontact31_i2==0 && npcontact12_i3==2 && npcontact23_i1==2 &&
                contact12 && contact23 ){

            return intersecArea(c1,r1,c2,r2)+intersecArea(c2,r2,c3,r3)-area(r2);
        }
        
        throw new StructureException("Intersect 3-Sphere - Unknown topological situation");
    }


    static private boolean contact(double r1, double r2, double dist12){

        double rmin=Math.min(r1,r2);
        double rmax=Math.max(r1,r2);

        if( Math.abs(rmax-rmin)>CalcGeom.MIN_DISTANCE && dist12<=rmax-rmin ){
            return false;
        }

        return true;
    }

    
    static private int pointsOfContact(Atom cA, double rA, Atom cB, double rB,
            Atom cC, double rC, double distAB, Atom normal,
            Atom [] pcontactAB, Boolean [] pcontactAB_inC )
            throws StructureException {

        double b=distInterSpheres( distAB, rA, rB);
        double radInter=Math.sqrt(rA*rA-b*b);

        Atom vb=CalcGeom.product(Calc.substract(cB, cA), 1.0/distAB);

        Atom va=Calc.vectorProduct(normal, vb);
        CalcGeom.product2(vb, b);
        vb=Calc.add(cA, vb);

        pcontactAB[0]=Calc.add(CalcGeom.product(va,radInter), vb);
        pcontactAB[1]=Calc.add(CalcGeom.product(va,-radInter), vb);

        //if the contact point are inside/outside the remaining sphere
        int npcontactAB_inC=0;
        double rC2=rC*rC;

        if( CalcGeom.squareDistance(pcontactAB[0], cC)<=rC2 ){
            pcontactAB_inC[0]=true;
            npcontactAB_inC++;
        }

        if( CalcGeom.squareDistance(pcontactAB[1], cC)<=rC2 ){
            pcontactAB_inC[1]=true;
            npcontactAB_inC++;
        }

        return npcontactAB_inC;
    }

    /**
     * calculates the triple intersection between 3 spheres
     *
     * @param c1 center 1
     * @param r1 radii 1
     * @param c2
     * @param r2
     * @param c3
     * @param r3
     * @param distances array of distances between centers 12,23,31
     * @return an array with t,tp and tm, t is the middle of tp-tm segment
     * @throws org.biojava.bio.structure.StructureException
     */
    static public Atom[] trilateration( Atom c1, double r1, Atom c2, double r2,
            Atom c3, double r3, double [] distances ) throws StructureException{

        double sdist12=distances[0]*distances[0];
        double sdist23=distances[1]*distances[1];
        double sdist31=distances[2]*distances[2];

        double sr1=r1*r1;
        double sr2=r2*r2;
        double sr3=r3*r3;

        //semiperimeter of the triangle c1-c2-c3
        double semiperim=0.5*(distances[0]+distances[1]+distances[2]);

        //square area of the triangle c1-c2-c3
        double ssurf123=semiperim*(semiperim-distances[0]) *
                (semiperim-distances[1]) * (semiperim-distances[2]);

        // barycentric t coordinates, relative to cl, c2, c3, are:
        double a1 = sr1*(-2*sdist23) +
                sr2*(sdist23+sdist31-sdist12) +
                sr3*(sdist23-sdist31+sdist12) +
                sdist23*(-sdist23+sdist31+sdist12);

        double a2 = sr2*(-2*sdist31) +
                sr3*(sdist31+sdist12-sdist23) +
                sr1*(sdist31-sdist12+sdist23) +
                sdist31*(-sdist31+sdist12+sdist23);

        double a3 = sr3*(-2*sdist12) +
                sr1*(sdist12+sdist23-sdist31) +
                sr2*(sdist12-sdist23+sdist31) +
                sdist12*(-sdist12+sdist23+sdist31);

        // calculate t
        Atom b1=CalcGeom.product(c1,a1);
        Atom b2=CalcGeom.product(c2,a2);
        Atom b3=CalcGeom.product(c3,a3);

        Atom t=Calc.add(b3, Calc.add(b1, b2) );
        a1=1.0/(16*ssurf123);
        CalcGeom.product2(t, a1);
        

        // distance t-tp
        a1=CalcGeom.squareDistance(c1, t);
        a1=Math.sqrt(sr1-a1);

        //calc tp & tm
        b3=CalcGeom.normalPlane(c1,c2,c3);

        CalcGeom.product2(b3, a1);
        Atom tp=Calc.add(t,b3);

        b3.getCoords()[0]=-b3.getCoords()[0];
        b3.getCoords()[1]=-b3.getCoords()[1];
        b3.getCoords()[2]=-b3.getCoords()[2];
        Atom tm=Calc.add(t,b3);

        //return t, tp & tm
        return new Atom [] {t,tp,tm};
    }

    /**
     * calculates the area and volume of the trihedron sphere intersection
     *
     * @param c  center of sphere
     * @param r  radii
     * @param z  trihedron origin
     * @param x1
     * @param x2
     * @param x3
     * @return  double array : array[0]=area, array[1]=volume
     * @throws org.biojava.bio.structure.StructureException
     */
    static private double [] trihedronSphereIntersect( Atom c, double r, Atom z,
            Atom x1, Atom x2, Atom x3, int maxDepth ) throws StructureException{

        //double invR=1.0/r;
        //normals to the sphere at Xi
        //Atom n1=CalcGeom.product( Calc.substract(x1, c), invR);
        //Atom n2=CalcGeom.product( Calc.substract(x2, c), invR);
        //Atom n3=CalcGeom.product( Calc.substract(x3, c), invR);
        
        Atom nd12=CalcGeom.normalPlane(z,x1,x2);
        Atom nd23=CalcGeom.normalPlane(z,x2,x3);
        Atom nd31=CalcGeom.normalPlane(z,x3,x1);

        //centers of the arcs
        Atom c12=CalcGeom.rayPlaneIntersec( c, nd12, nd12, z);
        Atom c23=CalcGeom.rayPlaneIntersec( c, nd23, nd23, z);
        Atom c31=CalcGeom.rayPlaneIntersec( c, nd31, nd31, z);

        //radii of circular arcs
        double sr12=CalcGeom.squareDistance(x1,c12);
        double sr23=CalcGeom.squareDistance(x2,c23);
        double sr31=CalcGeom.squareDistance(x3,c31);

        double r12=Math.sqrt(sr12);
        double r23=Math.sqrt(sr23);
        double r31=Math.sqrt(sr31);

        double ir12=1.0/r12;
        double ir23=1.0/r23;
        double ir31=1.0/r31;

        //angle of circular arcs
        double b12=Math.acos( Calc.skalarProduct(
                CalcGeom.product( Calc.substract(x1, c12), ir12),
                CalcGeom.product( Calc.substract(x2, c12), ir12)) );

        double b23=Math.acos( Calc.skalarProduct(
                CalcGeom.product( Calc.substract(x2, c23), ir23),
                CalcGeom.product( Calc.substract(x3, c23), ir23)) );
        
        double b31=Math.acos( Calc.skalarProduct(
                CalcGeom.product( Calc.substract(x3, c31), ir31),
                CalcGeom.product( Calc.substract(x1, c31), ir31)) );

      
        double area=trihedronSphereIntersect( c, r, x1, x2, x3,
            c12, c23, c31, nd12, nd23, nd31, b12, b23, b31,
            false, false, false, 0, maxDepth );

        return new double [] {area, 0};
    }

    /**
     * Calculates the spherical excess of a spherical triangle delimited by the
     * unit vectors A,B,C over a sphere of radii 1.
     *
     * @param vA
     * @param vB
     * @param vC
     * @return
     */
    static public double sphericalExcess(Atom vA, Atom vB, Atom vC){

        double cosa=Calc.skalarProduct(vB,vC);
        double cosb=Calc.skalarProduct(vA,vC);
        double cosc=Calc.skalarProduct(vA,vB);

        if( Math.abs(cosa-1.0)<1e-15 ){
            return 0.0;
        }
        else if( Math.abs(cosb-1.0)<1e-15 ){
            return 0.0;
        }
        else if( Math.abs(cosb-1.0)<1e-15 ){
            return 0.0;
        }

         if( Math.abs(CalcGeom.determinant(vA,vB,vC))<1e-15 ){
            //if coplanar
            return C_2_PI;
        }

        double sina=Math.sin(Math.acos(cosa));
        double sinb=Math.sin(Math.acos(cosb));
        double sinc=Math.sin(Math.acos(cosc));

        double A=Math.acos( (cosa-cosb*cosc)/(sinb*sinc) );
        double B=Math.acos( (cosb-cosa*cosc)/(sina*sinc) );
        double C=Math.acos( (cosc-cosb*cosa)/(sinb*sina) );

        return A+B+C-Math.PI;
    }
    
    /**
     * Volume of 3-spheres intersection, with bounding box window
     * 
     * @param c1
     * @param r1
     * @param c2
     * @param r2
     * @param c3
     * @param r3
     * @param N : iterations to do
     * @return
     */
    static public double monteCarloVol( Atom c1, double r1, Atom c2, double r2,
        Atom c3, double r3, long N){

        //calculate the bounding box
        double [] max ={0,0,0};
        double [] min ={0,0,0};

        min[0]=Math.min( c1.getX()-r1, c2.getX()-r2 );
        min[0]=Math.min( min[0], c3.getX()-r3) ;

        min[1]=Math.min( c1.getY()-r1, c2.getY()-r2 );
        min[1]=Math.min( min[1], c3.getY()-r3);

        min[2]=Math.min( c1.getZ()-r1, c2.getZ()-r2 );
        min[2]=Math.min( min[2], c3.getZ()-r3) ;

        max[0]=Math.max( c1.getX()+r1, c2.getX()+r2 );
        max[0]=Math.max( max[0], c3.getX()+r3) ;

        max[1]=Math.max( c1.getY()+r1, c2.getY()+r2 );
        max[1]=Math.max( max[1], c3.getY()+r3);

        max[2]=Math.max( c1.getZ()+r1, c2.getZ()+r2 );
        max[2]=Math.max( max[2], c3.getZ()+r3) ;


        long Ni=0;

        double [] l={(max[0]-min[0]),(max[1]-min[1]),(max[2]-min[2])};
        double W=l[0]*l[1]*l[2];//volume of window

        double r12=r1*r1;
        double r22=r2*r2;
        double r32=r3*r3;

        Atom p=new AtomImpl();

        for(long i=0;i<N;i++){

            //generate a random point inside the window
            p.setX( min[0]+l[0]*Math.random() );
            p.setY( min[1]+l[1]*Math.random() );
            p.setZ( min[2]+l[2]*Math.random() );

            if( CalcGeom.squareDistance(p,c1)<r12 && CalcGeom.squareDistance(p,c2)<r22
                && CalcGeom.squareDistance(p,c3)<r32 ){
                Ni++;
            }
        }

        return W*Ni/((double)N);
    }

    /**
     * Volume of 3-spheres intersection, with sphere centered at baricenter window
     *
     * @param c1
     * @param r1
     * @param c2
     * @param r2
     * @param c3
     * @param r3
     * @param center : baricenter coordinates
     * @param radii : radii of sphere centered at baricenter
     * @param N : montecarlo iterations
     * @return
     */
    static private double monteCarloVol( Atom c1, double r1, Atom c2, double r2,
        Atom c3, double r3, Atom center, double radii, long N){

        //calculate the bounding box
        double [] max ={0,0,0};
        double [] min ={0,0,0};

        long Ni=0;

        min[0]=center.getX()-radii;
        min[1]=center.getY()-radii;
        min[2]=center.getZ()-radii;

        max[0]=center.getX()+radii;
        max[1]=center.getY()+radii;
        max[2]=center.getZ()+radii;

        double [] l={(max[0]-min[0]),(max[1]-min[1]),(max[2]-min[2])};
        double W=l[0]*l[1]*l[2];//volume of window

        double r12=r1*r1;
        double r22=r2*r2;
        double r32=r3*r3;

        Atom p=new AtomImpl();

        for(long i=0;i<N;i++){

            //generate a random point inside the window
            p.setX( min[0]+l[0]*Math.random() );
            p.setY( min[1]+l[1]*Math.random() );
            p.setZ( min[2]+l[2]*Math.random() );

            if( CalcGeom.squareDistance(p,c1)<r12 && CalcGeom.squareDistance(p,c2)<r22
                && CalcGeom.squareDistance(p,c3)<r32 ){
                Ni++;
            }
        }

        return W*Ni/((double)N);
    }

    /**
     * Surface of 3-spheres intersection
     * @param c1
     * @param r1
     * @param c2
     * @param r2
     * @param c3
     * @param r3
     * @param N
     * @return
     */
    static public double monteCarloSurf( Atom c1, double r1, Atom c2, double r2,
        Atom c3, double r3, long N){

        long Ni=0;

        double T=area(r1)+area(r2)+area(r3);
        double [] S={area(r1),area(r1)+area(r2),T};

        Atom [] centers={c1,c2,c3};
        double [] radii={r1,r2,r3};
        double [] radii2={r1*r1,r2*r2,r3*r3};

        Atom p;

        for(long i=0;i<N;i++){

            double val=Math.random()*T;
            int sphere=2;

            if( val < S[0] )
                sphere=0;
            else if( val < S[1])
                sphere=1;

            //generate a random point in the surface of the sphere
            p=CalcGeom.marsagliaUnitVectGen(radii[sphere]);

            p=Calc.add(centers[sphere], p );

            int a=(sphere+1)%3;
            int b=(sphere+2)%3;

            //if the point is inside the remaining spheres
            if( CalcGeom.squareDistance(p,centers[a])<radii2[a]
                    && CalcGeom.squareDistance(p,centers[b])<radii2[b] ){
                Ni++;
            }
        }

        return T*Ni/((double)N);
    }


    /**
     *
     * divides the trihedron sphere intersection surface recursively and measures
     * the area of interior triangles as spherical triangles
     *
     * @param n : depth
     * @param max : max depth
     * @return  area of triangle
     * @throws org.biojava.bio.structure.StructureException
     */
    static private double trihedronSphereIntersect( Atom c, final double R,
        Atom x1, Atom x2, Atom x3,
        Atom c12, Atom c23, Atom c31,
        Atom nd12, Atom nd23, Atom nd31,
        double b12, double b23, double b31,
        boolean int12, boolean int23, boolean int31,
        int n, final int max ) throws StructureException {

        double invR=1.0/R;

        if(n==max || (int12 && int23 && int31) ){
            //base case
            return R*R*sphericalExcess( CalcGeom.product( Calc.substract(x1,c),invR),
                    CalcGeom.product( Calc.substract(x2,c), invR),
                    CalcGeom.product( Calc.substract(x3,c), invR) );

        }

        //calc midpoints

        //angles
        b12*=0.5;
        b23*=0.5;
        b31*=0.5;

        Atom m12=multPointRotatMatrix(x2, c12, nd12, b12);
        Atom m23=multPointRotatMatrix(x3, c23, nd23, b23);
        Atom m31=multPointRotatMatrix(x1, c31, nd31, b31);

        Atom nm12=CalcGeom.product( Calc.substract(m12,c),invR);
        Atom nm23=CalcGeom.product( Calc.substract(m23,c),invR);
        Atom nm31=CalcGeom.product( Calc.substract(m31,c),invR);

        double area=R*R*sphericalExcess( nm12, nm23, nm31 );

        Atom ndtr1=CalcGeom.normalPlane(c,m12,m31);
        Atom ndtr2=CalcGeom.normalPlane(c,m23,m12);
        Atom ndtr3=CalcGeom.normalPlane(c,m31,m23);

        double ang1=Math.acos( Calc.skalarProduct(nm12, nm31) );
        double ang2=Math.acos( Calc.skalarProduct(nm12, nm23) );
        double ang3=Math.acos( Calc.skalarProduct(nm31, nm23) );

        return area + trihedronSphereIntersect( c, R,
            x1, m12, m31,
            c12, c, c31,
            nd12, ndtr1, nd31,
            b12, ang1, b31,
            int12, true, int31,
            n+1, max) +

        trihedronSphereIntersect( c, R,
            m12, x2, m23,
            c12, c23, c,
            nd12, nd23, ndtr2,
            b12, b23,  ang2,
            int12, int23, true,
            n+1, max) +

        trihedronSphereIntersect( c, R,
            m31, m23, x3,
            c, c23, c31,
            ndtr3, nd23, nd31,
            ang3, b23, b31,
            true, int23, int31,
            n+1, max);

    }

    /**
     * rotate a point around a center
     *
     * @param point
     * @param center
     * @param axis
     * @param angle
     * @return
     */
    static private Atom multPointRotatMatrix( Atom point, Atom center, Atom axis, double angle ){

        Matrix t1=Matrix.identity(4,4);
        Matrix r = new Matrix( CalcTransform.matrixFromAxisAngle(axis, angle) );
        Matrix t2=Matrix.identity(4,4);

        t1.set(0, 3, -center.getX() );
        t1.set(1, 3, -center.getY() );
        t1.set(2, 3, -center.getZ() );

        t2.set(0, 3, center.getX() );
        t2.set(1, 3, center.getY() );
        t2.set(2, 3, center.getZ() );

        Atom v=new AtomImpl();

        v.setX( point.getX() );
        v.setY( point.getY() );
        v.setZ( point.getZ() );

        CalcTransform.applyTransform(v,t1);
        CalcTransform.applyTransform(v,r);
        CalcTransform.applyTransform(v,t2);

        return v;
    }
    
}
