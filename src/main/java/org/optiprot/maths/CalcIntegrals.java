/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.optiprot.potential.IForceField;

/**
 *
 * @author victor
 */
public class CalcIntegrals {

    private static double GRID_RESOL=0.5;//in angstroms
    private static double GRID_RESOL_INC=0.2;//in angstroms

    private static double MAX_ATOM_RADII=2.0;//in angstroms
    private static double WATER_RADII=1.4;

    static double _1_4_PI=1.0/(4.0*Math.PI);
    static double _4_PI=4.0*Math.PI;

   /**
    * integrates r-4 outside the rectangular region x1 < x < x2,
    * y1< y < y2, z1 < z < z2
    * from the paper:
    * "Potential energy functions for protein design"
    *
    * @param x array x1,x2
    * @param y
    * @param z
    * @return
    */
    static private double intOutRectangularR4( double [] x, double []  y,
            double []  z){

        double sum=0;

        sum+= -auxOutRectangularG4( x[0], y[0], z[0], z[1]);
        sum+= auxOutRectangularG4( x[0], y[1], z[0], z[1]);

        sum+= auxOutRectangularG4( x[1], y[0], z[0], z[1]);
        sum+= -auxOutRectangularG4( x[1], y[1], z[0], z[1]);

        sum+= -auxOutRectangularG4( x[0], z[0], y[0], y[1]);
        sum+= auxOutRectangularG4( x[0], z[1], y[0], y[1]);

        sum+= auxOutRectangularG4( x[1], z[0], y[0], y[1]);
        sum+= -auxOutRectangularG4( x[1], z[1], y[0], y[1]);

        sum+= -auxOutRectangularG4( y[0], z[0], x[0], x[1]);
        sum+= auxOutRectangularG4( y[0], z[1], x[0], x[1]);

        sum+= auxOutRectangularG4( y[1], z[0], x[0], x[1]);
        sum+= -auxOutRectangularG4( y[1], z[1], x[0], x[1]);
           

        return sum;
    }


    static private double auxOutRectangularG4( double a, double b,
            double c1, double c2 ){

        double ab2=Math.sqrt(a*a+b*b);

        return ab2*((Math.atan(c1/ab2))-(Math.atan(c2/ab2)))/(2*a*b) ;

    }

    /**
     * integrates r-5 outside the rectangular region x1 < x < x2,
     * y1< y < y2, z1 < z < z2
     *from the paper:
    * "Potential energy functions for protein design"
     *
     * @param x  array x1,x2
     * @param y
     * @param z
     * @return
     */
    static private double intOutRectangularR5( double [] x, double []  y,
            double []  z){

        double sum=0;

        sum+= -Math.sqrt(x[0]*x[0]+y[0]*y[0]+z[0]*z[0])/
                (x[0]*y[0]*z[0]);

        sum+= Math.sqrt(x[0]*x[0]+y[0]*y[0]+z[1]*z[1])/
                (x[0]*y[0]*z[1]);

        sum+= Math.sqrt(x[0]*x[0]+y[1]*y[1]+z[0]*z[0])/
                (x[0]*y[1]*z[0]);

        sum+= -Math.sqrt(x[0]*x[0]+y[1]*y[1]+z[1]*z[1])/
                (x[0]*y[1]*z[1]);

        sum+= Math.sqrt(x[1]*x[1]+y[0]*y[0]+z[0]*z[0])/
                (x[1]*y[0]*z[0]);

        sum+= -Math.sqrt(x[1]*x[1]+y[0]*y[0]+z[1]*z[1])/
                (x[1]*y[0]*z[1]);

        sum+= -Math.sqrt(x[1]*x[1]+y[1]*y[1]+z[0]*z[0])/
                (x[1]*y[1]*z[0]);

        sum+= Math.sqrt(x[1]*x[1]+y[1]*y[1]+z[1]*z[1])/
                (x[1]*y[1]*z[1]);
        

        double sum1=sum/6.0;
        sum=0;
        

        sum+= -auxOutRectangularG5( x[0], y[0], z[0]);
        sum+= auxOutRectangularG5( x[0], y[0], z[1]);

        sum+= -auxOutRectangularG5( x[0], z[0], y[0]);
        sum+= auxOutRectangularG5( x[0], z[0], y[1]);

        sum+= -auxOutRectangularG5( y[0], z[0], x[0]);
        sum+= +auxOutRectangularG5( y[0], z[0], x[1]);


        sum+= auxOutRectangularG5( x[0], y[1], z[0]);
        sum+= -auxOutRectangularG5( x[0], y[1], z[1]);

        sum+= auxOutRectangularG5( x[0], z[1], y[0]);
        sum+= -auxOutRectangularG5( x[0], z[1], y[1]);

        sum+= auxOutRectangularG5( y[0], z[1], x[0]);
        sum+= -auxOutRectangularG5( y[0], z[1], x[1]);


        sum+= auxOutRectangularG5( x[1], y[0], z[0]);
        sum+= -auxOutRectangularG5( x[1], y[0], z[1]);

        sum+= auxOutRectangularG5( x[1], z[0], y[0]);
        sum+= -auxOutRectangularG5( x[1], z[0], y[1]);

        sum+= auxOutRectangularG5( y[1], z[0], x[0]);
        sum+= -auxOutRectangularG5( y[1], z[0], x[1]);


        sum+= -auxOutRectangularG5( x[1], y[1], z[0]);
        sum+= auxOutRectangularG5( x[1], y[1], z[1]);

        sum+= -auxOutRectangularG5( x[1], z[1], y[0]);
        sum+= auxOutRectangularG5( x[1], z[1], y[1]);

        sum+= -auxOutRectangularG5( y[1], z[1], x[0]);
        sum+= auxOutRectangularG5( y[1], z[1], x[1]);

        sum/=6.0;

        return sum1+sum;
    }

    
    private static double auxOutRectangularG5(double a, double b, double c) {
        
        return (Math.atan((a*b)/(c*Math.sqrt(a*a+b*b+c*c))))/(c*c);
    }

    /**
     * generates a grid for the volume integral calculation
     * 
     * @param volume
     * @return
     */
    public static AtomGrid generateGrid(BSPTree volume, IForceField ffield, AtomGrid grid){

        Atom [] bbox=volume.getBBox();

        double gridResol=GRID_RESOL;

        int npointx=(int) ((bbox[1].getX() - bbox[0].getX() + 2*(MAX_ATOM_RADII+WATER_RADII)) / gridResol);
        int npointy=(int) ((bbox[1].getY() - bbox[0].getY() + 2*(MAX_ATOM_RADII+WATER_RADII)) / gridResol);
        int npointz=(int) ((bbox[1].getZ() - bbox[0].getZ() + 2*(MAX_ATOM_RADII+WATER_RADII)) / gridResol);

        //adjust grid resolution if the dimension is bigger
        while( npointx*npointy*npointz > grid.size() ){

            gridResol+=GRID_RESOL_INC;
            npointx=(int) ((bbox[1].getX() - bbox[0].getX() + 2*(MAX_ATOM_RADII+WATER_RADII)) / gridResol);
            npointy=(int) ((bbox[1].getY() - bbox[0].getY() + 2*(MAX_ATOM_RADII+WATER_RADII)) / gridResol);
            npointz=(int) ((bbox[1].getZ() - bbox[0].getZ() + 2*(MAX_ATOM_RADII+WATER_RADII)) / gridResol);
        }

        double addLength=MAX_ATOM_RADII+WATER_RADII;

        grid.clear();
        grid.setDimension(npointx,npointy,npointz);
        grid.setResol(gridResol);

        double minX=bbox[0].getX()-addLength;
        double minY=bbox[0].getY()-addLength;
        double minZ=bbox[0].getZ()-addLength;
        Atom point=new AtomImpl();

        for(int i=0;i<npointx;i++){
            for(int j=0;j<npointy;j++){
                for(int k=0;k<npointz;k++){
                    
                    int idx=grid.getIndex(i, j, k);

                    point.setX( minX+i*gridResol);
                    point.setY( minY+j*gridResol);
                    point.setZ( minZ+k*gridResol);

                    //if grid point is outside the protein volume
                    if( !volume.isInsideVolume( point, ffield, WATER_RADII) ){
                        grid.setPoint(idx, point, false);
                    }
                    else{
                        grid.setPoint(idx, point, true);
                    }
                }
            }
        }

        return grid;
    }

    /**
     * integrates r-4 inside the VDW volume defined by the protein, excluding
     * the sphere centered in atom (with radius=VDW radius)
     *
     * Method : see NR in C "Monte Carlo integration"
     *
     * @param atom
     * @param grid
     * @param volume  BSP indexing the protein
     * @param ffield
     * @return
     */
    private static double intInsideVolumeR4( Atom atom, AtomGrid grid, BSPTree volume, IForceField ffield){

        double sum=0;

        Atom [] bbox=volume.getBBox();

        //calculates bbox volume
        double boxvolume=(bbox[1].getX()-bbox[0].getX()+2*(MAX_ATOM_RADII+WATER_RADII)) *
                (bbox[1].getY()-bbox[0].getY()+2*(MAX_ATOM_RADII+WATER_RADII)) *
                (bbox[1].getZ()-bbox[0].getZ()+2*(MAX_ATOM_RADII+WATER_RADII));

        double radius2=ffield.getVDWRadius(atom);
        radius2*=radius2;
        

        for(int i=0;i<grid.size();i++){

            if( !grid.isInVolume(i) )
                continue;


            double dist2=CalcGeom.squareDistance( grid.getPoint(i), atom );

            //if grid point inside the atom's sphere
            if( dist2<=radius2 ){
                continue;
            }

            sum+=1.0/(dist2*dist2);
        }

        return (boxvolume*sum)/(grid.size());
    }

    /**
     * calculates born radii with the method of integral inside volume (slower)
     *
     * @param atom
     * @param grid
     * @param volume
     * @param ffield
     * @return
     */
    public static double calcBornRadiiVolInside(Atom atom, AtomGrid grid, BSPTree volume, IForceField ffield){

        double intr4=CalcIntegrals.intInsideVolumeR4(atom, grid, volume, ffield);
        return 1.0/(ffield.getInvVDWRadius(atom)-_1_4_PI*intr4);
    }


    /**
     * calculates born radii with the method of integral outside volume
     *
     * @param atom
     * @param grid
     * @param volume
     * @param accurate : true for accurate radii (slower)
     * @return
     */
    public static double calcBornRadiiVolOutside( Atom atom, AtomGrid grid, BSPTree volume, boolean accurate){

        
        Atom [] bbox=volume.getBBox();

        double addLength=MAX_ATOM_RADII+WATER_RADII;

        int ci=(int) ((atom.getX() - (bbox[0].getX() - addLength)) / grid.getResol());
        int cj=(int) ((atom.getY() - (bbox[0].getY() - addLength)) / grid.getResol());
        int ck=(int) ((atom.getZ() - (bbox[0].getZ() - addLength)) / grid.getResol());

        double [] x=new double [2];
        double [] y=new double [2];
        double [] z=new double [2];

        x[0]=bbox[0].getX() - addLength;
        y[0]=bbox[0].getY() - addLength;
        z[0]=bbox[0].getZ() - addLength;
        x[1]=bbox[1].getX() + addLength;
        y[1]=bbox[1].getY() + addLength;
        z[1]=bbox[1].getZ() + addLength;

        int idx=0;

        for(int i=ci;i>=0;i--){
            idx=grid.getIndex(i, cj, ck);

            if( !grid.isInVolume(idx) ){
                x[0]=grid.getPoint(idx).getX();
                break;
            }
        }

        for(int i=ci;i<grid.getNpointx();i++){
            idx=grid.getIndex(i, cj, ck);

            if( !grid.isInVolume(idx) ){
                x[1]=grid.getPoint(idx).getX();
                break;
            }
        }

        for(int i=cj;i>=0;i--){
            idx=grid.getIndex(ci, i, ck);

            if( !grid.isInVolume(idx) ){
                y[0]=grid.getPoint(idx).getY();
                break;
            }
        }

        for(int i=cj;i<grid.getNpointy();i++){
            idx=grid.getIndex(ci, i, ck);

            if( !grid.isInVolume(idx) ){
                y[1]=grid.getPoint(idx).getY();
                break;
            }
        }

        for(int i=ck;i>=0;i--){
            idx=grid.getIndex(ci, cj, i);

            if( !grid.isInVolume(idx) ){
                z[0]=grid.getPoint(idx).getZ();
                break;
            }
        }

        for(int i=ck;i<grid.getNpointz();i++){
            idx=grid.getIndex(ci, cj, i);

            if( !grid.isInVolume(idx) ){
                z[1]=grid.getPoint(idx).getZ();
                break;
            }
        }

        //center x,y,x
        x[0]-=atom.getX();
        x[1]-=atom.getX();
        y[0]-=atom.getY();
        y[1]-=atom.getY();
        z[0]-=atom.getZ();
        z[1]-=atom.getZ();

        
        double intr4=intOutRectangularR4(x,y,z);

        if( accurate ){
            double intr5=intOutRectangularR5(x,y,z);
            return _4_PI/( -intr4 + 3.0*Math.sqrt(_4_PI*intr5) );
        }

        return _4_PI/intr4;
    }

    
}
