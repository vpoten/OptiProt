/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.optiprot.cgal.*;
import org.optiprot.potential.IForceField;

/**
 *
 * @author victor
 */
public class CalcSurfVol {

    static final double RADII_WATER=1.4;//radii of a water molecule

    /**
     * calculates the solvent accesible surface of a molecule
     *
     * @param molecule
     * @param ff
     * @param maxDepth control recursive calc of 3 spheres intersect. (use 3 or 4)
     * @return
     */
    public static double SAsurface( Atom [] molecule, IForceField ff, int maxDepth)
            throws StructureException{

        return surface(molecule, ff, RADII_WATER, maxDepth);
    }

    public static double VDWsurface( Atom [] molecule, IForceField ff, int maxDepth)
            throws StructureException{

        return surface(molecule, ff, 0, maxDepth);
    }

    
    private static double surface( Atom [] molecule, IForceField ff, double probe,
            int maxDepth )
            throws StructureException{

        if( molecule.length==1 ){
            return CalcSphere.area( ff.getVDWRadius(molecule[0])+probe );
        }
        else if( molecule.length==2 ){
            return CalcSphere.intersecArea( molecule[0],
                    ff.getVDWRadius(molecule[0])+probe,
                    molecule[1],
                    ff.getVDWRadius(molecule[1])+probe );
        }
        else if( molecule.length==3 ){
            return CalcSphere.intersecArea( molecule[0],
                    ff.getVDWRadius(molecule[0])+probe,
                    molecule[1],
                    ff.getVDWRadius(molecule[1])+probe,
                    molecule[2],
                    ff.getVDWRadius(molecule[2])+probe, maxDepth );
        }
       
        double [] radii=CalcSurfVol.createArrayRadius(molecule,ff);

        AlphaShape as=new AlphaShape(molecule,radii,probe);
       
        //short inclusion/exclusion method, using the outside-fringe
        double area1=0;
        double area2=0;
        double area3=0;

        Atom [] centers=null;
        double [] sph_radii=null;
        double out_angle=0;

        for(BaseSimplex simplex : as.getSimplices() ){

            centers=simplex.getVertices();
            sph_radii=simplex.getWeights();
            out_angle=simplex.getOutAngle();

            //if( out_angle<0.0 ){
            //    throw new StructureException("Bad out angle");
            //}
            
            if( simplex instanceof _0Simplex ){

                area1 += out_angle * CalcSphere.area(sph_radii[0]);
            }
            else if( simplex instanceof _1Simplex ){

                area2 -= out_angle * CalcSphere.intersecArea(
                        centers[0],sph_radii[0],
                        centers[1],sph_radii[1]);
            }
            else if( simplex instanceof _2Simplex ){
                
                try {
                    area3 += out_angle * CalcSphere.intersecArea(
                            centers[0], sph_radii[0],
                            centers[1], sph_radii[1],
                            centers[2], sph_radii[2], maxDepth );
                } catch (StructureException ex) {
                    Logger.getLogger(CalcSurfVol.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

        }

        return area1+area2+area3;
    }


    private static double[] createArrayRadius(Atom[] molecule, IForceField ff) {

        double [] radii=new double [molecule.length];

        int i=0;

        for( Atom atom : molecule){
            radii[i++]=ff.getVDWRadius(atom);
        }

        return radii;
    }
}
