/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking;

import Jama.Matrix;
import org.biojava.bio.structure.Atom;
import org.optiprot.maths.CalcGeom;
import org.optiprot.potential.docking.element.DockingAtomClassify.simpletype;

/**
 * a simple function for docking score
 *
 * @author victor
 */
public class SimpleDockScore implements IDockingScore  {

    private int nevaluations=0;

    // types of interactions
    static private int num_inter = 5;
    static private int int_unk = -1; //unknown
    static private int int_l_l = 0; //hydrophobic - hydrophobic
    static private int int_l_h = 1; //hydrophobic - hydrophilic
    static private int int_h_h = 2; //hydrogen bond
    static private int int_me = 3; //ligand - metal ion
    static private int int_re = 4; //repulsive

    //parameters for interactions

    static private double [] opt_distances = new double [num_inter];
    static private double [] max_distances = new double [num_inter];

    static private double [] eps_LJ = new double [num_inter];//lennard jones
    static private double [] energ_base = new double [num_inter];

    static private double [][] polin_params = new double [num_inter][4];

    static private double kRot = 0.33;//entropy loss (rotating a bond) (kcal/mol)

    //constants for overlap penalties
    static private double F_PENAL1 = 0;
    static private double F_PENAL2 = 0;

    static private int[][] interTable=null;

    //initialize parameters
    static{

        calcPenalFactor( 1000, 7500);

        interTable=createInterTable();

        //in kcal/mol
        eps_LJ[int_l_l] = 0.06;
        eps_LJ[int_l_h] = 0.3;
        eps_LJ[int_h_h] = 0.6;
        eps_LJ[int_me] = 0.6;
        eps_LJ[int_re] = 0.06;

        //in kcal/mol
        energ_base[int_l_l] = -0.1;
        energ_base[int_l_h] = -0.05;
        energ_base[int_h_h] = -1.5;
        energ_base[int_me] = -1.5;
        energ_base[int_re] = 0.10;

        //in A
        opt_distances[int_l_l] = 4.1;
        opt_distances[int_l_h] = 3.6;
        opt_distances[int_h_h] = 2.8;
        opt_distances[int_me] = 2.0;
        opt_distances[int_re] = 3.2;

        max_distances[int_l_l] = 6.5;
        max_distances[int_l_h] = 5.5;
        max_distances[int_h_h] = 3.5;
        max_distances[int_me] = 3.5;
        max_distances[int_re] = 5.0;

        Matrix params=null;

        params = calcParameters(opt_distances[int_l_l],max_distances[int_l_l],energ_base[int_l_l]);
        polin_params[int_l_l][0]=params.get(0,0);
        polin_params[int_l_l][1]=params.get(1,0);
        polin_params[int_l_l][2]=params.get(2,0);
        polin_params[int_l_l][3]=params.get(3,0);

        params = calcParameters(opt_distances[int_l_h],max_distances[int_l_h],energ_base[int_l_h]);
        polin_params[int_l_h][0]=params.get(0,0);
        polin_params[int_l_h][1]=params.get(1,0);
        polin_params[int_l_h][2]=params.get(2,0);
        polin_params[int_l_h][3]=params.get(3,0);

        params = calcParameters(opt_distances[int_h_h],max_distances[int_h_h],energ_base[int_h_h]);
        polin_params[int_h_h][0]=params.get(0,0);
        polin_params[int_h_h][1]=params.get(1,0);
        polin_params[int_h_h][2]=params.get(2,0);
        polin_params[int_h_h][3]=params.get(3,0);

        params = calcParameters(opt_distances[int_me],max_distances[int_me],energ_base[int_me]);
        polin_params[int_me][0]=params.get(0,0);
        polin_params[int_me][1]=params.get(1,0);
        polin_params[int_me][2]=params.get(2,0);
        polin_params[int_me][3]=params.get(3,0);

        params = calcParameters(opt_distances[int_re],max_distances[int_re],energ_base[int_re]);
        polin_params[int_re][0]=params.get(0,0);
        polin_params[int_re][1]=params.get(1,0);
        polin_params[int_re][2]=params.get(2,0);
        polin_params[int_re][3]=params.get(3,0);

    }





    public double calcScore(Docked dock) {

        double energy=0;
        nevaluations++;

        for( Atom atL : dock.getLigandAtoms() ){
            for( Atom atP : dock.getProteinAtoms() ){

                //find out interaction type
                simpletype t1 =
                        dock.getLigand().getSimpletype( atL );

                simpletype t2 =
                        dock.getProtein().getSimpletype( atP );
               
                int interact_type = getInteracType(t1, t2);
               
                double sqdist = CalcGeom.squareDistance(atL, atP);

                //choose formula according distance and interaction
                energy += calcEnergy( interact_type, sqdist );

            }
        }

        //internal energy
        double energy2 = 0;
        
        if( dock.getLigand().getInternalEnergy()!=null )
            energy2 = dock.getLigand().getInternalEnergy();

        //entropy loss
        double energy3 = dock.getLigand().getNumRotableBonds()*kRot;

        double penalty = calcPenalty( dock.calcOverlapFactor() );
        
        energy += energy2 + energy3 + penalty;

        return energy;
    }

    /**
     * 
     * @param t1 : type of ligand atom
     * @param t2 : type of protein atom
     * @return
     */
    private int getInteracType( simpletype t1, simpletype t2){

        return interTable[getIndex(t1)][getIndex(t2)];
    }

    /**
     *
     * @param overlap
     * @return
     */
    private double calcPenalty( double overlap ){

        double val=1.0-overlap;
        return F_PENAL1*val*val+F_PENAL2*val;
    }

    /**
     * 
     * @param interact_type
     * @param sqdist
     * @return
     */
    private double calcEnergy( int interact_type, double sqdist ){

        if( interact_type==int_unk )
            return 0.0;

        double r2=max_distances[interact_type];

        if( sqdist>r2*r2 )
            return 0.0;

        double e=energ_base[interact_type];
        double eps=eps_LJ[interact_type];
        double r1=opt_distances[interact_type];
        double sqr1=r1*r1;

        if( sqdist<sqr1 ){
            double temp = sqr1/sqdist;
            temp = temp*temp*temp;
            return (e + eps*( temp*temp-2.0*temp+1.0 ));
        }

        double dist=Math.sqrt(sqdist);

        return polin_params[interact_type][0]*dist*sqdist +
                polin_params[interact_type][1]*sqdist +
                polin_params[interact_type][2]*dist +
                polin_params[interact_type][3];
    }

    /**
     * calc paramaters a,b,c,d for ar³+br²+cr+d
     * 
     * with f(r1)=fr1, f(r2)=0, f'(r1)=0, f'(r2)=0
     * 
     * @param r1
     * @param r2
     * @param fr1
     * @return
     */
    static private Matrix calcParameters( double r1, double r2, double fr1 ){

        Matrix mat=new Matrix(4,4);

        mat.set(0, 0, r1*r1*r1);
        mat.set(0, 1, r1*r1);
        mat.set(0, 2, r1);
        mat.set(0, 3, 1.0);

        mat.set(1, 0, r2*r2*r2);
        mat.set(1, 1, r2*r2);
        mat.set(1, 2, r2);
        mat.set(1, 3, 1.0);

        mat.set(2, 0, 3.0*r1*r1);
        mat.set(2, 1, 2.0*r1);
        mat.set(2, 2, 1.0);
        mat.set(2, 3, 0.0);

        mat.set(3, 0, 3.0*r2*r2);
        mat.set(3, 1, 2.0*r2);
        mat.set(3, 2, 1.0);
        mat.set(3, 3, 0.0);


        Matrix b=new Matrix(4,1);
        b.set(0, 0, fr1);

        return mat.solve(b);
    }

    /**
     * 
     * @param pen1 : penalty at 0.1
     * @param pen2 : penalty at 0.2
     */
    static private void calcPenalFactor( double pen1, double pen2){
        
        F_PENAL1=pen1/0.01;
        F_PENAL2=0.0;
    }

    /**
     * @return the nevaluations
     */
    public int getNevaluations() {
        return nevaluations;
    }

    /**
     *
     * @param type
     * @return
     */
    static private int getIndex( simpletype type ){

        switch( type ){
            case donorh: return 0;
            case donor: return 1;
            case donor_accept: return 2;
            case acceptor: return 3;
            case metal: return 4;
            case nonpolar: return 5;
            case nonpolarh: return 6;
        }

        return -1;
    }

    /**
     * creates the table of interactions, the row are the ligand types and
     * the column are the protein types
     *
     * @return
     */
    private static int[][] createInterTable() {

        int [][] table=new int[simpletype.values().length][simpletype.values().length];

        table[getIndex(simpletype.donorh)][getIndex(simpletype.donorh)]=int_re;
        table[getIndex(simpletype.donorh)][getIndex(simpletype.donor)]=int_re;
        table[getIndex(simpletype.donorh)][getIndex(simpletype.donor_accept)]=int_h_h;
        table[getIndex(simpletype.donorh)][getIndex(simpletype.acceptor)]=int_h_h;
        table[getIndex(simpletype.donorh)][getIndex(simpletype.metal)]=int_re;
        table[getIndex(simpletype.donorh)][getIndex(simpletype.nonpolar)]=int_l_h;
        table[getIndex(simpletype.donorh)][getIndex(simpletype.nonpolarh)]=int_re;

        table[getIndex(simpletype.metal)][getIndex(simpletype.donorh)]=int_unk;
        table[getIndex(simpletype.metal)][getIndex(simpletype.donor)]=int_unk;
        table[getIndex(simpletype.metal)][getIndex(simpletype.donor_accept)]=int_unk;
        table[getIndex(simpletype.metal)][getIndex(simpletype.acceptor)]=int_unk;
        table[getIndex(simpletype.metal)][getIndex(simpletype.metal)]=int_unk;
        table[getIndex(simpletype.metal)][getIndex(simpletype.nonpolar)]=int_unk;
        table[getIndex(simpletype.metal)][getIndex(simpletype.nonpolarh)]=int_unk;

        table[getIndex(simpletype.donor)][getIndex(simpletype.donorh)]=int_re;
        table[getIndex(simpletype.donor)][getIndex(simpletype.donor)]=int_re;
        table[getIndex(simpletype.donor)][getIndex(simpletype.donor_accept)]=int_l_h;
        table[getIndex(simpletype.donor)][getIndex(simpletype.acceptor)]=int_l_h;
        table[getIndex(simpletype.donor)][getIndex(simpletype.metal)]=int_re;
        table[getIndex(simpletype.donor)][getIndex(simpletype.nonpolar)]=int_l_h;
        table[getIndex(simpletype.donor)][getIndex(simpletype.nonpolarh)]=int_l_h;

        table[getIndex(simpletype.donor_accept)][getIndex(simpletype.donorh)]=int_h_h;
        table[getIndex(simpletype.donor_accept)][getIndex(simpletype.donor)]=int_l_h;
        table[getIndex(simpletype.donor_accept)][getIndex(simpletype.donor_accept)]=int_l_h;
        table[getIndex(simpletype.donor_accept)][getIndex(simpletype.acceptor)]=int_l_h;
        table[getIndex(simpletype.donor_accept)][getIndex(simpletype.metal)]=int_me;
        table[getIndex(simpletype.donor_accept)][getIndex(simpletype.nonpolar)]=int_l_h;
        table[getIndex(simpletype.donor_accept)][getIndex(simpletype.nonpolarh)]=int_l_h;

        table[getIndex(simpletype.acceptor)][getIndex(simpletype.donorh)]=int_h_h;
        table[getIndex(simpletype.acceptor)][getIndex(simpletype.donor)]=int_l_h;
        table[getIndex(simpletype.acceptor)][getIndex(simpletype.donor_accept)]=int_l_h;
        table[getIndex(simpletype.acceptor)][getIndex(simpletype.acceptor)]=int_re;
        table[getIndex(simpletype.acceptor)][getIndex(simpletype.metal)]=int_me;
        table[getIndex(simpletype.acceptor)][getIndex(simpletype.nonpolar)]=int_l_h;
        table[getIndex(simpletype.acceptor)][getIndex(simpletype.nonpolarh)]=int_l_h;

        table[getIndex(simpletype.nonpolar)][getIndex(simpletype.donorh)]=int_l_h;
        table[getIndex(simpletype.nonpolar)][getIndex(simpletype.donor)]=int_l_h;
        table[getIndex(simpletype.nonpolar)][getIndex(simpletype.donor_accept)]=int_l_h;
        table[getIndex(simpletype.nonpolar)][getIndex(simpletype.acceptor)]=int_l_h;
        table[getIndex(simpletype.nonpolar)][getIndex(simpletype.metal)]=int_l_h;
        table[getIndex(simpletype.nonpolar)][getIndex(simpletype.nonpolar)]=int_l_l;
        table[getIndex(simpletype.nonpolar)][getIndex(simpletype.nonpolarh)]=int_l_h;

        table[getIndex(simpletype.nonpolarh)][getIndex(simpletype.donorh)]=int_re;
        table[getIndex(simpletype.nonpolarh)][getIndex(simpletype.donor)]=int_l_h;
        table[getIndex(simpletype.nonpolarh)][getIndex(simpletype.donor_accept)]=int_l_h;
        table[getIndex(simpletype.nonpolarh)][getIndex(simpletype.acceptor)]=int_l_h;
        table[getIndex(simpletype.nonpolarh)][getIndex(simpletype.metal)]=int_l_h;
        table[getIndex(simpletype.nonpolarh)][getIndex(simpletype.nonpolar)]=int_l_h;
        table[getIndex(simpletype.nonpolarh)][getIndex(simpletype.nonpolarh)]=int_l_h;

        return table;
    }
    
}
