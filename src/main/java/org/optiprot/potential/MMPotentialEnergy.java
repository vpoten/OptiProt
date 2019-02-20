/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

import org.biojava.bio.structure.Chain;
import org.optiprot.OptiProtParameters;
import org.optiprot.potential.element.*;

/**
 * Molecular mechanics potential energy
 * @author victor
 */
public class MMPotentialEnergy {

   
    private static double KCOULOMB=332;
    private final static double KGBORN=332;

    private static double istrength=0;// ionic strength in Mol/l

    private static double idebye=0;//  inverse Debye-Hückel length

    private static double iEin=0;// inverse protein and ligand dielectric constant

    static {
        //calc. constants for GBorn solvation energy
        IForceField ffield=new TestForceField();

        istrength=0.1;//ionic strength
        idebye=Math.sqrt(istrength/ffield.getKWaterDielectric())/0.343;
        iEin=1.0/ffield.getKProtLigDielectric();
    }

    private static double VDW_CAP=2000;//Kcal/mol/pair
    private static double ELECTROST_CAP=1000;//Kcal/mol/pair

    private static double DIST_MAX_VDW = 9.0;//in angstroms
    
    
    static public Double calcEnergy( IMolecularElements molElements,
            OptiProtParameters parameters ){

        double energy=0;
        IForceField ffield=parameters.getForceField();

        energy+=calcBondLengthEnergy( molElements, ffield );

        energy+=calcAngleEnergy( molElements, ffield );

        energy+=calcTorsionEnergy( molElements, ffield );

        energy+=calcCoulombVDWEnergy( molElements, ffield );
       
        energy+=calcGBornSolvationEnergy( molElements, ffield );

        energy+=calcSASA( molElements.getChain(), ffield, parameters );

        return energy;
    }


    static private Double calcBondLengthEnergy( IMolecularElements molElements,
            IForceField ffield ){

        double energy=0.0;
        
        if( !molElements.hasBonds() ){
            return 0.0;
        }

        for( MolecularBond bond : molElements.getBonds() ){
            energy+=calcBondLengthEnergy( bond, ffield);
        }
        

        return energy;
    }


    static public Double calcBondLengthEnergy( MolecularBond bond,
            IForceField ffield ){

        //double kb, double b0
        double [] params=ffield.getParBond(bond);

        if(params==null)
            return 0.0;

        double temp=bond.getDistance()-params[1];
        return params[0]*temp*temp;
    }



    static private Double calcAngleEnergy( IMolecularElements molElements,
            IForceField ffield ){

        double energy=0.0;
       
        if( !molElements.hasAngles() ){
            return 0.0;
        }

        for( MolecularAngle angle : molElements.getAngles() ){
            energy+=calcAngleEnergy( angle, ffield);
        }
        

        return energy;
    }



    static public Double calcAngleEnergy( MolecularAngle angle,
            IForceField ffield ){

        double energy=0.0;
        double temp=0.0;

        //double ktheta, double theta0, double kub, double s0
        double [] params=ffield.getParAngle(angle);

        if(params==null)
            return 0.0;

        temp=angle.getDistance()-params[3];
        energy+=params[2]*temp*temp;

        temp=angle.getAngle()-params[1];
        energy+=params[0]*temp*temp;

        return energy;
    }



    static private Double calcTorsionEnergy( IMolecularElements molElements,
            IForceField ffield ){

        double energy=0.0;
        if( molElements.hasDihedrals() ){

            for( MolecularDihedral dihedral : molElements.getDihedrals() ){
                energy+=calcTorsionEnergy( dihedral, ffield);
            }
        }

        if( molElements.hasImpropers() ){

            for( MolecularImproper improper : molElements.getImpropers() ){
                energy+=calcTorsionEnergy( improper, ffield);
            }
        }

        return energy;
    }



     static public Double calcTorsionEnergy( MolecularDihedral dihedral,
            IForceField ffield ){

         //double kchi, double n, double delta
        double [] params=ffield.getParDihedral(dihedral);

        if(params==null)
            return 0.0;

         double temp=Math.toRadians(params[1]*dihedral.getAngle() - params[2]);
         return params[0]*(1+Math.cos(temp));
     }



     static public Double calcTorsionEnergy( MolecularImproper improper,
            IForceField ffield ){

         //double kpsi, double psi0
         double [] params=ffield.getParImproper(improper);

        if(params==null)
            return 0.0;

         double temp=improper.getAngle()-params[1];
         return params[0]*temp*temp;
     }



    static private Double calcCoulombVDWEnergy( IMolecularElements molElements,
            IForceField ffield ){

        double energy=0.0;
       
       if( !molElements.hasNonbondeds() ){
            return 0.0;
        }

        for( MolecularNonbonded nonbonded : molElements.getNonbondeds() ){
            energy+=calcCoulombVDWEnergy( nonbonded, ffield);
        }
       
        return energy;
    }



    static public Double calcCoulombVDWEnergy(  MolecularNonbonded nonbonded,
            IForceField ffield ){

        double temp=0.0;
        double coulomb=0.0;
        double lennardj=0.0;

        if( nonbonded.getDistance()<DIST_MAX_VDW ){
            
            //double eps, double rmin2, double eps14, double rmin2_14
            double [] params=ffield.getParVDW(nonbonded);

            if(params!=null){
                temp = params[1]*params[1]/nonbonded.getSqrDistance();
                temp = temp*temp*temp;
                lennardj = params[0]*
                        ( temp*temp-2.0*temp );
            }

            //capping VdW
            if( lennardj>VDW_CAP ){
                lennardj=VDW_CAP;
            }
        }

        coulomb = ffield.getChargeAB(nonbonded) /
                (ffield.getKProtLigDielectric()*nonbonded.getDistance());

        //capping coulomb
        if( coulomb>ELECTROST_CAP ){
            coulomb=ELECTROST_CAP;
        }
        else if( coulomb<-ELECTROST_CAP ){
            coulomb=-ELECTROST_CAP;
        }

        return lennardj+KCOULOMB*coulomb;
    }



    static private Double calcGBornSolvationEnergy( IMolecularElements molElements,
            IForceField ffield ){

        double energy=0.0;

        if( !molElements.hasPairs() ){
            return 0.0;
        }

        for( MolecularPair pair : molElements.getPairs() ){
            energy+=calcGBornSolvationEnergy( pair, ffield );
        }

        return energy;
    }



    static public Double calcGBornSolvationEnergy( MolecularPair pair,
            IForceField ffield ){


        double radix=Math.sqrt( pair.getSqrDistance() + pair.getBornRadiusAB()*
                    Math.exp( -pair.getSqrDistance()/(4.0*pair.getBornRadiusAB())) );

        double gborn = (iEin-(Math.exp(-idebye*radix)/ffield.getKWaterDielectric())) *
                   (ffield.getChargeAB(pair)/radix) ;

        //capping gborn
        if( gborn>ELECTROST_CAP ){
            gborn=ELECTROST_CAP;
        }
        else if( gborn<-ELECTROST_CAP ){
            gborn=-ELECTROST_CAP;
        }

        return -KGBORN*gborn;
    }


    
    static public Double calcSASA( Chain chain, IForceField ffield,
            OptiProtParameters parameters ){

        double area=NVSASACalc.calc( chain, parameters.getRotLib() );
        return ffield.getKSurfTensionWater()*area;
    }

    
}
