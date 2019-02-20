/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential;

/**
 * interface for force fields tables of parameters
 * @author victor
 */
public interface IForceFieldTable {

////    public void addBond( String at1, String at2, double kb, double b0 );
////
////    public void addAngle( String at1, String at2, String at3, double ktheta, double theta0,
////            double kub, double s0);
////
////    public void addDihedral( String at1, String at2, String at3, String at4,
////            double kchi, double n, double delta);
////
////    public void addImproper( String at1, String at2, String at3, String at4,
////            double kpsi, double psi0);
////
////    public void addNonbonded( String at1, double eps, double rmin2, double eps14,
////            double rmin2_14 );
////
////    public double [] getBond( String at1, String at2 );
////
////    public double [] getAngle( String at1, String at2, String at3 );
////
////    public double [] getDihedral( String at1, String at2, String at3, String at4 );
////
////    public double [] getImproper( String at1, String at2, String at3, String at4 );
////
////    public double [] getNonbonded( String at1 );

    public void addBond( Integer at1, Integer at2, double kb, double b0 );

    public void addAngle( Integer at1, Integer at2, Integer at3, double ktheta, double theta0,
            double kub, double s0);

    public void addDihedral( Integer at1, Integer at2, Integer at3, Integer at4,
            double kchi, double n, double delta);

    public void addImproper( Integer at1, Integer at2, Integer at3, Integer at4,
            double kpsi, double psi0);

    public void addNonbonded( Integer at1, double eps, double rmin2, double eps14,
            double rmin2_14 );

    public double [] getBond( Integer at1, Integer at2 );

    public double [] getAngle( Integer at1, Integer at2, Integer at3 );

    public double [] getDihedral( Integer at1, Integer at2, Integer at3, Integer at4 );

    public double [] getImproper( Integer at1, Integer at2, Integer at3, Integer at4 );

    public double [] getNonbonded( Integer at1 );

   
}
