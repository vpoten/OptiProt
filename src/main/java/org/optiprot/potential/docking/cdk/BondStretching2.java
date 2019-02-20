/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.cdk;

/**
 *
 * @author victor
 */
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.vecmath.GMatrix;
import javax.vecmath.GVector;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.modeling.builder3d.MMFF94ParametersCall;
import org.openscience.cdk.modeling.forcefield.BondStretching;
import org.openscience.cdk.modeling.forcefield.ForceFieldTools;


/**
 *  Bond Stretching calculator for the potential energy function.
 * Include function and derivatives.
 *
 * @author     victor
 */
public class BondStretching2 extends BondStretching {

	String functionShape = " Bond Stretching ";


	double mmff94SumEB = 0;
	GVector gradientMMFF94SumEB = null;
	GMatrix hessianMMFF94SumEB = null;
	double[] forHessian = null;
	GMatrix order2ndErrorApproximateHessianMMFF94SumEB = null;
	double[] forOrder2ndErrorApproximateHessian = null;

	int bondsNumber;
	int[][] bondAtomPosition = null;

	double[] r0 = null;	// Force field parameters
	double[] k2 = null;
	double[] k3 = null;
	double[] k4 = null;
	double cs = -2;

	double[] r = null;	// The actual bond lengths
	double[] deltar = null;	// The difference between actual and reference bond lengths

	double[][] dDeltar = null;
	double[][][] ddDeltar = null;

	


	/**
	 *  Set MMFF94 reference bond lengths r0IJ and the constants k2, k3, k4 for
	 *  each i-j bond in the molecule.
	 *
	 *@param  molecule       The molecule like an AtomContainer object.
	 *@param  parameterSet   MMFF94 parameters set
	 *@exception  Exception  Description of the Exception
	 */
    @Override
	public void setMMFF94BondStretchingParameters(IAtomContainer molecule, Map parameterSet) throws Exception {

		//logger.debug("setMMFF94BondStretchingParameters");

		bondsNumber = molecule.getBondCount();
		//logger.debug("bondsNumber = " + bondsNumber);
		bondAtomPosition = new int[molecule.getBondCount()][];

		List bondData = null;
		MMFF94ParametersCall pc = new MMFF94ParametersCall();
		pc.initialize(parameterSet);

		r0 = new double[molecule.getBondCount()];
		k2 = new double[molecule.getBondCount()];
		k3 = new double[molecule.getBondCount()];
		k4 = new double[molecule.getBondCount()];

        String bondType;
        Iterator bonds = molecule.bonds().iterator();
        int i = 0;
        while (bonds.hasNext()) {
            IBond bond = (IBond) bonds.next();

            //atomsInBond = bonds[i].getatoms();

            bondType = bond.getProperty("MMFF94 bond type").toString();
            //logger.debug("bondType " + i + " = " + bondType);

            bondAtomPosition[i] = new int[bond.getAtomCount()];

            for (int j = 0; j < bond.getAtomCount(); j++) {
                bondAtomPosition[i][j] = molecule.getAtomNumber(bond.getAtom(j));
            }


            bondData = pc.getBondData(bondType, bond.getAtom(0).getAtomTypeName(), bond.getAtom(1).getAtomTypeName());

            if( bondData!=null ){
                r0[i] = ((Double) bondData.get(0)).doubleValue();
                k2[i] = ((Double) bondData.get(1)).doubleValue();
                k3[i] = ((Double) bondData.get(2)).doubleValue();
                k4[i] = ((Double) bondData.get(3)).doubleValue();
            }
            else{
                r0[i] = 0;
                k2[i] = 0;
                k3[i] = 0;
                k4[i] = 0;
            }

            i++;
        }

        r = new double[molecule.getBondCount()];
		deltar = new double[molecule.getBondCount()];


	}


    /**
	 *  Evaluate the MMFF94 bond stretching term for the given atoms cartesian coordinates.
	 *
	 *@param  coord3d  Current molecule coordinates.
	 *@return        bond stretching value
	 */
    @Override
	public double functionMMFF94SumEB(GVector coord3d) {

		calculateDeltar(coord3d);


		mmff94SumEB = 0;

		for (int i = 0; i < bondsNumber; i++) {
			mmff94SumEB = mmff94SumEB + k2[i] * Math.pow(deltar[i],2)
							+ k3[i] * Math.pow(deltar[i],3) + k4[i] * Math.pow(deltar[i],4);
		}

		
		return mmff94SumEB;
	}

    /**
	 *  Calculate the actual bond distance rij and the difference with the reference bond distances.
	 *
	 *@param  coord3d  Current molecule coordinates.
	 */
    @Override
	public void calculateDeltar(GVector coord3d) {

		for (int i = 0; i < bondAtomPosition.length; i++) {
			r[i] = ForceFieldTools.distanceBetweenTwoAtomsFrom3xNCoordinates(coord3d, bondAtomPosition[i][0], bondAtomPosition[i][1]);
			deltar[i] = r[i] - r0[i];
		}
	}

}
