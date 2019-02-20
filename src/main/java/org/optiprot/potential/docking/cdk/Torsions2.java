/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.cdk;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.vecmath.GMatrix;
import javax.vecmath.GVector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.modeling.builder3d.MMFF94ParametersCall;
import org.openscience.cdk.modeling.forcefield.ForceFieldTools;
import org.openscience.cdk.modeling.forcefield.Torsions;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;

/**
 *
 * @author victor
 */
class Torsions2 extends Torsions {

    String functionShape = " Torsions ";

	double mmff94SumET = 0;
	GVector gradientMMFF94SumET = new GVector(3);
	GVector dPhi = new GVector(3);

	GVector order2ndErrorApproximateGradientMMFF94SumET = new GVector(3);
	GVector order5thErrorApproximateGradientMMFF94SumET = new GVector(3);
	GVector xplusSigma = null;
	GVector xminusSigma = null;
	double sigma = Math.pow(0.000000000000001,0.33);

	GMatrix hessianMMFF94SumET = null;
	double[] forHessian = null;
	GMatrix order2ndErrorApproximateHessianMMFF94SumET = null;
	double[] forOrder2ndErrorApproximateHessian = null;

	int torsionNumber = 0;
	int[][] torsionAtomPosition = null;

	double[] v1 = null;
	double[] v2 = null;
	double[] v3 = null;
	double[] phi = null;

	IBond[] bond = null;
	IAtom[] atomInBond = null;
	IBond[] bondConnectedBefore = null;
	IBond[] bondConnectedAfter = null;

	GVector moleculeCurrentCoordinates = null;
	boolean[] changeAtomCoordinates = null;
	int changedCoordinates;

    /**
	 *  Set MMFF94 constants V1, V2 and V3 for each i-j, j-k and k-l bonded pairs in the molecule.
	 *
	 *
	 *@param  molecule       The molecule like an AtomContainer object.
	 *@param  parameterSet   MMFF94 parameters set
	 *@exception  Exception  Description of the Exception
	 */
    @Override
	public void setMMFF94TorsionsParameters(IAtomContainer molecule, Map parameterSet) throws Exception {

        // looks like we need the bonds in an array for the rest of the class
        bond = new IBond[molecule.getBondCount()];
        int counter = 0;
        Iterator bonds = molecule.bonds().iterator();
        while (bonds.hasNext()) {
            IBond aBond = (IBond) bonds.next();
            bond[counter] = aBond;
            counter++;
        }

		for (int b=0; b<bond.length; b++) {
			atomInBond = BondManipulator.getAtomArray(bond[b]);
			bondConnectedBefore = AtomContainerManipulator.getBondArray(molecule.getConnectedBondsList(atomInBond[0]));
			if (bondConnectedBefore.length > 1) {
				bondConnectedAfter = AtomContainerManipulator.getBondArray(molecule.getConnectedBondsList(atomInBond[1]));
				if (bondConnectedAfter.length > 1) {
					for (int bb=0; bb<bondConnectedBefore.length; bb++) {
						if (bondConnectedBefore[bb].compare(bond[b])) {}
						else {
							for (int ba=0; ba<bondConnectedAfter.length; ba++) {
								if (bondConnectedAfter[ba].compare(bond[b])) {}
								else {
									if (bondConnectedBefore[bb].isConnectedTo(bondConnectedAfter[ba])) {}
									else {
										torsionNumber += 1;
									}
								}
							}
						}
					}
				}
			}
		}
		//logger.debug("torsionNumber = " + torsionNumber);

		List torsionsData = null;
		MMFF94ParametersCall pc = new MMFF94ParametersCall();
		pc.initialize(parameterSet);

		v1 = new double[torsionNumber];
		v2 = new double[torsionNumber];
		v3 = new double[torsionNumber];

		torsionAtomPosition = new int[torsionNumber][];

		String torsionType;
		int m = -1;
		for (int b=0; b<bond.length; b++) {
			atomInBond = BondManipulator.getAtomArray(bond[b]);
			bondConnectedBefore = AtomContainerManipulator.getBondArray(molecule.getConnectedBondsList(atomInBond[0]));
			if (bondConnectedBefore.length > 1) {
				bondConnectedAfter = AtomContainerManipulator.getBondArray(molecule.getConnectedBondsList(atomInBond[1]));
				if (bondConnectedAfter.length > 1) {
					for (int bb=0; bb<bondConnectedBefore.length; bb++) {
						if (bondConnectedBefore[bb].compare(bond[b])) {}
						else {
							for (int ba=0; ba<bondConnectedAfter.length; ba++) {
								if (bondConnectedAfter[ba].compare(bond[b])) {}
								else {
									if (bondConnectedBefore[bb].isConnectedTo(bondConnectedAfter[ba])) {}
									else {
										m += 1;
										torsionAtomPosition[m] = new int[4];
										torsionAtomPosition[m][0] = molecule.getAtomNumber(bondConnectedBefore[bb].getConnectedAtom(atomInBond[0]));
										torsionAtomPosition[m][1] = molecule.getAtomNumber(atomInBond[0]);
										torsionAtomPosition[m][2] = molecule.getAtomNumber(atomInBond[1]);
										torsionAtomPosition[m][3] = molecule.getAtomNumber(bondConnectedAfter[ba].getConnectedAtom(atomInBond[1]));

										

										torsionType = "0";
										if (bond[b].getProperty("MMFF94 bond type").toString().equals("1")) {
											torsionType = "1";
										}
										else if ((bond[b].getProperty("MMFF94 bond type").toString().equals("0")) &
												((bondConnectedBefore[bb].getProperty("MMFF94 bond type").toString().equals("1")) |
												(bondConnectedAfter[ba].getProperty("MMFF94 bond type").toString().equals("1")))) {
											torsionType = "2";
										}

										
										torsionsData = pc.getTorsionData(torsionType, bondConnectedBefore[bb].getConnectedAtom(atomInBond[0]).getAtomTypeName(),
												atomInBond[0].getAtomTypeName(), atomInBond[1].getAtomTypeName(), bondConnectedAfter[ba].getConnectedAtom(atomInBond[1]).getAtomTypeName());

										if( torsionsData!=null ){
                                            v1[m] = ((Double) torsionsData.get(0)).doubleValue();
                                            v2[m] = /*(-1) * */((Double) torsionsData.get(1)).doubleValue();
                                            v3[m] = ((Double) torsionsData.get(2)).doubleValue();
                                        }
                                        else{
                                            v1[m] = 0;
                                            v2[m] = 0;
                                            v3[m] = 0;
                                        }

									}
								}
							}
						}
					}
				}
			}
		}

		phi = new double[torsionNumber];

		this.moleculeCurrentCoordinates = new GVector(3 * molecule.getAtomCount());
		for (int i=0; i<moleculeCurrentCoordinates.getSize(); i++) {
			this.moleculeCurrentCoordinates.setElement(i,1E10);
		}

		this.changeAtomCoordinates = new boolean[molecule.getAtomCount()];

	}

    /**
	 *  Evaluate the MMFF94 torsions term.
	 *
	 *@param  coords3d  Current molecule coordinates.
	 *@return        MMFF94 torsions term value.
	 */
    @Override
	public double functionMMFF94SumET(GVector coords3d) {

		setPhi(coords3d);
		mmff94SumET = 0;
		double torsionEnergy=0;
		for (int m = 0; m < torsionNumber; m++) {

            torsionEnergy = v1[m] * (1 + Math.cos(phi[m])) + v2[m] * (1 - Math.cos(2 * phi[m])) + v3[m] * (1 + Math.cos(3 * phi[m]));

			mmff94SumET = mmff94SumET + torsionEnergy;

		}

		return mmff94SumET;
	}

    /**
	 *  Calculate the actual phi
	 *
	 *@param  coords3d  Current molecule coordinates.
	 */
    @Override
	public void setPhi(GVector coords3d) {
		changedCoordinates = 0;
		//logger.debug("Setting Phi");
		for (int i=0; i < changeAtomCoordinates.length; i++) {
			this.changeAtomCoordinates[i] = false;
		}
		this.moleculeCurrentCoordinates.sub(coords3d);
		for (int i = 0; i < this.moleculeCurrentCoordinates.getSize(); i++) {
			//logger.debug("moleculeCurrentCoordinates " + i + " = " + this.moleculeCurrentCoordinates.getElement(i));
			if (Math.abs(this.moleculeCurrentCoordinates.getElement(i)) > 0) {
				changeAtomCoordinates[i/3] = true;
				changedCoordinates = changedCoordinates + 1;
				//logger.debug("changeAtomCoordinates[" + i/3 + "] = " + changeAtomCoordinates[i/3]);
				i = i + (2 - i % 3);
			}
		}

		for (int m = 0; m < torsionNumber; m++) {
			if ((changeAtomCoordinates[torsionAtomPosition[m][0]] == true) |
					(changeAtomCoordinates[torsionAtomPosition[m][1]] == true) |
					(changeAtomCoordinates[torsionAtomPosition[m][2]] == true) |
					(changeAtomCoordinates[torsionAtomPosition[m][3]] == true))		{

				phi[m] = ForceFieldTools.torsionAngleFrom3xNCoordinates(coords3d, torsionAtomPosition[m][0], torsionAtomPosition[m][1],
							torsionAtomPosition[m][2], torsionAtomPosition[m][3]);
			}
		}
		
		moleculeCurrentCoordinates.set(coords3d);

	}
    
}
