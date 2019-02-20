/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.cdk;

import java.util.List;
import java.util.Map;
import javax.vecmath.GMatrix;
import javax.vecmath.GVector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.modeling.builder3d.MMFF94ParametersCall;
import org.openscience.cdk.modeling.forcefield.ForceFieldTools;
import org.openscience.cdk.modeling.forcefield.StretchBendInteractions;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PeriodicTablePositionDescriptor;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author victor
 */
class StretchBendInteractions2 extends StretchBendInteractions {

    String functionShape = " Stretch-Bend Interactions ";

	double mmff94SumEBA = 0;
	GVector gradientMMFF94SumEBA = null;
	GVector order2ndErrorApproximateGradientMMFF94SumEBA = null;
	GVector order5thErrorApproximateGradientMMFF94SumEBA = null;
	GMatrix hessianMMFF94SumEBA = null;
	GVector currentCoordinates = null;
	GVector gradientCurrentCoordinates = null;

	double[][] dDeltarij = null;
	double[][] dDeltarkj = null;
	double[][] dDeltav = null;

	int[][] bondijAtomPosition = null;
	int[][] bondkjAtomPosition = null;
	double[] r0IJ = null;
	double[] r0KJ = null;
	double[] kbaIJK = null;
	double[] kbaKJI = null;
	double[] rij = null;
	double[] rkj = null;
	double[] deltarij = null;
	double[] deltarkj = null;

	BondStretching2 bs = new BondStretching2();
	AngleBending2 ab = new AngleBending2();


	GVector moleculeCurrentCoordinates = null;
	boolean[] changeAtomCoordinates = null;
	int changedCoordinates;

   /**
	 *  Set MMFF94 reference bond lengths r0IJ and r0JK and stretch-bend
	 *  interaction constants kbaIJK and kbaKJI for each i-j-k angle in the
	 *  molecule.
	 *
	 *@param  molecule       The molecule like an AtomContainer object.
	 *@param  parameterSet   MMFF94 parameters set
	 *@exception  Exception  Description of the Exception
	 */
    @Override
	public void setMMFF94StretchBendParameters(IAtomContainer molecule, Map parameterSet, boolean angleBendingFlag) throws Exception {

		
		ab.setMMFF94AngleBendingParameters(molecule, parameterSet, angleBendingFlag);

		IAtom[] atomConnected = null;

		List stretchBendInteractionsData = null;
		List bondData = null;
		MMFF94ParametersCall pc = new MMFF94ParametersCall();
		pc.initialize(parameterSet);

		bondijAtomPosition = new int[ab.angleNumber][];
		bondkjAtomPosition = new int[ab.angleNumber][];
		r0IJ = new double[ab.angleNumber];
		r0KJ = new double[ab.angleNumber];
		kbaIJK = new double[ab.angleNumber];
		kbaKJI = new double[ab.angleNumber];

		String strbndType;
		String angleType;
		String bondIJType;
		String bondKJType;

		IBond bondIJ = null;
		IBond bondKJ = null;

		IAtomicDescriptor descriptor  = new PeriodicTablePositionDescriptor();
		int iR = 0;
		int jR = 0;
		int kR = 0;


		int l = -1;
		for (int j = 0; j < molecule.getAtomCount(); j++) {

			atomConnected = AtomContainerManipulator.getAtomArray(molecule.getConnectedAtomsList(molecule.getAtom(j)));

			if (atomConnected.length > 1) {

				for (int i = 0; i < atomConnected.length; i++) {

					for (int k = i + 1; k < atomConnected.length; k++) {

						l += 1;

						bondIJ = molecule.getBond(atomConnected[i], molecule.getAtom(j));
						bondIJType = bondIJ.getProperty("MMFF94 bond type").toString();

						bondKJ = molecule.getBond(atomConnected[k], molecule.getAtom(j));
						bondKJType = bondKJ.getProperty("MMFF94 bond type").toString();

						angleType = "0";
						if (bondIJType.equals("1") | bondKJType.equals("1")) {
							angleType = "1";
						}
						if (bondIJType.equals("1") & bondKJType.equals("1")) {
							angleType = "2";
						}

						//logger.debug("bondIJType = " + bondIJType + ", bondKJType = " + bondKJType + ", angleType = " + angleType);

						strbndType = "0";
						if (angleType.equals("0") & bondIJType.equals("0") & bondKJType.equals("0")) {strbndType = "0";}
						else if (angleType.equals("1") & bondIJType.equals("1") & bondKJType.equals("0")) {strbndType = "1";}
						else if (angleType.equals("1") & bondIJType.equals("0") & bondKJType.equals("1")) {strbndType = "2";}
						else if (angleType.equals("2") & bondIJType.equals("1") & bondKJType.equals("1")) {strbndType = "3";}
						else if (angleType.equals("4") & bondIJType.equals("0") & bondKJType.equals("0")) {strbndType = "4";}
						else if (angleType.equals("3") & bondIJType.equals("0") & bondKJType.equals("0")) {strbndType = "5";}
						else if (angleType.equals("5") & bondIJType.equals("1") & bondKJType.equals("0")) {strbndType = "6";}
						else if (angleType.equals("5") & bondIJType.equals("0") & bondKJType.equals("1")) {strbndType = "7";}
						else if (angleType.equals("6") & bondIJType.equals("1") & bondKJType.equals("1")) {strbndType = "8";}
						else if (angleType.equals("7") & bondIJType.equals("1") & bondKJType.equals("0")) {strbndType = "9";}
						else if (angleType.equals("7") & bondIJType.equals("0") & bondKJType.equals("1")) {strbndType = "10";}
						else if (angleType.equals("8") & bondIJType.equals("1") & bondKJType.equals("1")) {strbndType = "11";}

						//logger.debug("strbnd: " + strbndType + ", " + atomConnected[i].getAtomTypeName() + "(" + molecule.getAtomNumber(atomConnected[i]) + "), " + molecule.getAtom(j).getAtomTypeName() + "(" + molecule.getAtomNumber(molecule.getAtom(j)) + "), " + ((IAtom)atomConnected.get(k)).getAtomTypeName() + "(" + molecule.getAtomNumber((IAtom)atomConnected.get(k)) + ")");
						stretchBendInteractionsData = pc.getBondAngleInteractionData(strbndType, atomConnected[i].getAtomTypeName(), molecule.getAtom(j).getAtomTypeName(), atomConnected[k].getAtomTypeName());

						if (stretchBendInteractionsData == null) {
							if (angleType.equals("1")) {
								if (strbndType.equals("1")) {strbndType = "2";}
								else {strbndType = "1";}
								//logger.debug("strbnd: " + strbndType + ", " + ((IAtom)atomConnected.get(i)).getAtomTypeName() + "(" + molecule.getAtomNumber((IAtom)atomConnected.get(i)) + "), " + molecule.getAtom(j).getAtomTypeName() + "(" + molecule.getAtomNumber(molecule.getAtom(j)) + "), " + ((IAtom)atomConnected.get(k)).getAtomTypeName() + "(" + molecule.getAtomNumber((IAtom)atomConnected.get(k)) + ")");
								stretchBendInteractionsData = pc.getBondAngleInteractionData(strbndType, atomConnected[i].getAtomTypeName(), molecule.getAtom(j).getAtomTypeName(), atomConnected[k].getAtomTypeName());
							}
                        }

						if (stretchBendInteractionsData == null) {
							iR = ((IntegerResult)descriptor.calculate(atomConnected[i],molecule).getValue()).intValue();
							jR = ((IntegerResult)descriptor.calculate(molecule.getAtom(j),molecule).getValue()).intValue();
							kR = ((IntegerResult)descriptor.calculate(atomConnected[k],molecule).getValue()).intValue();
							stretchBendInteractionsData = pc.getDefaultStretchBendData(iR, jR, kR);

						}

						if (stretchBendInteractionsData != null) {
                            kbaIJK[l] = ((Double) stretchBendInteractionsData.get(0)).doubleValue();
                            kbaKJI[l] = ((Double) stretchBendInteractionsData.get(1)).doubleValue();
                        }
                        else{
                            kbaIJK[l] = 0;
                            kbaKJI[l] = 0;
                        }
                        

						
                        r0IJ[l] = 0;
                        r0KJ[l] = 0;
						bondData = pc.getBondData(bondIJType, atomConnected[i].getAtomTypeName(), molecule.getAtom(j).getAtomTypeName());

                        if( bondData!=null )
                            r0IJ[l] = ((Double) bondData.get(0)).doubleValue();

                        bondData = pc.getBondData(bondKJType, atomConnected[k].getAtomTypeName(), molecule.getAtom(j).getAtomTypeName());
						
                        if( bondData!=null )
                            r0KJ[l] = ((Double) bondData.get(0)).doubleValue();

						bondijAtomPosition[l] = new int[2];
						bondijAtomPosition[l][0] = molecule.getAtomNumber(atomConnected[i]);
						bondijAtomPosition[l][1] = j;

						bondkjAtomPosition[l] = new int[2];
						bondkjAtomPosition[l][0] = molecule.getAtomNumber(atomConnected[k]);
						bondkjAtomPosition[l][1] = j;
					}
				}
			}
		}
		rij = new double[ab.angleNumber];
		rkj = new double[ab.angleNumber];
		deltarij = new double[ab.angleNumber];
		deltarkj = new double[ab.angleNumber];
		currentCoordinates = new GVector(3 * molecule.getAtomCount());
		gradientCurrentCoordinates = new GVector(3 * molecule.getAtomCount());
		gradientMMFF94SumEBA = new GVector(3 * molecule.getAtomCount());
		dDeltarij = new double[3 * molecule.getAtomCount()][];
		dDeltarkj = new double[3 * molecule.getAtomCount()][];
		dDeltav = new double[3 * molecule.getAtomCount()][];
		hessianMMFF94SumEBA = new GMatrix(3 * molecule.getAtomCount(), 3 * molecule.getAtomCount());
		for (int i = 0; i < 3 * molecule.getAtomCount(); i++) {
			dDeltarij[i] = new double[ab.angleNumber];
			dDeltarkj[i] = new double[ab.angleNumber];
			dDeltav[i] = new double[ab.angleNumber];
		}

		this.moleculeCurrentCoordinates = new GVector(3 * molecule.getAtomCount());
		for (int i=0; i<moleculeCurrentCoordinates.getSize(); i++) {
			this.moleculeCurrentCoordinates.setElement(i,1E10);
		}

		this.changeAtomCoordinates = new boolean[molecule.getAtomCount()];

	}

    /**
	 *  Set the MMFF94 stretch-bend interaction term given the atoms cartesian
	 *  coordinates.
	 *
	 *@param  coords3d  Current molecule coordinates.
	 */
    @Override
	public void setFunctionMMFF94SumEBA(GVector coords3d) {
		
		if (currentCoordinates.equals(coords3d)) {
		}
        else {
			setDeltarijAndDeltarkj(coords3d);
			ab.setDeltav(coords3d);
			mmff94SumEBA = 0;
			for (int j = 0; j < ab.angleNumber; j++) {

				mmff94SumEBA = mmff94SumEBA + 2.51210 * (kbaIJK[j] * deltarij[j] + kbaKJI[j] * deltarkj[j]) * ab.deltav[j];

			}
			
			currentCoordinates.set(coords3d);
		}
	}

    /**
	 *  Calculate the current bond distances rij and rkj for each angle j, and the
	 *  difference with the reference bonds.
	 *
	 *@param  coords3d  Current molecule coordinates.
	 */
    @Override
	public void setDeltarijAndDeltarkj(GVector coords3d) {

		changedCoordinates = 0;
		//logger.debug("Setting Deltarij and Deltarkj");
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

		for (int i = 0; i < ab.angleNumber; i++) {
			if ((changeAtomCoordinates[ab.angleAtomPosition[i][0]] == true) |
					(changeAtomCoordinates[ab.angleAtomPosition[i][1]] == true))		{

				rij[i] = ForceFieldTools.distanceBetweenTwoAtomsFrom3xNCoordinates(coords3d, ab.angleAtomPosition[i][1], ab.angleAtomPosition[i][0]);
				deltarij[i] = rij[i] - r0IJ[i];
				//logger.debug("deltarij[" + i + "] = " + deltarij[i]);
			}
			//else {System.out.println("deltarij[" + i + "] was no recalculated");}
			if ((changeAtomCoordinates[ab.angleAtomPosition[i][1]] == true) |
					(changeAtomCoordinates[ab.angleAtomPosition[i][2]] == true))		{

				rkj[i] = ForceFieldTools.distanceBetweenTwoAtomsFrom3xNCoordinates(coords3d, ab.angleAtomPosition[i][1], ab.angleAtomPosition[i][2]);
				deltarkj[i] = rkj[i] - r0KJ[i];
				//logger.debug("deltarkj[" + i + "] = " + deltarkj[i]);
			}
		}

		moleculeCurrentCoordinates.set(coords3d);
	}
    
}
