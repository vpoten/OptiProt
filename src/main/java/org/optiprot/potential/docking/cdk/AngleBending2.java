/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking.cdk;

import java.util.List;
import java.util.Map;
import javax.vecmath.GVector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.modeling.builder3d.MMFF94ParametersCall;
import org.openscience.cdk.modeling.forcefield.AngleBending;
import org.openscience.cdk.modeling.forcefield.ForceFieldTools;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author victor
 */
class AngleBending2 extends AngleBending {

    boolean angleBending;
    int angleNumber = 0;
    int[][] angleAtomPosition = null;

    double[] v0 = null;
	double[] k2 = null;
	double[] k3 = null;
	double[] k4 = null;

    double[] v = null;
	double[] deltav = null;

    GVector moleculeCurrentCoordinates = null;
	boolean[] changeAtomCoordinates = null;
    int changedCoordinates;

    /**
     *  Set MMFF94 reference angle v0IJK and the constants k2, k3, k4 for each
     *  i-j-k angle in the molecule.
     *
     *@param  molecule       The molecule like an AtomContainer object.
     *@param  parameterSet   MMFF94 parameters set
     *@exception  Exception  Description of the Exception
     */
    @Override
    public void setMMFF94AngleBendingParameters(IAtomContainer molecule, Map parameterSet, boolean angleBendingFlag ) throws Exception {

            //logger.debug("setMMFF94AngleBendingParameters");

            IAtom[] atomConnected = null;
            angleBending=angleBendingFlag;
            for (int i = 0; i < molecule.getAtomCount(); i++) {
                    atomConnected = AtomContainerManipulator.getAtomArray(molecule.getConnectedAtomsList(molecule.getAtom(i)));
                    if (atomConnected.length > 1) {
                            for (int j = 0; j < atomConnected.length; j++) {
                                    for (int k = j+1; k < atomConnected.length; k++) {
                                            angleNumber += 1;
                                    }
                            }
                    }
            }
            //logger.debug("angleNumber = " + angleNumber);

            List angleData = null;
            MMFF94ParametersCall pc = new MMFF94ParametersCall();
            pc.initialize(parameterSet);

            v0 = new double[angleNumber];
            k2 = new double[angleNumber];
            k3 = new double[angleNumber];
            k4 = new double[angleNumber];

            angleAtomPosition = new int[angleNumber][];

            String angleType;
            IBond bondIJ = null;
            IBond bondKJ = null;
            String bondIJType;
            String bondKJType;
            int l = -1;
            for (int i = 0; i < molecule.getAtomCount(); i++) {
                    atomConnected = AtomContainerManipulator.getAtomArray(molecule.getConnectedAtomsList(molecule.getAtom(i)));
                    if (atomConnected.length > 1) {
                            for (int j = 0; j < atomConnected.length; j++) {
                                    for (int k = j+1; k < atomConnected.length; k++) {
                                            l += 1;
                                            bondIJ = molecule.getBond(atomConnected[j], molecule.getAtom(i));
                                            bondIJType = bondIJ.getProperty("MMFF94 bond type").toString();
                                            //logger.debug("bondIJType = " + bondIJType);

                                            bondKJ = molecule.getBond(atomConnected[k], molecule.getAtom(i));
                                            bondKJType = bondKJ.getProperty("MMFF94 bond type").toString();
                                            //logger.debug("bondKJType = " + bondKJType);

                                            angleType = "0";
                                            if ( bondIJType.equals("1") | (bondKJType.equals("1")) ){
                                                    angleType = "1";
                                            }
                                            if ( bondIJType.equals("1") & bondKJType.equals("1")) {
                                                    angleType = "2";
                                            }

                    angleData = pc.getAngleData(angleType, atomConnected[j].getAtomTypeName(),
                            molecule.getAtom(i).getAtomTypeName(), atomConnected[k].getAtomTypeName());

                    if( angleData!=null ){
                        v0[l] = ((Double) angleData.get(0)).doubleValue();
                        k2[l] = ((Double) angleData.get(1)).doubleValue();
                        k3[l] = ((Double) angleData.get(2)).doubleValue();
                        //k4[l] = ((Double) angleData.get(3)).doubleValue();
                    }
                    else{
                        v0[l] = 0;
                        k2[l] = 0;
                        k3[l] = 0;
                    }


                                            angleAtomPosition[l] = new int[3];
                                            angleAtomPosition[l][0] = molecule.getAtomNumber(atomConnected[j]);
                                            angleAtomPosition[l][1] = i;
                                            angleAtomPosition[l][2] = molecule.getAtomNumber(atomConnected[k]);

                                    }
                            }
                    }
            }

            v = new double[angleNumber];
            deltav = new double[angleNumber];

            this.moleculeCurrentCoordinates = new GVector(3 * molecule.getAtomCount());
            for (int i=0; i<moleculeCurrentCoordinates.getSize(); i++) {
                    this.moleculeCurrentCoordinates.setElement(i,1E10);
            }

            this.changeAtomCoordinates = new boolean[molecule.getAtomCount()];
            for (int i=0; i < molecule.getAtomCount(); i++) {
                    this.changeAtomCoordinates[i] = false;
            }

    }


    /**
     *  Calculate the actual bond angles vijk and the difference with the reference angles.
     *
     *@param  coord3d  Current molecule coordinates.
     */
    @Override
    public void setDeltav(GVector coord3d) {
            changedCoordinates = 0;
            //logger.debug("Setting Deltav");
            for (int i=0; i < changeAtomCoordinates.length; i++) {
                    this.changeAtomCoordinates[i] = false;
            }
            this.moleculeCurrentCoordinates.sub(coord3d);
            for (int i = 0; i < this.moleculeCurrentCoordinates.getSize(); i++) {
                    //logger.debug("this.moleculeCurrentCoordinates.getElement(i) = " + this.moleculeCurrentCoordinates.getElement(i));
                    if (Math.abs(this.moleculeCurrentCoordinates.getElement(i)) > 0) {
                            changeAtomCoordinates[i/3] = true;
                            changedCoordinates = changedCoordinates + 1;
                            //logger.debug("changeAtomCoordinates[" + i/3 + "] = " + changeAtomCoordinates[i/3]);
                            i = i + (2 - i % 3);
                    }
            }
            //logger.debug("currentCoordinates_deltav.length = " + currentCoordinates_deltav.length);

            for (int i = 0; i < angleNumber; i++) {
                    if ((changeAtomCoordinates[angleAtomPosition[i][0]] == true) |
                            (changeAtomCoordinates[angleAtomPosition[i][1]] == true) |
                            (changeAtomCoordinates[angleAtomPosition[i][2]] == true))		{

                            v[i] = ForceFieldTools.angleBetweenTwoBondsFrom3xNCoordinates(coord3d,angleAtomPosition[i][0],angleAtomPosition[i][1],angleAtomPosition[i][2]);
                            //logger.debug("currentCoordinates_v[" + i + "] = " + currentCoordinates_v[i]);
                            //logger.debug("v0[" + i + "] = " + v0[i]);
                            deltav[i] = v[i] - v0[i];
                            if (deltav[i] > 0 & angleBending) {
                                    deltav[i]= (-1) * deltav[i];
                            }else if (deltav[i] < 0 & !angleBending){
                                    deltav[i]= (-1) * deltav[i];
                            }

                    }

            }
            moleculeCurrentCoordinates.set(coord3d);
    }


}
