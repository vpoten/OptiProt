/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

/**
 * regular fuzzy set defined in normalized domain [0,1]
 *
 * @author victor
 */
public class RegFuzzySet {

    static private double [] auxVector=new double [3];

    /**
     * Calculates the value of menbership for all outputs.length different
     * triangular functions over [0,1]
     *
     * @param in : input value
     * @param outputs : I/O
     */
    static public void menbership(double in, double [] outputs){

        double range=1.0/(outputs.length-1.0);

        auxVector[0]=0;
        auxVector[1]=0;
        auxVector[2]=range;
        outputs[0]=membership( in, auxVector );

        for(int i=1; i<(outputs.length-1); i++){
            auxVector[0]=range*(i-1);
            auxVector[1]=range*i;
            auxVector[2]=range*(i+1);
            outputs[i]=membership( in, auxVector );
        }

        auxVector[0]=1-range;
        auxVector[1]=1;
        auxVector[2]=1;
        outputs[outputs.length-1]=membership( in, auxVector );
    }


    /**
     *
     * @param in
     * @param parameters min, mid, and max of a triangular function
     * @return
     */
    static private double membership(double in, double [] parameters) {
		// Outside range? => membership is 0
		if( (in < parameters[0]) || (in > parameters[2]) ) return 0;

		// Middle point of the triangle? => membership is 1.0
		if( in == parameters[1] ) return 1.0;

		// Value between 'min' and 'mid'
		if( in < parameters[1] ) return ((in - parameters[0]) / (parameters[1] - parameters[0]));

		// Value between 'mid' and 'max'
		return 1 - ((in - parameters[1]) / (parameters[2] - parameters[1]));
	}

}
