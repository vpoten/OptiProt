/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.util.HashMap;
import java.util.Iterator;
import org.optiprot.maths.RegFuzzySet;
import org.optiprot.neural.NNPattern;

/**
 *
 * @author victor
 */
public class SecStructPattGen1 implements ISecStructPatternGen {

    //COD1: polar, charge, vdw vol (double)
    public static final int MODE_COD1 = 1;

    //COD2: vdw vol, hydrophobic (double)
    public static final int MODE_COD2 = 2;

    //COD3: vdw vol, hydrophobic (binary-fuzzy) 5 sets
    public static final int MODE_COD3 = 3;

    //COD3: vdw vol, hydrophobic (binary-fuzzy) 3 sets
    public static final int MODE_COD4 = 4;

    public static final double ZERO_VALUE=1e-4;

    private int mode = MODE_COD1;
    private int numInputs=4;//num of inputs per residue
    private int numOutputs=3;

    //table for amino acids features
    private static HashMap<String, double []> aminoAFeat =
            new HashMap<String, double[]>();

    //list of amino acid codes
    private static String [] AminoAcids = {
                            "ala", "arg", "asn", "asp", "cys",
                            "gln", "glu", "gly", "his", "ile",
                            "leu", "lys", "met", "phe", "pro",
                            "ser", "thr", "trp", "tyr", "val",
                            "sec" };

    static {
        //create table of features

        // features: double [] (3)
        // 1 - polar : 0 nonpolar, 1 polar
        // 2 - charge : 0 negative, 0.5 neutral, 1 positive
        // 3 - VDW volume : O little, 0.5 average, 1 large
        // 4 - hydrophobic value

        double max_vol=163.0;//max VdW volume
        double min_vol=48.0;//min VdW volume

        double max_hyd=0.73;//max hydrophobic
        double min_hyd=-1.8;//min hydrophobic
        

        aminoAFeat.put( AminoAcids[0],
                new double[] { 0, 0.5, 67.0, 0.25 } );//ala
        aminoAFeat.put( AminoAcids[1], 
                new double[] { 1, 1, 148.0, -1.8 } );//arg
        aminoAFeat.put( AminoAcids[2], 
                new double[] { 1, 0.5, 96.0, -0.64 } );//asn
        aminoAFeat.put( AminoAcids[3], 
                new double[] { 1, 0, 91.0, -0.72 } );//asp
        aminoAFeat.put( AminoAcids[4], 
                new double[] { 0, 0.5, 86.0, 0.04 } );//cys
        aminoAFeat.put( AminoAcids[5], 
                new double[] { 1, 0.5, 114.0, -0.69 } );//gln
        aminoAFeat.put( AminoAcids[6], 
                new double[] { 1, 0, 109.0, -0.62 } );//glu
        aminoAFeat.put( AminoAcids[7], 
                new double[] { 0, 0.5, 48.0, 0.16 } );//gly
        aminoAFeat.put( AminoAcids[8], 
                new double[] { 1, 0.9, 118.0, -0.4 } );//his
        aminoAFeat.put( AminoAcids[9], 
                new double[] { 0, 0.5, 124.0, 0.73 } );//ile
        aminoAFeat.put( AminoAcids[10], 
                new double[] { 0, 0.5, 124.0, 0.53 } );//leu
        aminoAFeat.put( AminoAcids[11], 
                new double[] { 1, 1, 135.0, -1.1 } );//lys
        aminoAFeat.put( AminoAcids[12], 
                new double[] { 0, 0.5, 124.0, 0.26 } );//met
        aminoAFeat.put( AminoAcids[13], 
                new double[] { 0, 0.5, 135.0, 0.25 } );//phe
        aminoAFeat.put( AminoAcids[14], 
                new double[] { 0, 0.5, 90.0, -0.07 } );//pro
        aminoAFeat.put( AminoAcids[15], 
                new double[] { 1, 0.5, 73.0, -0.26 } );//ser
        aminoAFeat.put( AminoAcids[16], 
                new double[] { 1, 0.5, 93.0, -0.18 } );//thr
        aminoAFeat.put( AminoAcids[17], 
                new double[] { 0, 0.5, 163.0, 0.37 } );//trp
        aminoAFeat.put( AminoAcids[18], 
                new double[] { 1, 0.5, 141.0, -0.02 } );//tyr
        aminoAFeat.put( AminoAcids[19], 
                new double[] { 0, 0.5, 105.0, 0.54 } );//val
        aminoAFeat.put( AminoAcids[20], 
                new double[] { 0, 0.5, 88.0, 0.04} );//sec

        //normalize volume and hydrophobity
        Iterator<double[]> it=aminoAFeat.values().iterator();

        while( it.hasNext() ){
            double [] array=it.next();

            array[2] = (array[2]-min_vol)/(max_vol-min_vol);
            array[3] = (array[3]-min_hyd)/(max_hyd-min_hyd);
        }
    }

    private double [] auxInputs=new double [0];
    private double [] auxOutputs=null;

    public SecStructPattGen1() {
        
    }

    private SecStructPattGen1(SecStructPattGen1 other) {
        
        this.setMode(other.getMode());

        if( other.auxInputs!=null )
            this.auxInputs=new double [other.auxInputs.length];

        this.auxOutputs=new double [other.numOutputs];
    }

    
    public NNPattern createPattern(int windowSize) {

        if( auxInputs.length != (windowSize*numInputs) ){
            auxInputs=new double [windowSize*numInputs];
            auxOutputs=new double [numOutputs];
        }

        return new NNPattern( new double [windowSize*numInputs], new double [numOutputs]);
    }


    public void genPattern(SecStructSubPattern[] window, int resPredPos, NNPattern pat) {

        double [] outputs=null;

        if( getMode()==MODE_COD3 || getMode()==MODE_COD4 )
            //auxiliary vector for fuzzy transformation
            outputs=new double[(numInputs-1)/2];


        for(int i=0; i<window.length; i++){

            double [] feats = null;

            if( window[i]!=null )
                feats=aminoAFeat.get( SecStructPattManager.getAAName(window[i].aminoacid) );

            if( feats==null ){
                for( int j=0; j<numInputs; j++)
                    auxInputs[i*numInputs+j]=0;
            }
            else{

                auxInputs[i*numInputs]=1;

                if( getMode()==MODE_COD1 ){
                    auxInputs[i*numInputs+1]=feats[0];
                    auxInputs[i*numInputs+2]=feats[1];
                    auxInputs[i*numInputs+3]=feats[2];
                }
                else if( getMode()==MODE_COD2 ){
                    auxInputs[i*numInputs+1]=feats[2];
                    auxInputs[i*numInputs+2]=feats[3];
                }
                else if( getMode()==MODE_COD3 || getMode()==MODE_COD4 ){

                    RegFuzzySet.menbership(feats[2], outputs);

                    for( int j=0; j<outputs.length; j++)
                        auxInputs[i*numInputs+1+j]=outputs[j];

                    RegFuzzySet.menbership(feats[3], outputs);

                    for( int j=0; j<outputs.length; j++)
                        auxInputs[i*numInputs+1+outputs.length+j]=outputs[j];
                }
            }

        }

        setPattValues( window[resPredPos], pat );
        
    }


    private void setPattValues( SecStructSubPattern predict, NNPattern pat){
        //outputs: 0 helix, 1 strand, 2 coil

        for(int i=0; i<auxInputs.length; i++){
            //prevent for values==0 (problem with PUNN)
            if( auxInputs[i]<ZERO_VALUE )
                auxInputs[i]=ZERO_VALUE;
        }

        auxOutputs[0]=0.1;
        auxOutputs[1]=0.1;
        auxOutputs[2]=0.1;

        auxOutputs[ predict.type ]=0.9;

        pat.setInputs(auxInputs);
        pat.setOutputs(auxOutputs);
    }


    public NNPattern genPattern(SecStructSubPattern[] window, int resPredPos) {
        NNPattern pat=createPattern(window.length);
        genPattern( window, resPredPos, pat);
        return pat;
    }

    
    /**
     * @return the mode
     */
    public int getMode() {
        return mode;
    }

    /**
     * @param mode the mode to set
     */
    public void setMode(int mode) {
        this.mode = mode;

        if( getMode()==MODE_COD1 ){
            numInputs=4;
            numOutputs=3;
        }
        else if( getMode()==MODE_COD2 ){
            numInputs=3;
            numOutputs=3;
        }
        else if( getMode()==MODE_COD3 ){
            numInputs=11;
            numOutputs=3;
        }
        else if( getMode()==MODE_COD4 ){
            numInputs=7;
            numOutputs=3;
        }
    }

    @Override
    public ISecStructPatternGen clone() {
        return new SecStructPattGen1(this);
    }



}
