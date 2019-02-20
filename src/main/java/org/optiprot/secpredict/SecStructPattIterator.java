/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

import java.util.Iterator;
import java.util.List;
import org.optiprot.neural.NNPattern;

/**
 * Iterator that generates the patterns on the fly using the list of chain residues
 *
 * @author victor
 */
public class SecStructPattIterator implements Iterator<NNPattern> {

    private List<SecStructSubPattern> listSubPat=null;

    private SecStructSubPattern [] window=null;
    private int resPredPos;

    private NNPattern currPatt=null;
    private ISecStructPatternGen patGen=null;

    private int counter=0;


    /**
     *
     * @param list : list of chain subpatterns
     * @param windowSize : size of the window
     * @param resPredPos : position of the predicted residue inside the window
     * @param patgen : pattern generator
     * @param reusePat : if true reuse the generated pattern
     */
    public SecStructPattIterator( List<SecStructSubPattern> list, int windowSize,
            int resPredPos, ISecStructPatternGen patgen, boolean reusePat ) {
        this.listSubPat=list;

        this.resPredPos=resPredPos;
        this.patGen=patgen;

        this.window=new SecStructSubPattern [windowSize];

        if( reusePat )
            this.currPatt=this.patGen.createPattern(windowSize);

        for(int i=0;i<this.window.length;i++)
            this.window[i]=null;

        //initialize window
        int j=0;
        for( int i=this.resPredPos+1; i<this.window.length; i++){
            if( j<listSubPat.size() )
                this.window[i]=listSubPat.get(j++);
        }
    }


    public boolean hasNext() {
        return ( counter<listSubPat.size() );
    }


    public NNPattern next() {
        SecStructSubPattern subPat=null;
        
        int pos=counter+this.window.length-this.resPredPos-1;
        
        if( pos<this.listSubPat.size() )
            subPat=this.listSubPat.get(pos);

        for(int i=1; i<this.window.length; i++)
            this.window[i-1]=this.window[i];

        this.window[this.window.length-1]=subPat;

        NNPattern pat=this.currPatt;

        if( currPatt!=null ){//if reuse pattern
            this.patGen.genPattern(window, resPredPos, currPatt);
        }
        else{
            pat=this.patGen.genPattern(window, resPredPos);
        }

        counter++;
        return pat;

    }

    
    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
