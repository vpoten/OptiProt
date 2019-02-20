/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de;

import java.util.ArrayList;
import org.optiprot.metaheuristic.fans.FANSDouble;

/**
 *
 * @author victor
 */
public class DEFansWork implements Runnable {

    private int period=0;
    private int index=-1;
    
    private IDESolutionEval evaluator=null;
    private FANSDouble fans=null;
    private int iterations=50;

    private ArrayList<DESolution> listSol=null;


    public DEFansWork( int index, int period, FANSDouble fans, IDESolutionEval eval, int iter) {
        this.period=period;
        this.index=index;
        this.fans=fans;
        this.evaluator=eval;
        this.iterations=iter;
    }


    
    public void run() {

        for( int i=this.index; i<listSol.size(); i+=this.period ){
            DESolution currSol=listSol.get(i);
            fans.reset();
            fans.setSolution( currSol );
            ((DESolution)fans.getSolution()).setEvaluator(evaluator);
            ((DESolution)fans.getBestSolution()).setEvaluator(evaluator);
            fans.setIterLimit( fans.getIterLimit()+iterations );
            fans.doSearch();
            currSol.copy( (DESolution) fans.getBestSolution());
        }
    }

    /**
     * @param listSol the listSol to set
     */
    public void setListSol(ArrayList<DESolution> listSol) {
        this.listSol = listSol;
    }

    

}
