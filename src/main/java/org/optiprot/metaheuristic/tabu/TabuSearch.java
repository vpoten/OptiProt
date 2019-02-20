/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.tabu;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Date;
import java.util.Deque;
import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.metaheuristic.MetaheuristicBase;

/**
 *
 * @author victor
 */
public class TabuSearch extends MetaheuristicBase {


    private Collection<IMetaheuristicSolution> tabuList =
            new ArrayDeque<IMetaheuristicSolution>();

    private int tenure=0;

    private ITabuNeighManager neighManager=null;



    public void doSearch() {

        if( getNeighManager()==null || getTenure()==0 ){
            throw new RuntimeException("Tabu bad initialization.");
        }

        setStartTime( new Date() );

        while( !endCondition() ){

            IMetaheuristicSolution newsol =
                    getNeighManager().generate(getSolution(), getTabuList());

            setSolution(newsol);
            updateTabuList(newsol);
        }

        setEndTime( new Date() );
    }

    public double getLambda() {
        return 1.0;//nothing to do
    }

    public void setLambda(double val) {
        return;//nothing to do
    }


    /**
     * @return the tabuList
     */
    protected Collection<IMetaheuristicSolution> getTabuList() {
        return tabuList;
    }

    /**
     *
     * @param newsol
     */
    private void updateTabuList(IMetaheuristicSolution newsol) {
        
        if( getTabuList().size()==getTenure() ){
            ((Deque)getTabuList()).poll();
        }

        getTabuList().add(newsol);
    }

    /**
     * @return the tabu tenure
     */
    public int getTenure() {
        return tenure;
    }

    /**
     * @param tenure the tabu tenure to set
     */
    public void setTenure(int tenure) {
        this.tenure = tenure;
    }

    /**
     * @return the neighManager
     */
    public ITabuNeighManager getNeighManager() {
        return neighManager;
    }

    /**
     * @param neighManager the neighManager to set
     */
    public void setNeighManager(ITabuNeighManager neighManager) {
        this.neighManager = neighManager;
    }

    @Override
    public void reset(){

        super.reset();

        if( !getTabuList().isEmpty() )
            getTabuList().clear();
    }

}
