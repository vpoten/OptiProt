/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans;

import java.util.Date;
import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.metaheuristic.MetaheuristicBase;

/**
 * FANS metaheuristic framework
 *
 * @author victor
 */
public class FANS extends MetaheuristicBase {

    
    private IFANSFuzzyEvaluator evaluator=null;
    private IFANSConditions conditions=null;
    private IFANSNeighManager neighManager=null;
    private IFANSOperatorManager operManager=null;
    private IFANSOperator operator=null;
    

    public FANS() {
    }

    
    public void doSearch(){

        if( getSolution()==null || getConditions()==null || getOperator()==null
                || getNeighManager()==null || getOperManager()==null || getEvaluator()==null ){
            throw new RuntimeException("FANS bad initialization.");
        }

        IMetaheuristicSolution newSolution=null;

        setStartTime( new Date() );

        while( !this.getConditions().endCondition(this) ){

            //execute neighbor manager
            newSolution = this.getNeighManager().generate( getOperator(),
                    getEvaluator(), getSolution() );

            // if the new solution is good accept it
            if( newSolution!=null ){

                this.setSolution(newSolution);
                this.getEvaluator().adjust(newSolution);
            }
            else{
                //if the new solution is not good => modify the operator
                this.getOperManager().modify( this.getOperator() );
            }

            // if we are working in cooperative mode
            if( this.getCoordinator()!=null ){
                this.getCoordinator().sendReport( getId(), getSolution(), new Date() );
            }

            this.getConditions().update( this );

            if( this.getConditions().searchBlocked(this) ){
                //if the search is blocked => do a restart/escape action
                this.getConditions().restart( this );
            }

        }//

        setEndTime( new Date() );
    }

   

    /**
     * @return the evaluator
     */
    public IFANSFuzzyEvaluator getEvaluator() {
        return evaluator;
    }

    /**
     * @param evaluator the evaluator to set
     */
    public void setEvaluator(IFANSFuzzyEvaluator evaluator) {
        this.evaluator = evaluator;
    }

    /**
     * @return the condition
     */
    public IFANSConditions getConditions() {
        return conditions;
    }

    /**
     * @param condition the condition to set
     */
    public void setConditions(IFANSConditions condition) {
        this.conditions = condition;
    }

    /**
     * @return the neighManager
     */
    public IFANSNeighManager getNeighManager() {
        return neighManager;
    }

    /**
     * @param neighManager the neighManager to set
     */
    public void setNeighManager(IFANSNeighManager neighManager) {
        this.neighManager = neighManager;
    }

    /**
     * @return the operManager
     */
    public IFANSOperatorManager getOperManager() {
        return operManager;
    }

    /**
     * @param operManager the operManager to set
     */
    public void setOperManager(IFANSOperatorManager operManager) {
        this.operManager = operManager;
    }

    /**
     * @return the operator
     */
    public IFANSOperator getOperator() {
        return operator;
    }

    /**
     * @param operator the operator to set
     */
    public void setOperator(IFANSOperator operator) {
        this.operator = operator;
    }


    public double getLambda() {
        return this.getEvaluator().getLambda();
    }

    public void setLambda(double val) {
        this.getEvaluator().setLambda(val);
    }

    /**
     * reset the solutions
     */
    @Override
    public void reset() {
        this.getConditions().reset(this);
        this.setSolution(null);
        this.setBestSolution(null);
    }
   
}
