/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans.impl;

import org.optiprot.aacids.AABasicFactory;
import org.optiprot.metaheuristic.fans.FANS;
import org.optiprot.metaheuristic.fans.IFANSConditions;
import org.optiprot.metaheuristic.fans.IFANSOperator;
import org.optiprot.metaheuristic.fans.IFANSOperatorManager;
import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.maths.CalcTransform;

/**
 * implements OperatorManager, FANSOperator and FANSConditions for proteins
 *
 * @author victor
 */
public class ProtFANSOperator implements IFANSOperatorManager, IFANSOperator,
    IFANSConditions {

    private int kmodify=0;///number of positions to modify
    private int withoutImprove=0;//iterations without improve

    private FANS fans=null;

    private int withoutImprLimit=0;//limit of iterations without improve

    
    static private double K_START_RATE = 0.25;


    protected enum OperationMode {SEGMENT,FLIP};

    /////////////////////////////

    /**
     *
     * @param length : length of the chain
     * @param p_iter : limit iteration (0 if not set)
     * @param p_time : limit seconds (0 if not set)
     */
    public ProtFANSOperator(FANS fans, int length, int p_iter, int p_time) {
        kmodify = (int) (length * K_START_RATE);
        fans.setIterLimit( p_iter );
        fans.setTimeLimit( p_time );
        this.fans=fans;
        withoutImprLimit=(int) (length * K_START_RATE);
    }

    
    /**
     * modify the operator (IFANSOperatorManager)
     *
     * @param operator to modify (self)
     */
    public void modify(IFANSOperator operator) {

        if( !(operator instanceof ProtFANSOperator) ){
            throw new ClassCastException("bad IFANSOperator.");
        }

        if(kmodify>0)
            kmodify--;
    }

    /**
     * modify the solution (IFANSOperator)
     *
     * @param : solution to modify
     * @return : a new solution
     */
    public IMetaheuristicSolution modify(IMetaheuristicSolution solution) {

        if( !(solution instanceof ProtFANSSolution) ){
            throw new ClassCastException("bad ProtFANSSolution.");
        }

        //choose operation mode
        OperationMode mode=OperationMode.SEGMENT;

        if( Math.random() > 0.5 )
            mode=OperationMode.FLIP;

        //create a new solution
        ProtFANSSolution currsol = (ProtFANSSolution) solution;

        ProtFANSSolution newsol = new ProtFANSSolution( 
                AABasicFactory.clone( currsol.getChain() ), currsol.getParameters() );

        int length = currsol.getChain().length;
        int aaindex = (int) Math.floor( Math.random()*length );
       
        double percent = 0;

        // the angle step and rotamer change rate is proportional
        // to modification granularity
        double angleStep = newsol.getParameters().getMutAngleStepMax() *
                (kmodify/(length*K_START_RATE));

        double rotamerProb = (kmodify/(length*K_START_RATE)) /
                (double)newsol.getParameters().getMutRateRotamer();

        //angle to modify
        boolean psi=false;
        if( Math.random()<0.5 )
            psi=true;
        
        //modify the new solution
        for(int i=0;i<kmodify;i++){

            if( mode==OperationMode.SEGMENT ){
                aaindex = (aaindex+i)%length;
            }
            else if( mode==OperationMode.FLIP ){
                aaindex = (int) Math.floor( Math.random()*length );
            }

            double sign=-1;
            if( Math.random()<0.5 )
                sign=1.0;

            if( psi)
                CalcTransform.changePsi( newsol.getChain(), aaindex, angleStep*sign);
            else
                CalcTransform.changePhi( newsol.getChain(), aaindex, angleStep*sign);
            

            if( Math.random()<rotamerProb ){
                percent = -1 + Math.random()*2;
                newsol.getChain()[aaindex].applyMutation(aaindex, percent);
            }
                
        }

        currsol=null;
        return newsol;
    }

    /**
     * checks end condition  (IFANSConditions)
     *
     * @return
     */
    public boolean endCondition( FANS fans ) {

        return fans.endCondition();
    }

    /**
     * restart the search (IFANSConditions)
     *
     * @param solution
     */
    public void restart( FANS fans ) {
        
        if( !(fans.getSolution() instanceof ProtFANSSolution) ){
            throw new ClassCastException("bad ProtFANSSolution.");
        }
        
        ProtFANSSolution protSol=(ProtFANSSolution) fans.getSolution();

        //create a random solution
        kmodify = (int) (protSol.getChain().length * 0.5);

        protSol=(ProtFANSSolution) this.modify( protSol );
        fans.setSolution(protSol);

        //restart the parameter k and others
        reset( fans );
    }

    public void reset( FANS fans ) {
        ProtFANSSolution protSol=(ProtFANSSolution) fans.getSolution();
        kmodify = (int) (protSol.getChain().length * K_START_RATE);
        withoutImprove = 0;
    }

    /**
     * checks if the search is blocked (IFANSConditions)
     * @return
     */
    public boolean searchBlocked( FANS fans ) {

        if( kmodify==0 )
            return true;

        if( withoutImprove>=withoutImprLimit )
            return true;

        return false;
    }

    /**
     * update the conditions internals (IFANSConditions)
     * @param fans
     */
    public void update(FANS fans) {
        
        if( !fans.getSolution().isBetter( fans.getBestSolution()) ){
            withoutImprove++;
        }
        else{
            withoutImprove=0;
        }
    }

    public int getIterations() {
        return fans.getIterations();
    }
    
}
