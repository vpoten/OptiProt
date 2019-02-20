/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot;

import java.util.Date;
import java.util.List;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.optiprot.aacids.AABasicFactory;
import org.optiprot.metaheuristic.MetaCooperative;
import org.optiprot.metaheuristic.fans.FANS;
import org.optiprot.metaheuristic.fans.impl.ProtFANSNeighbor;
import org.optiprot.metaheuristic.fans.impl.ProtFANSOperator;
import org.optiprot.metaheuristic.fans.impl.ProtFANSSolution;

/**
 *
 * @author victor
 */
public class OptiProtFANSRun {

    /**
     * single fans run
     * 
     * @param sequence
     * @param initChain : Chain to optimize (random if null)
     * @param p_parameters
     * @param p_secondsLimit
     * @param p_iterLimit
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
     static public OptiProtResults run( final String sequence, Chain initChain,
             final OptiProtParameters p_parameters, final int p_secondsLimit,
             final int p_iterLimit )
            throws StructureException {

         OptiProtResults res=new OptiProtResults();
         ProtFANSSolution initSol=null;

         if( initChain==null ){
            //create a random init solution
            initSol = new ProtFANSSolution(
                AABasicFactory.create( sequence, p_parameters),
                p_parameters );
         }
         else{
             initSol = new ProtFANSSolution(
                AABasicFactory.create(initChain, p_parameters),
                p_parameters );
         }

         //set the init solution in results
        res.setInitSolution( AABasicFactory.toChain( initSol.getChain(), false) );
        res.setInitFitness( initSol.getFitness() );

        FANS instance = new FANS();

        instance.setSolution(initSol);

        int maxTrials=5;
        
        ProtFANSNeighbor neigFuzzMang =
                new ProtFANSNeighbor( 0.2, maxTrials, initSol.getFitness() );
        neigFuzzMang.setLambda( 0.9 );

        instance.setEvaluator(neigFuzzMang);
        instance.setNeighManager(neigFuzzMang);

        int length=initSol.getChain().length;

        ProtFANSOperator operMang =
                new ProtFANSOperator(instance, length , p_iterLimit, p_secondsLimit);

        instance.setConditions(operMang);
        instance.setOperator(operMang);
        instance.setOperManager(operMang);

        instance.doSearch();

        ProtFANSSolution bestSol=(ProtFANSSolution) instance.getBestSolution();
        bestSol.center();

        //set the results
        res.setBestSolution( AABasicFactory.toChain( bestSol.getChain(), false) );
        res.setBestFitness( bestSol.getFitness() );

        res.setStartTime( instance.getStartTime() );
        res.setEndTime( instance.getEndTime() );

        res.setIterations( operMang.getIterations() );

        return res;

     }

     /**
      * fans cooperative run
      *
      * @param sequence
      * @param initChain : Chain to optimize (random if null)
      * @param listPar : parameters for each solver
      * @param nthreads : number of solvers
      * @param p_secondsLimit
      * @param p_iterLimit
      * @return
      * @throws org.biojava.bio.structure.StructureException
      * @throws java.lang.Exception
      */
     static public OptiProtResults run( final String sequence, Chain initChain,
             List<OptiProtParameters> listPar, final int nthreads, final int p_secondsLimit,
             final int p_iterLimit )
            throws StructureException, Exception {
         
         MetaCooperative coordinator =  new MetaCooperative();

         OptiProtResults res=new OptiProtResults();
         res.setStartTime( new Date() );


         for( int i=0; i<nthreads; i++ ){

             ProtFANSSolution initSol=null;

             OptiProtParameters par = listPar.get(i);

             if( initChain==null ){
                //create a random init solution
                initSol = new ProtFANSSolution(
                    AABasicFactory.create( sequence, par),
                    par );
             }
             else{
                 initSol = new ProtFANSSolution(
                    AABasicFactory.create(initChain, par),
                    par );
             }

             //set the init solution in results
             if( res.getInitSolution()==null ){
                res.setInitSolution( AABasicFactory.toChain( initSol.getChain(), false) );
                res.setInitFitness( initSol.getFitness() );
             }

             //create and configure each solver
            FANS instance = new FANS();

            instance.setSolution(initSol);

            int maxTrials=5;

            ProtFANSNeighbor neigFuzzMang =
                    new ProtFANSNeighbor( 0.2, maxTrials, initSol.getFitness() );
            
            //set differents lambda values
            neigFuzzMang.setLambda( 0.95-0.05*i );

            instance.setEvaluator(neigFuzzMang);
            instance.setNeighManager(neigFuzzMang);

            int length=initSol.getChain().length;

            ProtFANSOperator operMang =
                    new ProtFANSOperator(instance, length , p_iterLimit, p_secondsLimit);

            instance.setConditions(operMang);
            instance.setOperator(operMang);
            instance.setOperManager(operMang);

            coordinator.addSolver(instance);

         }

         //run the cooperative threads
         MetaCooperative.runSolvers(coordinator, 10000);

         ProtFANSSolution bestSol=(ProtFANSSolution) coordinator.getBestSolution();
         bestSol.center();

         //set the results
         res.setBestSolution( AABasicFactory.toChain( bestSol.getChain(), false) );
         res.setBestFitness( bestSol.getFitness() );
         res.setEndTime( new Date() );

         return res;
     }
}
