/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.tabu.impl;

import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.jama.Matrix;
import org.optiprot.maths.CalcRmsd;
import org.optiprot.metaheuristic.IMetaheuristicSolution;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.de.DifferentialEvolution;
import org.optiprot.metaheuristic.de.impl.BinCrossover;
import org.optiprot.metaheuristic.de.impl.DockScoreEval;
import org.optiprot.metaheuristic.de.impl.RandMutation;
import org.optiprot.metaheuristic.meme.IMemeSolution;
import org.optiprot.potential.docking.Docked;

/**
 * docked conformation solution
 *
 * @author victor
 */
public class DockConfSolution implements IMemeSolution {

    private int [] parameters=null;
    private Double fitness=null;
    private Docked docked=null;
    private int modeOpt=0;

    DifferentialEvolution diffEvolution=null;

    //DE conf. parameters
    private static double DE_CR=0.9;
    private static double DE_F=0.8;
    private static int DE_NP=40;

    public static final int MODE_OPT_POSE=1;//pose optimization

    //DE operators
    private static BinCrossover binCross=new BinCrossover();
    private static RandMutation randMut=new RandMutation();

    //DE parameters bounds
    private static double [] lowersLim=new double [] { -1, -1, -1, -1, -0.2, -0.2, -0.2 };
    private static double [] uppersLim=new double [] { 1, 1, 1, 1, 0.2, 0.2, 0.2};


    /**
     * create a new random solution
     *
     * @param dock
     * @param mode_opt
     * @param evolutions : initial evolutions of DE
     */
    public DockConfSolution( Docked dock, int mode_opt, int evolutions ) {
        setDocked(dock);
        modeOpt=mode_opt;

        createParameters();

        optimize(evolutions);
    }

    private DockConfSolution(DockConfSolution init) {
        setDocked(init.getDocked());
        this.modeOpt=init.getModeOpt();

        //copy parameters
        this.parameters=new int [init.size()];
        System.arraycopy(init.parameters, 0, this.parameters, 0, this.parameters.length );
    }

    /**
     * copy the parameters of other to this.
     * note, this method dont reset the internals
     *
     * @param other
     */
    public void copyParameters( DockConfSolution other ){
        fitness=null;
        System.arraycopy(other.parameters, 0, this.parameters, 0, this.parameters.length );
    }

    /**
     * reset with another dock and initializes randomly
     *
     * @param dock
     * @param evolutions : initial evolutions of DE
     */
    public void reset( Docked dock, int evolutions ) {
        setDocked(dock);

        createParameters();

        if( diffEvolution!=null )
            diffEvolution.reset();

        optimize(evolutions);
    }

    /**
     * create a new solution in the neighborhood of this
     *
     * @param init
     * @return
     */
    public DockConfSolution newNeighbor(){

        DockConfSolution newsol=new DockConfSolution(this);

        mutate(newsol);

        return newsol;
    }

    private static void mutate( DockConfSolution sol ){

        if( sol.getModeOpt()==MODE_OPT_POSE ){

            Integer id_act=null;

            // fix active atom or active point randomly
            int fixed=0;

            if( Math.random()<0.5)
                fixed=1;

            if( fixed==0 )
                id_act=sol.getDocked().getCompActPoint( sol.getParameter(0) );
            else
                id_act=sol.getDocked().getCompActAtom( sol.getParameter(1) );


            sol.setParameter( (fixed+1)%2, id_act);
        }
    }


    /**
     * mutate the current solution,
     * 
     * @return : mutated solution (this)
     */
    public IMemeSolution mutate() {
        fitness=null;
        mutate(this);
        return this;
    }


    /**
     * do more iterations of Diff. Evolution
     *
     * @param evolutions
     */
    public void optimize(int evolutions){

        if( evolutions<1 )
            return;

        if( diffEvolution==null )
            diffEvolution=createDE();

        fitness=null;
        getDocked().assignAtomToActivePoint( getParameter(0), getParameter(1) );
        diffEvolution.setIterLimit( diffEvolution.getIterLimit()+evolutions );
        diffEvolution.doSearch();
    }

    public Double getFitness() {

        if( fitness==null ){
            getDocked().assignAtomToActivePoint( getParameter(0), getParameter(1) );
            fitness=diffEvolution.getBestSolution().getFitness();
        }

        return fitness;
    }

    /**
     * force calc of fitness
     *
     * @return
     */
    public Double recalcFitness() {

        getDocked().assignAtomToActivePoint( getParameter(0), getParameter(1) );
        fitness = ((DESolution)diffEvolution.getBestSolution()).recalcFitness();

        return fitness;
    }

    public boolean isBetter(IMetaheuristicSolution other) {
        
        if( this.getFitness()<other.getFitness() )
            return true;

        return false;
    }

    public int getParameter( int idx ){
        return parameters[idx];
    }

    protected void setParameter( int idx, int val ){

        if( parameters[idx]==val )
            return;
        
        fitness=null;
        parameters[idx]=val;

        if( diffEvolution!=null )
            diffEvolution.invalidate();
    }

    public int size(){
        return parameters.length;
    }

    @Override
    public boolean equals(Object obj) {

        DockConfSolution other=(DockConfSolution) obj;

        for(int i=0;i<parameters.length;i++){
            if( other.parameters[i]!=this.parameters[i] )
                return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 71 * hash + (this.parameters != null ? this.parameters.hashCode() : 0);
        return hash;
    }

    /**
     * @return the docked
     */
    public Docked getDocked() {
        return docked;
    }

    /**
     * @param docked the docked to set
     */
    protected void setDocked(Docked docked) {
        this.docked = docked;
    }

    /**
     * 
     * @return
     */
    private DifferentialEvolution createDE(){

        DifferentialEvolution instance = new DifferentialEvolution();

        instance.setCR( DE_CR );
        instance.setF( DE_F );
        instance.setNP( DE_NP );
        instance.setCrossover( binCross );
        instance.setMutation( randMut );
        instance.setBounds(lowersLim, uppersLim);

        instance.createRandomPopulation( DockScoreEval.NUM_PARAMETERS, getDocked().getEvaluator() );

        return instance;

    }

    /**
     * @return the modeOpt
     */
    public int getModeOpt() {
        return modeOpt;
    }

    /**
     * create the vector of parameters, the length differs with mode of
     * optimization, after initializes it qith random values
     */
    private void createParameters() {

        if( this.modeOpt==MODE_OPT_POSE ){

            if( parameters==null)
                parameters=new int [2];

            Integer id_actp=null;

            do{
                parameters[0] =
                    (int)Math.floor( Math.random()*getDocked().getNumActAtoms() );

                id_actp=getDocked().getCompActPoint( parameters[0] );
            }
            while( id_actp==null );

            parameters[1] = id_actp;
        }
    }

    /**
     * gets the active point fixed by the solution
     *
     * @return
     */
    public int getActivePoint(){
        return getParameter(1);
    }

    /**
     * gets the ligand active atom fixed by the solution
     * 
     * @return
     */
    public int getLigActiveAtom(){
        return getParameter(0);
    }

    /**
     * gets the matrix calculated by the DE
     *
     * @return
     */
    public Matrix getMTrans() {

        DESolution desol=(DESolution) diffEvolution.getBestSolution();
        DockScoreEval evaluator=(DockScoreEval) desol.getEvaluator();

        return evaluator.getMTrans(desol);
    }

    /**
     * checks if the ligand encoded by other is equals to the ligand
     * encoded by this using a rmsd threshold
     * 
     * @param other
     * @param rmsdLim
     * @return
     */
    public boolean equalConformation( DockConfSolution other, double rmsdLim ){

        getDocked().assignAtomToActivePoint( this.getLigActiveAtom(), this.getActivePoint());
        getDocked().setTransform( this.getMTrans() );
        List<Atom> l1=getDocked().getLigandAtomsCopy1();

        getDocked().assignAtomToActivePoint( other.getLigActiveAtom(), other.getActivePoint());
        getDocked().setTransform( other.getMTrans() );
        List<Atom> l2=getDocked().getLigandAtomsCopy2();

        double rmsd=CalcRmsd.bruteRmsd( l1, l2);

        if( rmsd < rmsdLim )
            return true;

        return false;
    }

}
