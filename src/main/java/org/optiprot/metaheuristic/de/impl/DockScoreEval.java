/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de.impl;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.jama.Matrix;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.Quaternion;
import org.optiprot.metaheuristic.de.DESolution;
import org.optiprot.metaheuristic.de.IDESolutionEval;
import org.optiprot.potential.docking.IDockingScore;
import org.optiprot.potential.docking.Docked;

/**
 * this class evaluates a differential evolution solution applied to docking
 *
 * q = quaternion that encodes a rotation
 * solution = { q0, q1, q2, q3, tx, ty, tz }
 * 
 * @author victor
 */
public class DockScoreEval implements IDESolutionEval {

    private IDockingScore scoreFunc=null;
    private Docked docked=null;

    private Atom axis=new AtomImpl();
    private Atom shift=new AtomImpl();
    private double angle=0;

    public static final int NUM_PARAMETERS = 7;

    public DockScoreEval() {
    }

    public DockScoreEval( IDockingScore scoreFunc ) {
        this.setScoreFunc(scoreFunc);
    }
    

    /**
     * extract the parameters encoded by the DE solution
     * 
     * @param sol
     */
    private void decodeSolution(DESolution sol){

        // solution = { q0, q1, q2, q3, tx, ty, tz }

        Quaternion q=new Quaternion(sol.getParameter(0), sol.getParameter(1),
                sol.getParameter(2), sol.getParameter(3));

        angle=q.getAngle();

        q.getAxis( axis.getCoords() );

        //normalize axis
        CalcGeom.product2(axis, 1.0/Calc.amount(axis));

        shift.setX( sol.getParameter(4) );
        shift.setY( sol.getParameter(5) );
        shift.setZ( sol.getParameter(6) );

    }
    
    public Double getFitness(DESolution sol) {

        decodeSolution(sol);

        getDocked().setTransform(axis, angle, shift);
        return getScoreFunc().calcScore( getDocked() );
    }

    /**
     * returns true if val1 is better than val2
     *
     * @param val1
     * @param val2
     * @return
     */
    public boolean isBetter(double val1, double val2) {

        if(val1<val2)
            return true;

        return false;
    }

    /**
     * @return the scoreFunc
     */
    public IDockingScore getScoreFunc() {
        return scoreFunc;
    }

    /**
     * @param scoreFunc the scoreFunc to set
     */
    public void setScoreFunc(IDockingScore scoreFunc) {
        this.scoreFunc = scoreFunc;
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
    public void setDocked(Docked docked) {
        this.docked = docked;
    }

    /**
     * get the Matrix encoded by the DE solution
     *
     * @param sol
     * @return a new Matrix
     */
    public Matrix getMTrans(DESolution sol) {
        
        decodeSolution(sol);

        Matrix mat=Matrix.identity(4, 4);
        Docked.setTransform( mat, axis, angle, shift);

        return mat;
    }

    public int getNevaluations() {
        return this.getScoreFunc().getNevaluations();
    }
}
