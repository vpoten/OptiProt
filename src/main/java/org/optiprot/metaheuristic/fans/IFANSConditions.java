/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.fans;

/**
 *
 * @author victor
 */
public interface IFANSConditions {

    public boolean endCondition( FANS fans );

    public void restart( FANS fans );

    public boolean searchBlocked( FANS fans );

    public void update( FANS fans );

    public int getIterations();

    public void reset( FANS fans );
}
