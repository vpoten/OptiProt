/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.docking;

import org.optiprot.potential.docking.Docked;

/**
 *
 * @author victor
 */
public interface IDockingScore {

    public double calcScore( Docked dock );

    public int getNevaluations();
}
