/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.metaheuristic.de;

import java.util.List;
import org.optiprot.metaheuristic.de.DESolution;

/**
 * interface for mutation operator of DE
 *
 * @author victor
 */
public interface IDEMutation {

    public void mutate( DESolution mutant, int idx, List<DESolution> population, double F);

}
