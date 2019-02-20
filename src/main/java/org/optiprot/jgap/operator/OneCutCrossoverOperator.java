/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.operator;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.jgap.*;
import org.jgap.impl.*;
import org.optiprot.jgap.chromosome.ProteinChromosome;

/**
 *
 * @author victor
 */
public class OneCutCrossoverOperator extends CrossoverOperator {

    public OneCutCrossoverOperator(final Configuration a_config )
        throws InvalidConfigurationException {

        super(a_config);
    }

    public OneCutCrossoverOperator(final Configuration a_config,
            int a_crossoverRate )
        throws InvalidConfigurationException {

        super(a_config, a_crossoverRate);
    }

//    @Override
//    public void operate(final Population a_population,
//                      final List a_candidateChromosomes) {
//
//    }

    @Override
    protected void doCrossover(IChromosome firstMate,
                           IChromosome secondMate,
                           java.util.List a_candidateChromosomes,
                           RandomGenerator generator){

        if( firstMate instanceof ProteinChromosome ){

            /* Find the cut place */
            int cutPoint = generator.nextInt( firstMate.size() - 1);

            // take a copy of both mates
            IChromosome copyOfMate1 = (IChromosome) ((ProteinChromosome )firstMate).clone();
            IChromosome copyOfMate2 = (IChromosome) ((ProteinChromosome )secondMate).clone();
            
            try {
                ((ProteinChromosome) copyOfMate1).oneCutCross(cutPoint,
                        (ProteinChromosome) secondMate);
                ((ProteinChromosome) copyOfMate2).oneCutCross(cutPoint,
                    (ProteinChromosome)firstMate);
            } catch (Exception ex) {
                Logger.getLogger(OneCutCrossoverOperator.class.getName()).log(Level.SEVERE, null, ex);
            }

            
            // add the offspring to the list of candidates
            a_candidateChromosomes.add(copyOfMate1);
            a_candidateChromosomes.add(copyOfMate2);
        }
        else{
              throw new ClassCastException();
        }
    }
    
}
