/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.operator;

import java.util.List;
import java.util.Vector;
import org.optiprot.jgap.chromosome.ProteinChromosome;
import org.jgap.*;
import org.jgap.impl.*;

/**
 *
 * @author victor
 */
public class AngleMutationOperator extends MutationOperator {

    public AngleMutationOperator(final Configuration a_config)
        throws InvalidConfigurationException {

        super(a_config);
    }

    public AngleMutationOperator(final Configuration a_config,
                           final int a_desiredMutationRate)
        throws InvalidConfigurationException {

        super(a_config, a_desiredMutationRate);
    }

    @Override
    public void operate(final Population a_population,
                      final List a_candidateChromosomes) {

        if (a_population == null || a_candidateChromosomes == null) {
          // Population or candidate chromosomes list empty:
          // nothing to do.
          // -----------------------------------------------
          return;
        }
        if (getMutationRate() == 0 && getMutationRateCalc() == null) {
          // If the mutation rate is set to zero and dynamic mutation rate is
          // disabled, then we don't perform any mutation.
          // ----------------------------------------------------------------
          return;
        }

        Configuration conf = getConfiguration();

        // Determine the mutation rate. If dynamic rate is enabled, then
        // calculate it using the IUniversalRateCalculator instance.
        // Otherwise, go with the mutation rate set upon construction.
        // -------------------------------------------------------------
        boolean mutate = false;
        RandomGenerator generator = conf.getRandomGenerator();

        // It would be inefficient to create copies of each Chromosome just
        // to decide whether to mutate them. Instead, we only make a copy
        // once we've positively decided to perform a mutation.
        // ----------------------------------------------------------------
        int size = Math.min(conf.getPopulationSize(),
                            a_population.size());
        IGeneticOperatorConstraint constraint = conf.getJGAPFactory().
            getGeneticOperatorConstraint();

        for (int i = 0; i < size; i++) {
          IChromosome chrom = a_population.getChromosome(i);
          Gene[] genes = chrom.getGenes();
          IChromosome copyOfChromosome = null;
          // For each Chromosome in the population...
          // ----------------------------------------

          /* Find a aminoacid to mutate */
          int target = generator.nextInt( chrom.size() - 1);

          if (getMutationRateCalc() != null) {
            // If it's a dynamic mutation rate then let the calculator decide
            // whether the current gene should be mutated.
            // --------------------------------------------------------------
            mutate = getMutationRateCalc().toBePermutated(chrom, target);
          }
          else {
            // Non-dynamic, so just mutate based on the the current rate.
            // In fact we use a rate of 1/m_mutationRate.
            // ----------------------------------------------------------
            mutate = (generator.nextInt(getMutationRate()) == 0);
          }

          if (mutate) {
            // Verify that crossover allowed.
            // ------------------------------
            
            if (constraint != null) {
              List v = new Vector();
              v.add(chrom);
              if (!constraint.isValid(a_population, v, this)) {
                continue;
              }
            }
            // Now that we want to actually modify the Chromosome,
            // let's make a copy of it (if we haven't already) and
            // add it to the candidate chromosomes so that it will
            // be considered for natural selection during the next
            // phase of evolution. Then we'll set the gene's value
            // to a random value as the implementation of our
            // "mutation" of the gene.
            // ---------------------------------------------------
            if (copyOfChromosome == null) {
              // ...take a copy of it...
              // -----------------------
              copyOfChromosome = (IChromosome) ((ProteinChromosome )chrom).clone();
              // ...add it to the candidate pool...
              // ----------------------------------
              a_candidateChromosomes.add(copyOfChromosome);
              // ...then mutate all its genes...
              // -------------------------------
              genes = copyOfChromosome.getGenes();
            }
            
            mutateGene(copyOfChromosome, target, genes[target], generator);
            
          }
        }
      }

      private void mutateGene(final IChromosome chromosome, int target, final Gene a_gene,
              final RandomGenerator a_generator)  throws ClassCastException {
          
          if( chromosome instanceof ProteinChromosome ){
          
              ProteinChromosome pChro=(ProteinChromosome)chromosome;
              
              double percentage = -1 + a_generator.nextDouble()*2;

              pChro.changePhi(target, percentage);
              
              percentage = -1 + a_generator.nextDouble()*2;

              pChro.changePsi(target, percentage);
          }
          else{
              throw new ClassCastException();
          }
          
      }

       
    

}
