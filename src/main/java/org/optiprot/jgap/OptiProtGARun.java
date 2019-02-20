/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap;

import org.optiprot.jgap.chromosome.ProteinChromFactory;
import org.optiprot.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.jgap.*;


/**
 *
 * @author victor
 */
public class OptiProtGARun {

    static public List<Chain> run( final List<Chain> p_chains,
            final OptiProtParameters p_parameters )
            throws InvalidConfigurationException, StructureException {

        OptiProtConfiguration conf=new OptiProtConfiguration(p_parameters);

        Iterator<Chain> it=p_chains.iterator();

        List<IChromosome> listChrom=new ArrayList<IChromosome>();

        while(it.hasNext()){

            listChrom.add( ProteinChromFactory.create(conf, it.next()) );
        }

        Population pop=new Population(
                conf,
                (IChromosome [])listChrom.toArray(new IChromosome [listChrom.size()])
                );

        Genotype geno=new Genotype(conf, pop);

        geno.evolve( p_parameters.getNumberOfEvolutions() );

        listChrom=(List<IChromosome>)geno.getFittestChromosomes( p_parameters.getResultSize() );

        return toChain(listChrom);

    }

    /**
     * converts a list of chromosomes to a list of Chain
     *
     * @param listChrom
     * @return
     */
    private static List<Chain> toChain(List<IChromosome> listChrom) {

        List<Chain> chains=new ArrayList<Chain>();

        Iterator<IChromosome> it=listChrom.iterator();

        while( it.hasNext() ){

            chains.add( ProteinChromFactory.toChain(it.next(),true) );
        }

        return chains;
    }

}
