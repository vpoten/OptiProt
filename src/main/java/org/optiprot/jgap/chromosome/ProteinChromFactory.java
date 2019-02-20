/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.jgap.chromosome;

import org.optiprot.jgap.genes.AABasicGene;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.biojava.bio.structure.*;
import org.jgap.*;
import org.optiprot.jgap.OptiProtConfiguration;
import org.optiprot.aacids.*;

/**
 *
 * @author victor
 */
public class ProteinChromFactory {


    /**
     * creates a chromosome for a given protein structure
     *
     * @param a_conf
     * @param a_protChain
     * @return
     * @throws org.jgap.InvalidConfigurationException
     * @throws org.biojava.bio.structure.StructureException
     */
    public static IChromosome create( final OptiProtConfiguration a_conf,
            final Chain a_protChain )
            throws InvalidConfigurationException, StructureException{

        ProteinChromosome chrom=new ProteinChromosome(a_conf);

        chrom.setGenes( createGenes(a_conf, a_protChain) );

        return chrom;
    }


    /**
     *
     * @param a_conf
     * @param aacid
     * @return
     * @throws org.jgap.InvalidConfigurationException
     * @throws org.biojava.bio.structure.StructureException
     */
    private static Gene createAAGen( OptiProtConfiguration a_conf, Group aacid)
            throws InvalidConfigurationException, StructureException {

        IAAGene gene=new AABasicGene(a_conf);

        gene.setAminoAcid( new AABasicCaCbNC(aacid, a_conf.getParameters() ) );

        return (Gene)gene;
    }


    /**
     * 
     * @param a_conf
     * @param a_protChain
     * @return
     * @throws org.jgap.InvalidConfigurationException
     * @throws org.biojava.bio.structure.StructureException
     */
    private static Gene[] createGenes( OptiProtConfiguration a_conf,
            Chain a_protChain )
            throws InvalidConfigurationException, StructureException {

        Iterator<Group> aacids=a_protChain.getAtomGroups().iterator();
        
        Gene [] genes=new Gene [a_protChain.getAtomGroups().size()];

        int i=0;

        while( aacids.hasNext() ){

            genes[i++]=createAAGen( a_conf, aacids.next() );
        }

        return genes;
    }

    /**
     * converts the given chromosome to structure chain
     *
     * @param chrom
     * @return
     */
    public static Chain toChain( IChromosome chrom, boolean generatesH ){

        Chain chain=new ChainImpl();

        Gene [] genes=chrom.getGenes();

        List<Group> list = new ArrayList<Group>();

        int contAtoms=1;
        
        for(int i=0;i<genes.length;i++){

            AminoAcid aa=((IAABasic)genes[i].getAllele()).toAminoAcid(generatesH);

            aa.setPDBCode( Integer.toString(i+1) );

            for( Atom atom : aa.getAtoms() ){
                atom.setPDBserial(contAtoms++);
            }

            list.add( aa  );
        }

        chain.setAtomGroups(list);
        return chain;
    }

    
    
}
