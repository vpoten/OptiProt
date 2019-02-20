/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.aacids;

import org.biojava.bio.structure.AminoAcid;
import org.optiprot.OptiProtParameters;
import org.biojava.bio.structure.Atom;

/**
 * Interface for amino acid alleles
 *
 * @author victor
 */
public interface IAABasic {


    public void applyMutation(int i, double v);

    public void applyMutationChangeRotamer();

    public Atom[] getAtoms();

    public int getNumAtoms();

    public Atom getVectorNCa();

    public Atom getVectorCaC();

    public Atom getCa();
    
    public Atom getCb();

    public Atom getN();

    public Atom getC();

    public Atom getAtom( String name );

    public OptiProtParameters getParameters();

    public void setParameters(OptiProtParameters parameters);

    public String getName();

    public Object clone();

    public AminoAcid toAminoAcid(boolean generatesH);

    public Integer getRotIdx();

    public void changeRotIdx(Integer rotIdx);
    
}
