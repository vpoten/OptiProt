/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.aacids;

/**
 * Interface for amino acid genes
 *
 * @author victor
 */
public interface IAAGene {

     public IAABasic getAminoAcid();
     public void setAminoAcid(IAABasic aa);
}
