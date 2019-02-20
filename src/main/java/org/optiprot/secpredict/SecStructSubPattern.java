/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.secpredict;

/**
 *
 * @author victor
 */
public class SecStructSubPattern {
    
    int aminoacid;//aa index
    int type;//type of secondary structure
    int subtype;//subtype of secondary structure
    double phi;
    double psi;

    public SecStructSubPattern(int aminoacid, int type, int subtype, double phi, double psi) {
        this.aminoacid = aminoacid;
        this.type = type;
        this.subtype = subtype;
        this.phi = phi;
        this.psi = psi;
    }

    @Override
    public SecStructSubPattern clone(){
        return new SecStructSubPattern(this.aminoacid, this.type, this.subtype,
                this.phi, this.psi);
    }
}
