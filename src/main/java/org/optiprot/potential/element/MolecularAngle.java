/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.potential.element;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;

/**
 *
 * @author victor
 */
public class MolecularAngle extends MolecularTorsional {

    double m_distance=-1;

    @Override
    public double getAngle() {

        if( m_angle>0 )
            return m_angle;
        
        Atom at1=this.getAtoms()[0];
        Atom at2=this.getAtoms()[1];
        Atom at3=this.getAtoms()[2];
        
        try {
            m_angle=Calc.angle(Calc.substract(at1, at2), Calc.substract(at3, at2));
            return m_angle;
        } catch (StructureException ex) {
            Logger.getLogger(MolecularAngle.class.getName()).log(Level.SEVERE, null, ex);
            return 0;
        }
    }

    public double getDistance() {

        if( m_distance>0 )
            return m_distance;

        Atom at1=this.getAtoms()[0];
        Atom at3=this.getAtoms()[2];
        try {
            m_distance=Calc.getDistance(at1, at3);
            return m_distance;
        } catch (StructureException ex) {
            Logger.getLogger(MolecularAngle.class.getName()).log(Level.SEVERE, null, ex);
            return 0;
        }
    }

}
