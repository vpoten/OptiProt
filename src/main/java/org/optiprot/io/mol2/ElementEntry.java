/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io.mol2;

/**
 * Columns represent:
 *   - atomic number (used as an index to the vector as well)
 *   - elemental symbol
 *   - Allred and Rochow electronegativity  0.0 if unknown
 *   - covalent radii (in Angstrom)         1.6 if unknown
 *       from http://dx.doi.org/10.1039/b801115j
 *   - "bond order" radii -- ignored, but included for compatibility
 *   - van der Waals radii (in Angstrom)    2.0 if unknown
 *       from http://dx.doi.org/10.1021/jp8111556
 *   - maximum bond valence                   6 if unknown
 *   - IUPAC recommended atomic masses (in amu)
 *   - Pauling electronegativity            0.0 if unknown
 *   - ionization potential (in eV)         0.0 if unknown
 *   - electron affinity (in eV)            0.0 if unknown
 *   - RGB values (defaults for visualization)
 *   - element name (in English)
 *
 *
 * @author victor
 */

public class ElementEntry {
    private int Num;
    private String Symb;
    private double ARENeg;
    private double RCov;
    private double RBO;
    private double RVdW;
    private int MaxBnd;
    private double Mass;
    private double ElNeg;
    private double Ionization;
    private double ElAffinity;
    private double Red;
    private double Green;
    private double Blue;
    private String Name;

    public ElementEntry(int Num, String Symb, double ARENeg, double RCov,
            double RBO, double RVdW, int MaxBnd, double Mass, double ElNeg,
            double Ionization, double ElAffinity,
            double Red, double Green, double Blue, String Name) {
        this.Num = Num;
        this.Symb = Symb;
        this.ARENeg = ARENeg;
        this.RCov = RCov;
        this.RBO = RBO;
        this.RVdW = RVdW;
        this.MaxBnd = MaxBnd;
        this.Mass = Mass;
        this.ElNeg = ElNeg;
        this.Ionization = Ionization;
        this.ElAffinity = ElAffinity;
        this.Red = Red;
        this.Green = Green;
        this.Blue = Blue;
        this.Name = Name;
    }

    /**
     * @return the Num
     */
    public int getNum() {
        return Num;
    }

    /**
     * @param Num the Num to set
     */
    public void setNum(int Num) {
        this.Num = Num;
    }

    /**
     * @return the Symb
     */
    public String getSymb() {
        return Symb;
    }

    /**
     * @param Symb the Symb to set
     */
    public void setSymb(String Symb) {
        this.Symb = Symb;
    }

    /**
     * @return the ARENeg
     */
    public double getARENeg() {
        return ARENeg;
    }

    /**
     * @param ARENeg the ARENeg to set
     */
    public void setARENeg(double ARENeg) {
        this.ARENeg = ARENeg;
    }

    /**
     * @return the RCov
     */
    public double getRCov() {
        return RCov;
    }

    /**
     * @param RCov the RCov to set
     */
    public void setRCov(double RCov) {
        this.RCov = RCov;
    }

    /**
     * @return the RBO
     */
    public double getRBO() {
        return RBO;
    }

    /**
     * @param RBO the RBO to set
     */
    public void setRBO(double RBO) {
        this.RBO = RBO;
    }

    /**
     * @return the RVdW
     */
    public double getRVdW() {
        return RVdW;
    }

    /**
     * @param RVdW the RVdW to set
     */
    public void setRVdW(double RVdW) {
        this.RVdW = RVdW;
    }

    /**
     * @return the MaxBnd
     */
    public int getMaxBnd() {
        return MaxBnd;
    }

    /**
     * @param MaxBnd the MaxBnd to set
     */
    public void setMaxBnd(int MaxBnd) {
        this.MaxBnd = MaxBnd;
    }

    /**
     * @return the Mass
     */
    public double getMass() {
        return Mass;
    }

    /**
     * @param Mass the Mass to set
     */
    public void setMass(double Mass) {
        this.Mass = Mass;
    }

    /**
     * @return the ElNeg
     */
    public double getElNeg() {
        return ElNeg;
    }

    /**
     * @param ElNeg the ElNeg to set
     */
    public void setElNeg(double ElNeg) {
        this.ElNeg = ElNeg;
    }

    /**
     * @return the Ionization
     */
    public double getIonization() {
        return Ionization;
    }

    /**
     * @param Ionization the Ionization to set
     */
    public void setIonization(double Ionization) {
        this.Ionization = Ionization;
    }

    /**
     * @return the ElAffinity
     */
    public double getElAffinity() {
        return ElAffinity;
    }

    /**
     * @param ElAffinity the ElAffinity to set
     */
    public void setElAffinity(double ElAffinity) {
        this.ElAffinity = ElAffinity;
    }

    /**
     * @return the Red
     */
    public double getRed() {
        return Red;
    }

    /**
     * @param Red the Red to set
     */
    public void setRed(double Red) {
        this.Red = Red;
    }

    /**
     * @return the Green
     */
    public double getGreen() {
        return Green;
    }

    /**
     * @param Green the Green to set
     */
    public void setGreen(double Green) {
        this.Green = Green;
    }

    /**
     * @return the Blue
     */
    public double getBlue() {
        return Blue;
    }

    /**
     * @param Blue the Blue to set
     */
    public void setBlue(double Blue) {
        this.Blue = Blue;
    }

    /**
     * @return the Name
     */
    public String getName() {
        return Name;
    }

    /**
     * @param Name the Name to set
     */
    public void setName(String Name) {
        this.Name = Name;
    }


}
