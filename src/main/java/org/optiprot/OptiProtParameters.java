/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot;

import java.util.ArrayList;
import java.util.List;
import org.optiprot.io.CharmmTopParReader;
import org.optiprot.maths.AtomGrid;
import org.optiprot.potential.CharmmForceField;
import org.optiprot.potential.IForceField;
import org.optiprot.potential.element.CharmmElement;
import org.optiprot.potential.element.CharmmResidue;
import org.optiprot.rotamer.RotamerLibrary;

/**
 * Encapsulates the constants used by the program
 *
 * @author victor
 */
public class OptiProtParameters implements java.io.Serializable {

    private double MutDespUnit=0.6d;

    private double MutAngleStepMax=Math.PI/10.0;

    private int MutRateAngle=2;//probability = 1/rate
    private int MutRateRotamer=8;
    private int CrossOneCutRate=4;

    private String TinkerPath="";
    private String TinkerForceField="";

    private String StructureWriteDir="";

    private RotamerLibrary rotLib=null;

    private int populationSize=20;

    private int numberOfEvolutions=100;

    private int resultSize=1;

    private String name="";

    private String workDir="";

    private boolean generatesH=false;

    private IForceField forceField=null;

    private List<CharmmResidue> topology=null;

    private String PDBDir="";

    //parameters for potential energy calculation
    private boolean calcMMechanics=true;
    private boolean calcGBorn=true;
    private boolean calcSASA=true;

    private AtomGrid grid=new AtomGrid(100,100,100);

    public OptiProtParameters(){

    }


    /**
     * @return the MutDespUnit
     */
    public double getMutDespUnit() {
        return MutDespUnit;
    }

    /**
     * @param MutDespUnit the MutDespUnit to set
     */
    public void setMutDespUnit(double MutDespUnit) {
        this.MutDespUnit = MutDespUnit;
    }

    
    public List<CharmmResidue> getTopology() {
        return topology;
    }

    public void setTopology( List<CharmmResidue> list ) {
        topology=list;
    }


    /**
     * @return the TinkerPath
     */
    public String getTinkerPath() {
        return TinkerPath;
    }

    /**
     * @param TinkerPath the TinkerPath to set
     */
    public void setTinkerPath(String TinkerPath) {
        this.TinkerPath = TinkerPath;
    }

    /**
     * @return the TinkerForceField
     */
    public String getTinkerForceField() {
        return TinkerForceField;
    }

    /**
     * @param TinkerForceField the TinkerForceField to set
     */
    public void setTinkerForceField(String ForceField) {
        this.TinkerForceField = ForceField;
    }

    /**
     * @return the StructureWriteDir
     */
    public String getStructureWriteDir() {
        return StructureWriteDir;
    }

    /**
     * @param StructureWriteDir the StructureWriteDir to set
     */
    public void setStructureWriteDir(String StructureWriteDir) {
        this.StructureWriteDir = StructureWriteDir;
    }

    /**
     * @return the rotLib
     */
    public RotamerLibrary getRotLib() {
        return rotLib;
    }

    /**
     * @param rotLib the rotLib to set
     */
    public void setRotLib(RotamerLibrary rotLib) {
        this.rotLib = rotLib;
    }

    /**
     * @return the MutAngleStepMax
     */
    public double getMutAngleStepMax() {
        return MutAngleStepMax;
    }

    /**
     * @param MutAngleStepMax the MutAngleStepMax to set
     */
    public void setMutAngleStepMax(double MutAngleStepMax) {
        this.MutAngleStepMax = MutAngleStepMax;
    }

    /**
     * @return the populationSize
     */
    public int getPopulationSize() {
        return populationSize;
    }

    /**
     * @param populationSize the populationSize to set
     */
    public void setPopulationSize(int populationSize) {
        this.populationSize = populationSize;
    }

    /**
     * @return the MutRateAngle
     */
    public int getMutRateAngle() {
        return MutRateAngle;
    }

    /**
     * @param MutRateAngle the MutRateAngle to set
     */
    public void setMutRateAngle(int MutRateAngle) {
        this.MutRateAngle = MutRateAngle;
    }

    /**
     * @return the MutRateRotamer
     */
    public int getMutRateRotamer() {
        return MutRateRotamer;
    }

    /**
     * @param MutRateRotamer the MutRateRotamer to set
     */
    public void setMutRateRotamer(int MutRateRotamer) {
        this.MutRateRotamer = MutRateRotamer;
    }

    /**
     * @return the CrossOneCutRate
     */
    public int getCrossOneCutRate() {
        return CrossOneCutRate;
    }

    /**
     * @param CrossOneCutRate the CrossOneCutRate to set
     */
    public void setCrossOneCutRate(int CrossOneCutRate) {
        this.CrossOneCutRate = CrossOneCutRate;
    }

    /**
     * @return the numberOfEvolutions
     */
    public int getNumberOfEvolutions() {
        return numberOfEvolutions;
    }

    /**
     * @param numberOfEvolutions the numberOfEvolutions to set
     */
    public void setNumberOfEvolutions(int numberOfEvolutions) {
        this.numberOfEvolutions = numberOfEvolutions;
    }

    /**
     * @return the resultSize
     */
    public int getResultSize() {
        return resultSize;
    }

    /**
     * @param resultSize the resultSize to set
     */
    public void setResultSize(int resultSize) {
        this.resultSize = resultSize;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the workDir
     */
    public String getWorkDir() {
        return workDir;
    }

    /**
     * @param workDir the workDir to set
     */
    public void setWorkDir(String workDir) {
        this.workDir = workDir;
    }

    /**
     * @return the generatesH
     */
    public boolean isGeneratesH() {
        return generatesH;
    }

    /**
     * @param generatesH the generatesH to set
     */
    public void setGeneratesH(boolean generatesH) {
        this.generatesH = generatesH;
    }

    /**
     * @return the forceField
     */
    public IForceField getForceField() {
        return forceField;
    }

    /**
     * @param forceField the forceField to set
     */
    public void setForceField(IForceField forceField) {
        this.forceField = forceField;
    }

    /**
     * @return the calcMMechanics
     */
    public boolean isCalcMMechanics() {
        return calcMMechanics;
    }

    /**
     * @param calcMMechanics the calcMMechanics to set
     */
    public void setCalcMMechanics(boolean calcMMechanics) {
        this.calcMMechanics = calcMMechanics;
    }

    /**
     * @return the calcGBorn
     */
    public boolean isCalcGBorn() {
        return calcGBorn;
    }

    /**
     * @param calcGBorn the calcGBorn to set
     */
    public void setCalcGBorn(boolean calcGBorn) {
        this.calcGBorn = calcGBorn;
    }

    /**
     * @return the calcSASA
     */
    public boolean isCalcSASA() {
        return calcSASA;
    }

    /**
     * @param calcSASA the calcSASA to set
     */
    public void setCalcSASA(boolean calcSASA) {
        this.calcSASA = calcSASA;
    }

    /**
     * @return the grid
     */
    public AtomGrid getGrid() {
        return grid;
    }

    /**
     * create a parameters object and loads force field and rotamer library
     *
     * @param workdir : working directory
     * @param rotdir : rotamer library directory
     * @return
     * @throws java.lang.Exception
     */
    public static OptiProtParameters createParameters(String workdir, String rotdir )
            throws Exception{

        OptiProtParameters parameters=new OptiProtParameters();

        parameters.setWorkDir( workdir );

        //load topology
        List<CharmmResidue> topology=null;
        ArrayList<CharmmElement> elements=new ArrayList<CharmmElement>();

        topology = CharmmTopParReader.parseTopFile(
                "org/optiprot/data/top_all22_prot.inp",
                elements);

        parameters.setTopology(topology);

        //load force field
        CharmmForceField ffield=new CharmmForceField(elements);


        CharmmTopParReader.parseParFile(
                "org/optiprot/data/par_all22_prot.inp",
                ffield, elements);

        parameters.setForceField(ffield);


        //create and load rotamer library
        RotamerLibrary lib=new RotamerLibrary();

        lib.loadLibrary( rotdir, "allH");

        parameters.setRotLib(lib);

        return parameters;

    }

    /**
     * @return the PDBDir
     */
    public String getPDBDir() {
        return PDBDir;
    }

    /**
     * @param PDBDir the PDBDir to set
     */
    public void setPDBDir(String PDBDir) {
        this.PDBDir = PDBDir;
    }
    
}
