/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import java.util.*;
import org.biojava.bio.structure.*;
import org.optiprot.potential.IForceField;

/**
 *
 * @author victor
 */
public class BSPTree {

    private BSPTree m_nodeLeft=null;
    private BSPTree m_nodeRight=null;
    private List<Atom> m_elements=null;
    private double m_median=0;
    private int m_axis = -1;
    private int m_depth=0;
    private int numAtoms=0;

    private Atom [] m_bBox=null;

    static private int STOP_LIST_SIZE=10;
    static private double MAX_ATOM_RADII=2.0;

    
    /**
     * class for results of nearestNeighbour method
     */
    class NNInfo {
        public Atom atom=null;
        public double sqDist=java.lang.Double.MAX_VALUE;

        private NNInfo() {
        }

        private NNInfo(Atom at, double sqdist) {
            atom=at;
            sqDist=sqdist;
        }
    }

    //////////////////////////

    public BSPTree() {
    }

    public BSPTree(Chain chain) {
        List<Atom> listAtoms = getListAtoms(chain);
        buildBSPTree(listAtoms);
    }

    public BSPTree(List<Atom> listAtoms) {
        buildBSPTree(listAtoms);
    }

    private void buildBSPTree(List<Atom> listAtoms) {

        numAtoms=listAtoms.size();

        Atom [] bbox = boundingBox( listAtoms );

        setBBox(bbox);
       
        BSPTree currentNode=null;
        ArrayDeque<BSPTree> fifo=new ArrayDeque<BSPTree>();
        ArrayDeque<List<Atom>> fifoList=new ArrayDeque<List<Atom>>();
        fifo.add(this);
        fifoList.add(listAtoms);

        while( !fifo.isEmpty() ){

            currentNode=fifo.poll();
            listAtoms=fifoList.poll();

            if( listAtoms.size()<=BSPTree.STOP_LIST_SIZE ){
                currentNode.setElements(listAtoms);
                continue;
            }

            int axis = 0;
            bbox=currentNode.getBBox();

            //use the more spread axis

            double spread=bbox[1].getX()-bbox[0].getX();

            if(spread<bbox[1].getY()-bbox[0].getY()){
                axis=1;
                spread=bbox[1].getY()-bbox[0].getY();
            }

            if(spread<bbox[1].getZ()-bbox[0].getZ()){
                axis=2;
                spread=bbox[1].getZ()-bbox[0].getZ();
            }

            List<Atom> listBefore=new ArrayList<Atom>();
            List<Atom> listAfter=new ArrayList<Atom>();

            double l_median=0;

            while(true){

                currentNode.setAxis(axis);

                listAfter.clear();
                listBefore.clear();

                //calculates the median at chosen axis
                l_median=0;

                for( Atom atom : listAtoms ){
                    l_median+=atom.getCoords()[axis];
                }

                l_median/=listAtoms.size();

                currentNode.setMedian(l_median);

                //divide the points after and before the median
                for( Atom atom : listAtoms ){

                    if( atom.getCoords()[axis]>l_median ){
                        listAfter.add(atom);
                    }
                    else{
                        listBefore.add(atom);
                    }
                }

                if( listBefore.size()==listAtoms.size() ||
                        listAfter.size()==listAtoms.size() ){
                    //if cannot separate the points choose another axis
                    axis=(axis+1)%3;
                }
                else{
                    break;
                }
            }

            if( !listBefore.isEmpty() ){

                Atom [] bbox2=new Atom [2];
                bbox2[0]=new AtomImpl();
                bbox2[1]=new AtomImpl();

                bbox2[0].setCoords( bbox[0].getCoords().clone() );
                bbox2[1].setCoords( bbox[1].getCoords().clone() );
                bbox2[1].getCoords()[axis]=l_median;

                currentNode.setNodeLeft( new BSPTree() );
                currentNode.getNodeLeft().setDepth( currentNode.getDepth()+1 );
                currentNode.getNodeLeft().setBBox( bbox2 );

                fifo.add( currentNode.getNodeLeft() );
                fifoList.add( listBefore );

            }

            if( !listAfter.isEmpty() ){

                Atom [] bbox1=new Atom [2];
                bbox1[0]=new AtomImpl();
                bbox1[1]=new AtomImpl();

                bbox1[0].setCoords( bbox[0].getCoords().clone() );
                bbox1[1].setCoords( bbox[1].getCoords().clone() );
                bbox1[0].getCoords()[axis]=l_median;

                currentNode.setNodeRight( new BSPTree() );
                currentNode.getNodeRight().setDepth( currentNode.getDepth()+1 );
                currentNode.getNodeRight().setBBox( bbox1 );

                fifo.add( currentNode.getNodeRight() );
                fifoList.add( listAfter );

            }

            listAtoms=null;
        }
        
        fifo=null;
        fifoList=null;
        listAtoms=null;
    }

    static public List<Atom> getListAtoms(Chain chain){

        ArrayList<Atom> list=new ArrayList<Atom>();

        for ( Group group : chain.getAtomGroups() ){
            for( Atom atom : group.getAtoms() ){
                list.add(atom);
            }
        }

        return list;
    }

    /**
     * @return the numAtoms
     */
    public int getNumAtoms() {
        return numAtoms;
    }

    /**
     * calcs the bounding box
     * 
     * @param chain
     * @return array of lenght 2. array[0]=min corner array[1]=max corner
     */
    static public Atom [] boundingBox( Chain chain ){

        return boundingBox(getListAtoms(chain));
    }

    static public Atom [] boundingBox( List<Atom> list ){

        Atom [] bbox=new Atom [2];
        bbox[0]=new AtomImpl();
        bbox[1]=new AtomImpl();
        bbox[0].setX(java.lang.Double.MAX_VALUE);
        bbox[0].setY(java.lang.Double.MAX_VALUE);
        bbox[0].setZ(java.lang.Double.MAX_VALUE);
        bbox[1].setX(java.lang.Double.MIN_VALUE);
        bbox[1].setY(java.lang.Double.MIN_VALUE);
        bbox[1].setZ(java.lang.Double.MIN_VALUE);

        for( Atom atom : list ){

            bbox[1].setX( Math.max(bbox[1].getX(), atom.getX()) );
            bbox[1].setY( Math.max(bbox[1].getY(), atom.getY()) );
            bbox[1].setZ( Math.max(bbox[1].getZ(), atom.getZ()) );

            bbox[0].setX( Math.min(bbox[0].getX(), atom.getX()) );
            bbox[0].setY( Math.min(bbox[0].getY(), atom.getY()) );
            bbox[0].setZ( Math.min(bbox[0].getZ(), atom.getZ()) );
        }

        return bbox;
    }

////    private void createBSPTree(List<Atom> listAtoms, Atom[] bbox,int depth) {
////
////        this.setDepth(depth);
////
////        int axis = 0;
////        double spread=bbox[1].getX()-bbox[0].getX();
////
////        if(spread<bbox[1].getY()-bbox[0].getY()){
////            axis=1;
////            spread=bbox[1].getY()-bbox[0].getY();
////        }
////
////        if(spread<bbox[1].getZ()-bbox[0].getZ()){
////            axis=2;
////            spread=bbox[1].getZ()-bbox[0].getZ();
////        }
////
////        this.setAxis(axis);
////
////
////        if( listAtoms.size()<=BSPTree.STOP_LIST_SIZE ){
////            setElements(listAtoms);
////            return;
////        }
////
////        double l_median=0;
////
////        for( Atom atom : listAtoms ){
////            l_median+=atom.getCoords()[axis];
////        }
////
////        l_median/=listAtoms.size();
////
////        this.setMedian(l_median);
////
////        List<Atom> listBefore=new ArrayList<Atom>();
////        List<Atom> listAfter=new ArrayList<Atom>();
////
////        for( Atom atom : listAtoms ){
////
////            if( atom.getCoords()[axis]>=l_median ){
////                listAfter.add(atom);
////            }
////            if( atom.getCoords()[axis]<l_median ){
////                listBefore.add(atom);
////            }
////        }
////
////        Atom [] bbox1=new Atom [2];
////        Atom [] bbox2=new Atom [2];
////        bbox1[0]=new AtomImpl();
////        bbox1[1]=new AtomImpl();
////        bbox2[0]=new AtomImpl();
////        bbox2[1]=new AtomImpl();
////
////        bbox1[0].setCoords( bbox[0].getCoords().clone() );
////        bbox1[1].setCoords( bbox[1].getCoords().clone() );
////        bbox1[0].getCoords()[axis]=l_median;
////
////        bbox2[0].setCoords( bbox[0].getCoords().clone() );
////        bbox2[1].setCoords( bbox[1].getCoords().clone() );
////        bbox2[1].getCoords()[axis]=l_median;
////
////        setNodeLeft( new BSPTree() );
////        setNodeRight( new BSPTree() );
////
////        getNodeRight().createBSPTree(listAfter, bbox1, depth+1);
////        getNodeLeft().createBSPTree(listBefore, bbox2, depth+1);
////
////    }

    /**
     * @return the m_nodeLeft
     */
    private BSPTree getNodeLeft() {
        return m_nodeLeft;
    }

    /**
     * @param m_nodeLeft the m_nodeLeft to set
     */
    private void setNodeLeft(BSPTree m_nodeLeft) {
        this.m_nodeLeft = m_nodeLeft;
    }

    /**
     * @return the m_nodeRight
     */
    private BSPTree getNodeRight() {
        return m_nodeRight;
    }

    /**
     * @param m_nodeRight the m_nodeRight to set
     */
    private void setNodeRight(BSPTree m_nodeRight) {
        this.m_nodeRight = m_nodeRight;
    }

    /**
     * @return the m_elements
     */
    private List<Atom> getElements() {
        return m_elements;
    }

    /**
     * @param m_elements the m_elements to set
     */
    private void setElements(List<Atom> m_elements) {
        this.m_elements = m_elements;
    }

    /**
     * @return the m_median
     */
    private double getMedian() {
        return m_median;
    }

     /**
     * @return the m_bBox
     */
    public Atom[] getBBox() {
        return m_bBox;
    }

    /**
     * @param m_bBox the m_bBox to set
     */
    private void setBBox(Atom[] m_bBox) {
        this.m_bBox = m_bBox;
    }

    /**
     * @return the m_axis
     */
    private int getAxis() {
        return m_axis;
    }

    /**
     * @param m_axis the m_axis to set
     */
    private void setAxis(int m_axis) {
        this.m_axis = m_axis;
    }


    /**
     * @param m_median the m_median to set
     */
    private void setMedian(double m_median) {
        this.m_median = m_median;
    }

    /**
     * @return the m_depth
     */
    private int getDepth() {
        return m_depth;
    }

    /**
     * @param m_depth the m_depth to set
     */
    private void setDepth(int m_depth) {
        this.m_depth = m_depth;
    }

    /**
     *
     * @param atom
     * @param ffield : implements getVDWRadius
     * @param p_radius : probe radii
     * @return
     */
    public boolean isInsideVolume(final Atom atom, IForceField ffield, final double p_radius){

////        int axis=getAxis();
////
////        if( getElements()!=null ){
////            //if is a leaf node do the test
////
////            for(Atom element : getElements() ){
////
////                double radius=ffield.getVDWRadius(element)+p_radius;
////
////                if ( CalcGeom.squareDistance(atom,element) <=
////                        radius*radius ) {
////                    return true;
////                }
////
////            }
////
////            return false;
////        }
////
////        if( atom.getCoords()[axis]<getMedian()-p_radius-MAX_ATOM_RADII ){
////            return getNodeLeft().isInsideVolume(atom,ffield,p_radius);
////        }
////        else if( atom.getCoords()[axis]>getMedian()+p_radius+MAX_ATOM_RADII ){
////            return getNodeRight().isInsideVolume(atom,ffield,p_radius);
////        }
////        else{
////            return getNodeLeft().isInsideVolume(atom,ffield,p_radius) ||
////                    getNodeRight().isInsideVolume(atom,ffield,p_radius);
////        }

        BSPTree currentNode=null;
        ArrayDeque<BSPTree> fifo=new ArrayDeque<BSPTree>();
        fifo.add(this);

        while( !fifo.isEmpty() ){

            currentNode=fifo.poll();

            int axis=currentNode.getAxis();

            if( currentNode.getElements()!=null ){
            //if is a leaf node do the test

                for(Atom element : currentNode.getElements() ){

                    double radius=ffield.getVDWRadius(element)+p_radius;

                    if ( CalcGeom.squareDistance(atom,element) <=
                            radius*radius ) {
                        return true;
                    }

                }
            }
            else if( atom.getCoords()[axis]<currentNode.getMedian()-p_radius-MAX_ATOM_RADII ){
                if( currentNode.getNodeLeft()!=null )
                    fifo.add( currentNode.getNodeLeft() );
            }
            else if( atom.getCoords()[axis]>currentNode.getMedian()+p_radius+MAX_ATOM_RADII ){
                if( currentNode.getNodeRight()!=null )
                    fifo.add( currentNode.getNodeRight() );
            }
            else{
                if( currentNode.getNodeLeft()!=null )
                    fifo.add( currentNode.getNodeLeft() );
                if( currentNode.getNodeRight()!=null )
                    fifo.add( currentNode.getNodeRight() );
            }

        }//

        fifo=null;
        return false;

    }

   

    /**
     * returns the neighbours in a specific distance
     *
     * @param atom
     * @param distance
     * @param list
     */
    public void neighbours(final Atom atom, final double distance, List<Atom> list){

////        if( getElements()!=null ){
////            //if is a leaf node do the test
////
////            double distance2=distance*distance;
////
////            for(Atom element : getElements() ){
////
////                if ( CalcGeom.squareDistance(atom,element) <=
////                        distance2 ) {
////                    list.add(element);
////                }
////            }
////
////            return;
////        }
////
////        double min=atom.getCoords()[getAxis()]-distance;
////        double max=atom.getCoords()[getAxis()]+distance;
////
////        if( max<=getMedian() ){
////            getNodeLeft().neighbours(atom, distance, list);
////        }
////        else if(min>=getMedian() ){
////            getNodeRight().neighbours(atom, distance, list);
////        }
////        else{
////            getNodeLeft().neighbours(atom, distance, list);
////            getNodeRight().neighbours(atom, distance, list);
////        }

        BSPTree currentNode=null;
        ArrayDeque<BSPTree> fifo=new ArrayDeque<BSPTree>();
        fifo.add(this);

        while( !fifo.isEmpty() ){

            currentNode=fifo.poll();

            if( currentNode.getElements()!=null ){
            //if is a leaf node do the test

                double distance2=distance*distance;

                for(Atom element : currentNode.getElements() ){

                    if ( CalcGeom.squareDistance(atom,element) <=
                            distance2 ) {
                        list.add(element);
                    }
                }

                continue;
            }

            double min=atom.getCoords()[currentNode.getAxis()]-distance;
            double max=atom.getCoords()[currentNode.getAxis()]+distance;

            if( max<=currentNode.getMedian() ){
                if( currentNode.getNodeLeft()!=null )
                    fifo.add( currentNode.getNodeLeft() );
            }
            else if( min>=currentNode.getMedian() ){
                if( currentNode.getNodeRight()!=null )
                    fifo.add( currentNode.getNodeRight() );
            }
            else{
                if( currentNode.getNodeLeft()!=null )
                    fifo.add( currentNode.getNodeLeft() );
                if( currentNode.getNodeRight()!=null )
                    fifo.add( currentNode.getNodeRight() );
            }

        }//

        fifo=null;

    }

    /**
     * returns if there are neighbours in a specific distance
     *
     * @param atom
     * @param distance
     * @return true if exist any atom nearer than distance
     */
    public boolean areNeighbours(final Atom atom, final double distance) {

        BSPTree currentNode=null;
        ArrayDeque<BSPTree> fifo=new ArrayDeque<BSPTree>();
        fifo.add(this);
        double distance2=distance*distance;

        while( !fifo.isEmpty() ){

            currentNode=fifo.poll();
            
            if( currentNode.getElements()!=null ){
            //if is a leaf node do the test

                for(Atom element : currentNode.getElements() ){

                    if ( CalcGeom.squareDistance(atom,element) <=
                            distance2 ) {
                        return true;
                    }
                }

                continue;
            }

            double min=atom.getCoords()[currentNode.getAxis()]-distance;
            double max=atom.getCoords()[currentNode.getAxis()]+distance;

            if( max<=currentNode.getMedian() ){
                if( currentNode.getNodeLeft()!=null )
                    fifo.add( currentNode.getNodeLeft() );
            }
            else if( min>=currentNode.getMedian() ){
                if( currentNode.getNodeRight()!=null )
                    fifo.add( currentNode.getNodeRight() );
            }
            else{
                if( currentNode.getNodeLeft()!=null )
                    fifo.add( currentNode.getNodeLeft() );
                if( currentNode.getNodeRight()!=null )
                    fifo.add( currentNode.getNodeRight() );
            }

        }//

        fifo=null;
        return false;
    }

    /**
     * returns the nearest neighbour of a point
     *
     * @param point
     * @return
     */

    public Atom nearestNeighbour(final Atom point){


        NNInfo result=nearestNeighbourInt(point,null);
        return result.atom;
    }

    private NNInfo nearestNeighbourInt(final Atom point, Double bestSqdist){

        if( this.getDepth()==0 ){
            bestSqdist=null;
        }

        if( getElements()!=null ){
            //if is a leaf node do the test

            double min_sqdist=java.lang.Double.MAX_VALUE;

            if(bestSqdist!=null)
                min_sqdist=bestSqdist;

            Atom nearest=null;

            for(Atom element : getElements() ){

                double sq_dist=CalcGeom.squareDistance(point,element);
                if ( sq_dist < min_sqdist ) {
                    min_sqdist=sq_dist;
                    nearest=element;
                }
            }

            return new NNInfo( nearest,min_sqdist);
        }

        NNInfo nearest=null;
        NNInfo nearest2=new NNInfo();

        double sqDistHiperplane=(point.getCoords()[getAxis()]-getMedian()) *
                (point.getCoords()[getAxis()]-getMedian());

        if( point.getCoords()[getAxis()] <= getMedian() ){
            nearest=this.getNodeLeft().nearestNeighbourInt(point,bestSqdist);

            if( nearest.sqDist > sqDistHiperplane ){
                nearest2=this.getNodeRight().nearestNeighbourInt(point,nearest.sqDist);
            }
        }
        else{
            nearest=this.getNodeRight().nearestNeighbourInt(point,bestSqdist);

            if( nearest.sqDist > sqDistHiperplane ){
                nearest2=this.getNodeLeft().nearestNeighbourInt(point,nearest.sqDist);
            }
        }

        return nearest2.atom==null ? nearest : nearest2;
    }

    
    

}
