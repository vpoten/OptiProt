/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.molsurface;

import org.biojava.bio.structure.Atom;

/**
 *
 * @author victor
 */
public class RSFace {

    private RSVertex vertex1=null;
    private RSVertex vertex2 = null;
    private RSVertex vertex3 = null;
    private RSEdge edge1 = null;
    private RSEdge edge2 = null;
    private RSEdge edge3 = null;

    private Atom normal=null;
    private Atom center=null;

    public RSFace() {
    }

    /**
     * @return the vertex1
     */
    public RSVertex getVertex1() {
        return vertex1;
    }

    /**
     * @param vertex1 the vertex1 to set
     */
    public void setVertex1(RSVertex vertex1) {
        this.vertex1 = vertex1;
    }

    /**
     * @return the vertex2
     */
    public RSVertex getVertex2() {
        return vertex2;
    }

    /**
     * @param vertex2 the vertex2 to set
     */
    public void setVertex2(RSVertex vertex2) {
        this.vertex2 = vertex2;
    }

    /**
     * @return the vertex3
     */
    public RSVertex getVertex3() {
        return vertex3;
    }

    /**
     * @param vertex3 the vertex3 to set
     */
    public void setVertex3(RSVertex vertex3) {
        this.vertex3 = vertex3;
    }

    /**
     * @return the normal
     */
    public Atom getNormal() {
        return normal;
    }

    /**
     * @param normal the normal to set
     */
    public void setNormal(Atom normal) {
        this.normal = normal;
    }

    /**
     * @return the center
     */
    public Atom getCenter() {
        return center;
    }

    /**
     * @param center the center to set
     */
    public void setCenter(Atom center) {
        this.center = center;
    }

    /**
     * @return the edge1
     */
    public RSEdge getEdge1() {
        return edge1;
    }

    /**
     * @param edge1 the edge1 to set
     */
    public void setEdge1(RSEdge edge1) {
        this.edge1 = edge1;
    }

    /**
     * @return the edge2
     */
    public RSEdge getEdge2() {
        return edge2;
    }

    /**
     * @param edge2 the edge2 to set
     */
    public void setEdge2(RSEdge edge2) {
        this.edge2 = edge2;
    }

    /**
     * @return the edge3
     */
    public RSEdge getEdge3() {
        return edge3;
    }

    /**
     * @param edge3 the edge3 to set
     */
    public void setEdge3(RSEdge edge3) {
        this.edge3 = edge3;
    }

    /**
     * get the vertex of the face that not is in the edge
     * 
     * @param edge
     */
    public RSVertex getOtherVertex( RSEdge edge){

         //
        if( !getVertex1().equals( edge.getVertex1() ) &&
                !getVertex1().equals( edge.getVertex2() ) ){
            return getVertex1();
        }
        else if( !getVertex2().equals( edge.getVertex1() ) &&
                !getVertex2().equals( edge.getVertex2() ) ){
            return getVertex2();
        }

        return getVertex3();
    }

}
