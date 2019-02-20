/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.molsurface;

/**
 *
 * @author victor
 */
public class RSEdge {
    private RSVertex vertex1=null;
    private RSVertex vertex2 = null;
    private RSFace face1 = null;
    private RSFace face2 = null;

    public RSEdge() {
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
     * @return the face1
     */
    public RSFace getFace1() {
        return face1;
    }

    /**
     * @param face1 the face1 to set
     */
    public void setFace1(RSFace face1) {
        this.face1 = face1;
    }

    /**
     * @return the face2
     */
    public RSFace getFace2() {
        return face2;
    }

    /**
     * @param face2 the face2 to set
     */
    public void setFace2(RSFace face2) {
        this.face2 = face2;
    }

    
}
