/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.pockets;

/**
 *
 * @author victor
 */
class IntCounter {

    private int value=0;

    public IntCounter(int value) {
        this.value = value;
    }

    public void decrement(){
        value--;
    }

    public int getValue(){
        return this.value;
    }
}
