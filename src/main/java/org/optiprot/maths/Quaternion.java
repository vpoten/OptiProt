/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.biojava.bio.structure.*;

/**
 *
 * @author victor
 */
public class Quaternion {

    private final double x0, x1, x2, x3;

    // create a new object with the given components
    public Quaternion(double x0, double x1, double x2, double x3) {
        this.x0 = x0;
        this.x1 = x1;
        this.x2 = x2;
        this.x3 = x3;
    }

    /** Axis does not need to be normalized but must not be the zero
      vector. Angle is in radians. */
    public Quaternion(Atom axis, double angle) {
        double halfTheta = angle * 0.5;
        x0 = Math.cos(halfTheta);
        double sinHalfTheta = Math.sin(halfTheta);
        Atom realAxis = Calc.unitVector(axis);
        
        x1 = realAxis.getX() * sinHalfTheta;
        x2 = realAxis.getY() * sinHalfTheta;
        x3 = realAxis.getZ() * sinHalfTheta;
    }

    public double getAngle() {
        return 2 * Math.acos(x0);
    }

    public void getAxis(double []vector){
        
        double temp=Math.sqrt(1.0-x0*x0);

        if (temp < 0.001) {
            vector[0]=x1;
            vector[1]=x2;
            vector[2]=x3;
            return;
        }

        temp=1.0/temp;

        vector[0]=x1*temp;
        vector[1]=x2*temp;
        vector[2]=x3*temp;
    }

    // return a string representation of the invoking object
    @Override
    public String toString() {
        return x0 + " + " + x1 + "i + " + x2 + "j + " + x3 + "k";
    }

    // return the quaternion norm
    public double norm() {
        return Math.sqrt(x0*x0 + x1*x1 +x2*x2 + x3*x3);
    }

    // return the quaternion conjugate
    public Quaternion conjugate() {
        return new Quaternion(x0, -x1, -x2, -x3);
    }

    // return a new Quaternion whose value is (this + b)
    public Quaternion plus(Quaternion b) {
        Quaternion a = this;
        return new Quaternion(a.x0+b.x0, a.x1+b.x1, a.x2+b.x2, a.x3+b.x3);
    }


    // return a new Quaternion whose value is (this * b)
    public Quaternion times(Quaternion b) {
        Quaternion a = this;
        double y0 = a.x0*b.x0 - a.x1*b.x1 - a.x2*b.x2 - a.x3*b.x3;
        double y1 = a.x0*b.x1 + a.x1*b.x0 + a.x2*b.x3 - a.x3*b.x2;
        double y2 = a.x0*b.x2 - a.x1*b.x3 + a.x2*b.x0 + a.x3*b.x1;
        double y3 = a.x0*b.x3 + a.x1*b.x2 - a.x2*b.x1 + a.x3*b.x0;
        return new Quaternion(y0, y1, y2, y3);
    }

    // return a new Quaternion whose value is the inverse of this
    public Quaternion inverse() {
        double d = x0*x0 + x1*x1 + x2*x2 + x3*x3;
        return new Quaternion(x0/d, -x1/d, -x2/d, -x3/d);
    }


    // return a / b
    public Quaternion divides(Quaternion b) {
         Quaternion a = this;
        return a.inverse().times(b);
    }

    /** Rotate a vector by this quaternion. Implementation is from
      Horn's <u>Robot Vision</u>. NOTE: src and dest must be different
      vectors. */
    public Atom rotateVector(Atom src) {
        Atom qVec = new AtomImpl();
        qVec.setX(x1);
        qVec.setY(x2);
        qVec.setZ(x3);

        Atom qCrossX = Calc.vectorProduct(qVec,src);
        Atom qCrossXCrossQ = Calc.vectorProduct(qCrossX,qVec);

        CalcGeom.product2(qCrossX, 2.0 * x0);

        CalcGeom.product2(qCrossXCrossQ, -2.0 );

        Atom dest=Calc.add(src, qCrossX);
        CalcGeom.addEquals(dest, qCrossXCrossQ);

        return dest;
    }


    /**
     * returns a normalized quaternion
     * @return
     */
    public Quaternion normalize() {

        double norm=this.norm();
        norm=1.0/norm;

        return new Quaternion(x0*norm,x1*norm,x2*norm,x3*norm);
    }

    
}
