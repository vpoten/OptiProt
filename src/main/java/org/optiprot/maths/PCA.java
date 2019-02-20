/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import Jama.EigenvalueDecomposition;

// static import of all array methods : linear algebra and statistics
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.math.array.StatisticSample;
import static org.math.array.LinearAlgebra.*;
import static org.math.array.StatisticSample.*;

/**
 * Principal Component Analysis
 * 
 * @author victor
 */
public class PCA {

    // initial datas : lines = events and columns = variables
    private double[][] X=null;

	private double[] meanX, stdevX;

	private double[][] Z=null; // X centered reduced

	private double[][] cov=null; // Z covariance matrix

	private double[][] U=null; // projection matrix

	private double[] info=null; // information matrix

    

	public PCA( double[][] _X ) {
		X = _X;

		stdevX = stddeviation(X);
		meanX = mean(X);

		Z = center_reduce(X);

		cov = StatisticSample.covariance(Z);

		EigenvalueDecomposition e = eigen(cov);
		U = e.getV().transpose().getArray();

        // covariance matrix is symetric, so only real eigenvalues...
		info = e.getRealEigenvalues(); 
	}


    /**
     * calculates PCA and transforms the coordinates in chain using
     * the eigenvectors of covarianze matrix.
     * Precondition : the centroid of the chain is (0,0,0)
     *
     * @param chain I/O
     * @return : the matrix of PCA eigenvectors (rows)
     * @throws org.biojava.bio.structure.StructureException
     */
    public static double [][] transformPCA( Chain chain ) throws StructureException {

        List<Atom> list=BSPTree.getListAtoms(chain);
        Atom [] atomSet=list.toArray(new Atom [list.size()]);

        double [][] _X=new double [atomSet.length][3];

        for( int i=0;i<atomSet.length;i++ ){
            _X[i][0]=atomSet[i].getX();
            _X[i][1]=atomSet[i].getY();
            _X[i][2]=atomSet[i].getZ();
        }

        double [][] _cov = PCA.covariance(_X);

        EigenvalueDecomposition e = eigen(_cov);
		double [][]_U = e.getV().transpose().getArray();

        // covariance matrix is symetric, so only real eigenvalues...
		double []_info = e.getRealEigenvalues();

        double [] temp=null;
        double itemp=0;

        //order the eigenvectors by eigenvalue
        if( _info[1]>_info[0] ){
            temp=_U[0];
            itemp=_info[0];
            _U[0]=_U[1];
            _info[0]=_info[1];
            _U[1]=temp;
            _info[1]=itemp;
        }

        if( _info[2]>_info[1] ){
            temp=_U[1];
            itemp=_info[1];
            _U[1]=_U[2];
            _info[1]=_info[2];
            _U[2]=temp;
            _info[2]=itemp;
        }

        if( _info[1]>_info[0] ){
            temp=_U[0];
            _U[0]=_U[1];
            _U[1]=temp;
        }

        CalcTransform.applyTransform(atomSet, _U);

        return _U;

    }

    /**
     * 
     * @param atoms : I/O
     * @param pca
     * @return
     */
    public static Atom[] inversePCA( Atom [] atoms, double [][] pca ){

        Atom c1=new AtomImpl();
        c1.setCoords(pca[0]);

        Atom c2=new AtomImpl();
        c2.setCoords(pca[1]);

        Atom c3=new AtomImpl();
        c3.setCoords(pca[2]);

        for( Atom at : atoms ){

            Atom vector=CalcGeom.product(c1, at.getX());

            CalcGeom.addEquals(vector, CalcGeom.product(c2, at.getY()));
            CalcGeom.addEquals(vector, CalcGeom.product(c3, at.getZ()));

            at.setX(vector.getX());
            at.setY(vector.getY());
            at.setZ(vector.getZ());
        }

        return atoms;
    }


	// normalization of x relatively to X mean and standard deviation
	protected double[][] center_reduce(double[][] x) {
		double[][] y = new double[x.length][x[0].length];
		for (int i = 0; i < y.length; i++)
			for (int j = 0; j < y[i].length; j++)
				y[i][j] = (x[i][j] - meanX[j]) / stdevX[j];
		return y;
	}

	// de-normalization of y relatively to X mean and standard deviation
	protected double[] inv_center_reduce(double[] y) {
		return inv_center_reduce(new double[][] { y })[0];
	}

	// de-normalization of y relatively to X mean and standard deviation
	protected double[][] inv_center_reduce(double[][] y) {
		double[][] x = new double[y.length][y[0].length];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				x[i][j] = (y[i][j] * stdevX[j]) + meanX[j];
		return x;
	}


    /**
     * calculates covarianze matrix, assumes data with mean 0
     *
     * @param v
     * @return
     */
    protected static double[][] covariance(double[][] v) {
        int m = v.length;
        int n = v[0].length;
        double[][] C = new double[n][n];
        double _inv_degrees = 1.0/(m - 1);
        double c;

        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                c = 0;

                for (int k = 0; k < m; k++)
                    c += v[k][i] * v[k][j];
                
                C[i][j] = c * _inv_degrees;
                C[j][i] = C[i][j];
            }
        }
        return C;
    }

}
