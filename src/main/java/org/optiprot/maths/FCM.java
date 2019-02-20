/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

/**
 * Fuzzy C-Means clustering algorythm
 * 
 * @author victor
 */
public class FCM {

    private int dimensions = 0; //<  Number of dimensions in data
    private int numClusters = 0; //<  The number of clusters to create
    private int setSize = 0; //<  Size of the data set to cluster
    private double fuzzification = FCM_FUZZIFICATION; //<  The fuzzification parameter
    private double [][]U = null;  //<  Fuzzy partition matrix
    private double [][]M = null; //<  The previous/current prototype matrix (clusters)


    static private final double FCM_FUZZIFICATION = 2;
    static private final double FCM_MINDIST = 1e-18;

    /**
     * @return the numClusters
     */
    public int getNumClusters() {
        return numClusters;
    }

    /**
     * @return the setSize
     */
    public int getSetSize() {
        return setSize;
    }

    /**
     * @return the dimensions
     */
    public int getDimensions() {
        return dimensions;
    }

    /**
     * 
     * @return True if data has been created
     */
    protected boolean wasCreated()
    {
        if( M!=null )
            return true;
        return false;
    }

    /**
     * @return the fuzzification
     */
    public double getFuzzification() {
        return fuzzification;
    }

    /**
     * @param fuzzification the fuzzification to set
     */
    public void setFuzzification(double fuzzification) {
        this.fuzzification = fuzzification;
    }
    
    /**
     * computes the square distance
     * 
     * @param vector
     * @param idxC cluster index
     * @return
     */
    protected double sqDistance( double [] vector, int idxC ){

        double sum = 0;

        for(int d = 0; d < getDimensions(); d++){
            double temp = vector[d] - M[idxC][d];
            sum += temp*temp;
        }

        return sum;
    }

    /**
     * Returns what cluster the data value d belongs to
     *
     * @param data The data to see which cluster it belongs to
     * @return null on fail or no clusters created, otherwise the cluster index/class
     */
    public Integer getCluster(double []data)
    {
        Integer type = null;

        if(!wasCreated() || data.length != getDimensions() )
            return null;

        double min = java.lang.Double.MAX_VALUE;
        double expn=2.0/(getFuzzification()-1.0);

        for(int j = 0; j < getNumClusters(); j++){
           
            double dclust = sqDistance( data, j);

            if( dclust<FCM_MINDIST )
                return j;
            
            double sum2 = 0;

            for(int k = 0; k < getNumClusters(); k++){

                double dclust2 = sqDistance( data, k);

                sum2 += Math.pow( dclust/dclust2, expn);
            }
            
            if(sum2 < min ){
                min = sum2;
                type = j;
            }
        }

        return type;
    }

    /**
     * Allocates memory for clusters and initializes with random data.
     *
     * @param c How many clusters
     * @param d Number of dimensions in data
     * @return
     */
    public boolean createClusters(int c, int d)
    {
        if(c < 0 || d < 0)
            return false;

        if( d!=getDimensions() || c!=getNumClusters() ){
            dimensions = d;
            numClusters = c;
            M = new double [numClusters][dimensions];
        }

        for(int i = 0; i < getNumClusters(); i++){
            for(int j = 0; j < getDimensions(); j++){
                M[i][j] = Math.random();
            }
        }

        return true;
    }

    /**
     * Gets the cluster at index i
     *
     * @param i
     * @param cluster I/O vector where to save the cluster
     * @return
     */
    public boolean getCluster(int i, double []cluster)
    {
        if(i < 0 || i > getNumClusters())
            return false;

        System.arraycopy(M[i], 0, cluster, 0, getDimensions());

        return true;
    }

    /**
     * sets the cluster at index i
     * 
     * @param i
     * @param cluster
     * @return
     */
    public boolean setCluster(int i, double []cluster)
    {
        if(i < 0 || i > getNumClusters())
            return false;

        System.arraycopy(cluster, 0, M[i], 0, getDimensions());

        return true;
    }


    /**
     * 
     * @param c number of clusters
     * @param size size of data set
     */
    protected void createPartitionMatrix(int c, int size)
    {
        
        setSize = size;
        U = new double[getSetSize()][getNumClusters()];

        //  Initialize the matrix U randomly
////
////        for(int j = 0; j < getSetSize(); j++){
////
////            double sum=0.0;
////
////            for(int k = 0; k < getNumClusters(); k++){
////                //  Initialize the data
////                U[j][k] = Math.random();
////                sum += U[j][k];
////            }
////
////            sum=1.0/sum;
////
////            //normalize U[j][k]
////            for(int k = 0; k < getNumClusters(); k++){
////                U[j][k] *= sum;
////            }
////        }
         
    }


    /**
     * Creates cluster prototypes (center means) based on the 1D data
     * set passed to it.  This only works if the clusters have been
     * created/initialized first.
     *
     * This function does not have any stopping criteria, this must be
     * done by external source.  The change in clusters is returned though.
     *
     * This function is used by the overloaded public cluster function
     * which initializes the cluster prototypes based on the size of the
     * data set being clustered.
     *
     * @param dataSet Data for creating clusters
     * @param MT Translated M array (prototype array)
     * @return The change in prototypes (used for stopping criteria)
     */
    protected double cluster( double [][]dataSet, double [][]MT)
    {
        
        double e = 0.0;
        double [] data=null;
        double [] M2 = new double [getDimensions()];

        if(!wasCreated() || dataSet.length != getSetSize())
            return 0.0;

        double expn=2.0/(getFuzzification()-1.0);

        //  Update membership matrix U
        for(int i = 0; i < getSetSize(); i++){
            data = dataSet[i];

            for(int j = 0; j < getNumClusters(); j++){

                double dclust = sqDistance( data, j);
                double sum = 0.0;

                if( dclust<FCM_MINDIST ){
                    U[i][j]=1.0-FCM_MINDIST*10;
                    continue;
                }

                for(int l = 0; l < getNumClusters(); l++){

                    double dclust2 = sqDistance( data, l);

                    if( dclust2<FCM_MINDIST ){
                        U[i][j]=FCM_MINDIST;
                        sum=-1;
                        break;
                    }

                    sum += Math.pow( dclust/dclust2, expn);
                }

                if( sum>0.0 )
                    U[i][j] = 1.0/sum;
            }
        }

        //  Now we must update the prototype matrix MT
        for(int i = 0; i < getNumClusters(); i++){

            double bottom = 0;

            //backup old MT
            System.arraycopy( MT[i], 0, M2, 0, getDimensions() );


            for(int d = 0; d < getDimensions(); d++)
                MT[i][d]=0.0;

            for(int j = 0; j < getSetSize(); j++){

                double u=Math.pow( U[j][i], getFuzzification());
                bottom += u;
                data = dataSet[j];

                for(int d = 0; d < getDimensions(); d++){
                    MT[i][d] += u*data[d];
                }
            }

            bottom = 1.0/bottom;

            for(int d = 0; d < getDimensions(); d++){
                MT[i][d]*=bottom;

                //calculates variation
                e += Math.abs( MT[i][d]-M2[d] );
            }
        }

        return e;
    }


    /**
     * Uses the training set to adjust the clusters.  Clustering
     * process ends when minimum error is achieved or we have performed
     * a set number of epochs.
     *
     * @param dataSet The set to create cluster prototypes from.
     * @param minError The minimum error to try to achieve
     * @param epochs The number of epochs to try cluster on data set.
     * @return the obtained error
     */
    public double cluster( double [][]dataSet, double minError, int epochs)
    {
        double e = 0.0;

        if(!wasCreated() || dataSet.length == 0)
            return 0.0;

        createPartitionMatrix( getNumClusters(), dataSet.length );

        //  Make a copy of the cluster templates
        double [][]MT = new double[getNumClusters()][getDimensions()];

        for(int c = 0; c < getNumClusters(); c++){
            System.arraycopy(M[c], 0, MT[c], 0, getDimensions() );
        }

        for(int i = 0; i < epochs; i++){
            e = cluster(dataSet, MT);
            if(e < minError)
                break;
        }

        //  Now copy the clusters
        for(int c = 0; c < getNumClusters(); c++){
            System.arraycopy(MT[c], 0, M[c], 0, getDimensions() );
        }

        return e;
    }

    
}
