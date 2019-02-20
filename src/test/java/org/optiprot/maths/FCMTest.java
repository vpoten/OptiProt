/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.maths;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author victor
 */
public class FCMTest {

    public FCMTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    
    /**
     * Test of wasCreated method, of class FCM.
     */
    @Test
    public void testWasCreated() {
        System.out.println("wasCreated");
        FCM instance = new FCM();
        boolean expResult = false;
        boolean result = instance.wasCreated();
        assertEquals(expResult, result);
    }


    /**
     * Test of sqDistance method, of class FCM.
     */
    @Test
    public void testSqDistance() {
        System.out.println("sqDistance");

        double[] vector = new double [2];
        vector[0] = vector[1] =0.0;

        int idxC = 0;
        FCM instance = new FCM();
        instance.createClusters(1, 2);
        
        double result = instance.sqDistance(vector, idxC);

        assertTrue( result<2.0 );
    }

    /**
     * Test of getCluster method, of class FCM.
     */
    @Test
    public void testGetCluster_doubleArr() {
        System.out.println("getCluster");

        double[] data = new double [2];
        data[0] = data[1] = 0.0;

        FCM instance = new FCM();
        instance.createClusters(1, 2);
        Integer expResult = 0;
        Integer result = instance.getCluster(data);

        assertEquals(expResult, result);

        instance.createClusters(2, 2);
        instance.setCluster(expResult, new double [] {1e-3,1e-2} );

        result = instance.getCluster(data);
        assertEquals(expResult, result);
    }

    /**
     * Test of createClusters method, of class FCM.
     */
    @Test
    public void testCreateClusters() {
        System.out.println("createClusters");
        int c = 2;
        int d = 2;
        FCM instance = new FCM();
        boolean expResult = true;
        boolean result = instance.createClusters(c, d);

        assertEquals(expResult, result);
    }

    /**
     * Test of getCluster method, of class FCM.
     */
    @Test
    public void testGetCluster_int_doubleArr() {
        System.out.println("getCluster");
        int i = 0;
        int d=2;

        double[] cluster = new double [d];
        
        FCM instance = new FCM();
        instance.createClusters(2, d);

        boolean expResult = false;
        boolean result = instance.getCluster(-1, cluster);
        assertEquals(expResult, result);

        instance.getCluster(i, cluster);

        assertTrue( instance.sqDistance(cluster, i)<1e-12 );
    }

    /**
     * Test of setCluster method, of class FCM.
     */
    @Test
    public void testSetCluster() {
        System.out.println("setCluster");
        int i = 0;
        int d=2;

        double[] cluster = new double [d];
        cluster[0] = cluster[1] = 2.5;

        FCM instance = new FCM();
        instance.createClusters(3, d);

        boolean expResult = true;
        boolean result = instance.setCluster(i, cluster);
        assertEquals(expResult, result);

        assertTrue( instance.sqDistance(cluster, i)<1e-12 );
    }

    /**
     * Test of createPartitionMatrix method, of class FCM.
     */
    @Test
    public void testCreatePartitionMatrix() {
        System.out.println("createPartitionMatrix");
        int c = 2;
        int size = 5;
        FCM instance = new FCM();
        instance.createPartitionMatrix(c, size);

        assertTrue( instance.getSetSize()==size );
    }

    

    /**
     * Test of cluster method, of class FCM.
     */
    @Test
    public void testCluster_3args() {
        System.out.println("cluster");
        
        int c=3;
        int d=3;
        int n=30;

        double[][] dataSet = new double [n][d];
        
        double [] center1=new double [] {-1,-1,1};
        double [] center2=new double [] {10,10,0};
        double [] center3=new double [] {4,4,-2};

        //generate the dataset
        for(int i=0;i<n;i++){

            double [] center=null;

            if( i%c==1 )
                center=center1;
            else if( i%c==2 )
                center=center2;
            else
                center=center3;
            
            for(int j=0;j<d;j++){
                dataSet[i][j]=Math.random()-0.5+center[j];
            }
        }
        
        double minError = 1e-3;
        int epochs = 100;

        FCM instance = new FCM();
        instance.createClusters(c, d);

        //for(int i=0;i<c;i++){
        //    instance.setCluster(i, dataSet[ (int)Math.floor(Math.random()*n)] );
        //}

        ////instance.setFuzzification(3);
        double result = instance.cluster(dataSet, minError, epochs);
        
        assertTrue( result<minError );


    }

    @Test
    public void testCluster_3args_2() {
        System.out.println("cluster2");

        int c=5;
        int d=5;
        int n=25;

        double[][] dataSet = new double [n][d];

        dataSet[0]=new double [] { 90.1, 2.6, 1.0, 6.9, 0.35};
        dataSet[1]=new double [] { 88.5, 1.4, 3.5, 6.0, 0.24};
        dataSet[2]=new double [] { 88.4, 2.2, 2.7, 6.4, 0.18};
        dataSet[3]=new double [] { 90.3, 1.7, 1.4, 6.2, 0.40};
        dataSet[4]=new double [] { 90.4, 0.6, 4.5, 4.4, 0.10};
        dataSet[5]=new double [] { 87.7, 3.5, 3.4, 4.8, 0.71};
        dataSet[6]=new double [] { 86.9, 4.8, 1.7, 5.7, 0.90};
        dataSet[7]=new double [] { 82.1, 5.9, 7.9, 4.7, 0.78};
        dataSet[8]=new double [] { 81.9, 7.4, 7.2, 2.7, 0.85};
        dataSet[9]=new double [] { 81.6, 10.1, 6.3, 4.4, 0.75};
        dataSet[10]=new double [] { 81.6, 6.6, 5.9, 4.9, 0.93};
        dataSet[11]=new double [] { 86.5, 3.9, 3.2, 5.6, 0.80};
        dataSet[12]=new double [] { 90.0, 2.0, 1.8, 5.5, 0.47};
        dataSet[13]=new double [] { 82.8, 7.1, 5.1, 3.7, 1.10};
        dataSet[14]=new double [] { 86.2, 3.0, 4.8, 5.3, 0.70};
        dataSet[15]=new double [] { 82.0, 5.6, 6.4, 4.7, 0.91};
        dataSet[16]=new double [] { 76.3, 9.3, 9.5, 3.0, 1.20};
        dataSet[17]=new double [] { 70.7, 3.6, 17.6, 5.6, 0.63};
        dataSet[18]=new double [] { 71.3, 12.3, 13.1, 1.9, 2.30};
        dataSet[19]=new double [] { 72.5, 9.2, 12.6, 3.3, 1.40};
        dataSet[20]=new double [] { 65.9, 10.4, 19.7, 2.6, 1.40};
        dataSet[21]=new double [] { 64.8, 10.7, 20.3, 2.5, 1.40};
        dataSet[22]=new double [] { 64.8, 11.1, 21.2, 1.6, 1.70};
        dataSet[23]=new double [] { 46.4, 9.7, 42.0, 0.0, 0.85};
        dataSet[24]=new double [] { 44.9, 10.6, 34.9, 0.9, 0.53};


        double minError = 1e-3;
        int epochs = 100;

        FCM instance = new FCM();
        instance.createClusters(c, d);

        for(int i=0;i<c;i++){
            instance.setCluster(i, dataSet[ (int)Math.floor(Math.random()*n)] );
        }

        double result = instance.cluster(dataSet, minError, epochs);

        assertTrue( result<minError );

    }

}