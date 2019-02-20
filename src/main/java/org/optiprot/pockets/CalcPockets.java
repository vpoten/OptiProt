/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.pockets;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.StructureException;
import org.optiprot.maths.AtomGrid;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcSphere;
import org.optiprot.potential.IForceField;
import org.optiprot.potential.TestForceField;

/**
 *
 * @author victor
 */
public class CalcPockets {

    private static double GRID_RESOL=1.0;//in angstroms
    private static double GRID_RESOL_INC=0.1;//in angstroms
    private static double OUT_DIST = 4.5;
    

    private static double BORDER_VALUE = 5.0;
    private static double SOLVENT_VALUE = -1.0;
    private static double VOID_VALUE = 1.0;

    private static double RAY_LENGTH = 10.0;
    private static double RAY_WIDTH = 2.0;

    //buried index limits 
    private static int HIGH_LIMIT = 26;
    private static int LOW_LIMIT = 15;

    private static Atom[] rays=null;

    static{
        rays=triangOctahedron();
    }



    /**
     *
     * @param chain
     * @param grid
     * @param npockets
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static GridProbeCluster[] calc( Chain chain, AtomGrid grid, int npockets )
            throws StructureException {

        GridProbeCluster[] pockets=null;


        IForceField ffield = new TestForceField();

        //generate the grid
        BSPTree btree=new BSPTree( chain );
        generateGrid(btree, ffield, grid);

        // get the grid points that are outside the VdW volume and its distance
        // to the protein are less than a value (OUT_DIST = 4.5A)
        int lim=grid.getNpointx()*grid.getNpointy()*grid.getNpointz();
        
        ArrayList<GridProbe> gridProbes=new ArrayList<GridProbe>();

        for(int i=0;i<lim;i++){

            Atom p=grid.getPoint(i);

            if( !grid.isInVolume(i) ){
                boolean areNeigh = btree.areNeighbours( p, OUT_DIST);

                if( areNeigh ){
                    gridProbes.add( new GridProbe(i, p, 0) );
                }
            }
        }

        //////test - return all probes
        ////pockets=new GridProbeCluster [1];
        ////pockets[0]=new GridProbeCluster();
        ////pockets[0].getProbes().addAll(gridProbes);
        ////return pockets;

        // Calculation of buriedness values of grid probes installed in areas
        // closely above the protein surface

        calcBuriedIndices( gridProbes, btree );

        // Clustering of adjoining grid probes indicating buried regions of the
        // structure to find potential binding-sites.
        ArrayList<GridProbeCluster> clusters =
                calcClusters( gridProbes, grid, LOW_LIMIT, HIGH_LIMIT );

        if( clusters==null )
            return null;
        

        // get the 3 largest clusters (sort in ascent order)
        Collections.sort( clusters, new GridProbeClusterComp() );

        if( clusters.size()>=npockets ){
            pockets=new GridProbeCluster[npockets];
        }
        else{
            pockets=new GridProbeCluster[clusters.size()];
        }


        for( int i=0; i<pockets.length; i++ ){
            pockets[i] = clusters.get(clusters.size()-i-1);
        }


        return pockets;
    }

    /**
     * Calculation of buriedness values of grid probes
     *
     * @param gridProbes
     * @param btree
     */
    private static void calcBuriedIndices( List<GridProbe> gridProbes, BSPTree btree ){

        double sqwidth=RAY_WIDTH*0.5;
        sqwidth*=sqwidth;

        ArrayList<Atom> neighbors=new ArrayList<Atom>();
        
        Atom p=new AtomImpl();
        Atom pq=new AtomImpl();
        Atom vr=new AtomImpl();

        //for all grid probes
        for( GridProbe probe : gridProbes ){

            //get the neigbors of probe i
            btree.neighbours( probe.getAtom(), RAY_LENGTH , neighbors);

            if( neighbors.isEmpty() )
                continue;

            p.setX( probe.getAtom().getX()*-1.0 );
            p.setY( probe.getAtom().getY()*-1.0 );
            p.setZ( probe.getAtom().getZ()*-1.0 );

            //scan in the directions stored in rays
            for (int i=0;i<rays.length;i++ ){
                for( Atom at : neighbors ){

                    pq.setX( at.getX() );
                    pq.setY( at.getY() );
                    pq.setZ( at.getZ() );

                    //calculates q-p
                    CalcGeom.addEquals(pq, p);

                    CalcGeom.vectorProduct( vr, pq, rays[i]);
                    double sqd = CalcGeom.squareLength(vr);
                    double t = Calc.skalarProduct(pq, rays[i]);

                    if( sqd<sqwidth && t>0.0 && t<RAY_LENGTH ){
                        //if a atom is encountered in the search ray
                        probe.incBuriedIdx(1); //increase buried index
                        break; // scan the next ray
                    }
                }

                if( (rays.length-i+probe.getBuriedIdx())<LOW_LIMIT ){
                    break;///if cannot reach the low limit
                }
            }

            neighbors.clear();
        }
    }

    /**
     * Clustering of adjoining grid probes
     *
     * @param gridProbes
     * @param grid
     * @param lowlim
     * @param highlim
     * @return
     */
    private static ArrayList<GridProbeCluster> calcClusters(
            ArrayList<GridProbe> gridProbes, AtomGrid grid, int lowlim, int highlim ){

        ArrayList<GridProbeCluster> clusters=new ArrayList<GridProbeCluster>();

        //for all grid probes
        for( GridProbe probe : gridProbes ){

            if( probe.getBuriedIdx()>highlim ||
                    probe.getBuriedIdx()<lowlim ){
                // this probe doesnt belong to any pocket
                continue;
            }

            boolean adjoin=false;

            for( GridProbeCluster cluster : clusters ){
                if( cluster.isAdjoin(probe, grid) ){
                    //add the probe to the cluster
                    cluster.getProbes().add(probe);
                    adjoin=true;
                    break;
                }
            }

            if( adjoin )
                continue;

            //create a new cluster with this probe
            GridProbeCluster cluster=new GridProbeCluster();
            cluster.getProbes().add(probe);
            clusters.add( cluster );

        }

        if( clusters.isEmpty() )
            return null;

        //join conected clusters

        boolean join=true;
        int lim=0;

        while( join ){

            GridProbeCluster cluster=clusters.get( clusters.size()-1 );
            lim=clusters.size()-1;
            join=false;

            for(int i=0;i<lim;i++){
                if( cluster.isAdjoin(clusters.get(i), grid) ){
                    //if the two clusters are conected join them
                    cluster.join( clusters.get(i) );
                    clusters.remove(i);
                    join=true;
                    break;
                }
            }
        }

        return clusters;
    }


    /**
     *
     * @param btree
     * @param ffield
     * @param grid
     * @return
     */
    private static AtomGrid generateGrid(BSPTree btree, IForceField ffield, AtomGrid grid) {

        Atom [] bbox=btree.getBBox();

        double gridResol=GRID_RESOL;

        int npointx=(int) Math.ceil((bbox[1].getX() - bbox[0].getX()) / gridResol);
        int npointy=(int) Math.ceil((bbox[1].getY() - bbox[0].getY()) / gridResol);
        int npointz=(int) Math.ceil((bbox[1].getZ() - bbox[0].getZ()) / gridResol);

        //adjust grid resolution if the dimension is bigger
        while( npointx*npointy*npointz > grid.size() ){

            gridResol+=GRID_RESOL_INC;
            npointx=(int) Math.ceil((bbox[1].getX() - bbox[0].getX()) / gridResol);
            npointy=(int) Math.ceil((bbox[1].getY() - bbox[0].getY()) / gridResol);
            npointz=(int) Math.ceil((bbox[1].getZ() - bbox[0].getZ()) / gridResol);
        }


        grid.clear();
        grid.setDimension(npointx,npointy,npointz);
        grid.setResol(gridResol);

        double minX=bbox[0].getX();
        double minY=bbox[0].getY();
        double minZ=bbox[0].getZ();
        Atom point=new AtomImpl();

        // construct the grid points
        for(int i=0;i<npointx;i++){
            for(int j=0;j<npointy;j++){
                for(int k=0;k<npointz;k++){

                    int idx=grid.getIndex(i, j, k);

                    point.setX( minX+i*gridResol);
                    point.setY( minY+j*gridResol);
                    point.setZ( minZ+k*gridResol);

                    grid.setPoint(idx, point, false);
                }
            }
        }

        //fill the grid starting in a corner, borders stop the propagation
        fillExterior( grid, btree, ffield, null );

        int lim=npointx*npointy*npointz;
        int value=0;

        //assign points inside atoms volume
        for(int i=0;i<lim;i++){

            double val=grid.getValue(i);

            if( val!=BORDER_VALUE && val!=SOLVENT_VALUE ){
                if( btree.isInsideVolume( grid.getPoint(i), ffield, 0.0) ){
                    //if exists an atom center nearer than atom radii + value
                    grid.setValue( i, BORDER_VALUE);
                    grid.setInVolume( i, true);
                }
                else{
                    value++;
                }
            }
        }

        //scan for interior pockets
        IntCounter counter=new IntCounter( value );

        while( counter.getValue()>0 ){

            for(int i=0;i<lim;i++){

                double val=grid.getValue(i);

                if( val!=BORDER_VALUE && val!=SOLVENT_VALUE
                        && val!=VOID_VALUE ){
                    counter.decrement();
                    fillVoids( i, grid, btree, ffield, counter );
                }
            }

        }//
        

        

        return grid;
    }

    /**
     * fill the interior voids, classifies into voids and borders
     *
     * @param seed
     * @param grid
     * @param btree
     * @param ffield
     * @param counter
     */
    private static void fillVoids( int seed, AtomGrid grid, BSPTree btree, 
            IForceField ffield, IntCounter counter ){

        ArrayList<Atom> listNeigh=new ArrayList<Atom>();
        ArrayDeque<Integer> fifo=new ArrayDeque<Integer>();
        
        // first classify seed
        int ninter = numInterSphereNeighbors( grid.getPoint(seed),
                        btree, ffield, listNeigh );

        if( ninter==0 ){
            // if not collision
            grid.setValue( seed, VOID_VALUE);
            fifo.add(seed);
        }
        else{
            grid.setValue( seed, BORDER_VALUE);
            grid.setInVolume( seed, true);
            listNeigh.clear();
            return;
        }


        int [] neighbors=new int [GridProbe.NUM_NEIGHBORS];
        

        while( !fifo.isEmpty() ) {

            int curridx=fifo.poll();

            GridProbe.getNeighbors(curridx, grid, neighbors);

            // for all the grid points connected with current
            for( int i=0; i<neighbors.length; i++){

                if( neighbors[i]<0 )
                    continue;

                double val=grid.getValue(neighbors[i]);

                if( val==BORDER_VALUE || val==SOLVENT_VALUE || val==VOID_VALUE )
                    continue;// if visited

                counter.decrement();

                ninter = numInterSphereNeighbors( grid.getPoint(neighbors[i]),
                        btree, ffield, listNeigh );

                if( ninter==0 ){
                    // if not collision
                    grid.setValue( neighbors[i], VOID_VALUE);
                    fifo.add( neighbors[i] );
                }
                else{
                    grid.setValue( neighbors[i], BORDER_VALUE);
                    grid.setInVolume( neighbors[i], true);
                    listNeigh.clear();
                }


            }//end loop neighbors

        }

    }

    /**
     * checks if a sphere with radii 1.4 is able to reach the
     * neighbor position from current one (without collision)
     *
     * @param point
     * @param btree
     * @param ffield
     * @param listNeigh
     * @return
     */
    private static int numInterSphereNeighbors(Atom point, BSPTree btree,
            IForceField ffield, ArrayList<Atom> listNeigh ){

        int ninter=0;

        btree.neighbours( point, 3.4, listNeigh);//1.4+2.0 A

        if( !listNeigh.isEmpty() ){

            for( Atom at : listNeigh ){
                if( CalcSphere.isIntersection(point, 1.4, at,
                        ffield.getVDWRadius(at)) ){
                        ninter++;
                        break;
                }
            }
        }

        return ninter;
    }

    /**
     * classifies the grid point into solvent or border (exterior)
     *
     * @param grid
     * @param btree
     * @param ffield
     * @param start = if not null, starts the procedure at this point
     */
    private static void fillExterior( AtomGrid grid, BSPTree btree, IForceField ffield, Atom start ){

        int seed=-1;

        if( start==null ){
            int [] seeds =
            {
                grid.getIndex(0, 0, 0), grid.getIndex(0, 0, grid.getNpointz()-1),
                grid.getIndex(0, grid.getNpointy()-1, 0), grid.getIndex(grid.getNpointx()-1,0,0),
                grid.getIndex(grid.getNpointx()-1, 0, grid.getNpointz()-1),
                grid.getIndex( 0, grid.getNpointy()-1, grid.getNpointz()-1),
                grid.getIndex(grid.getNpointx()-1, grid.getNpointy()-1, 0),
                grid.getIndex(grid.getNpointx()-1, grid.getNpointy()-1, grid.getNpointz()-1)
            };



            for(int i=0;i<seeds.length;i++){
                if( !btree.isInsideVolume( grid.getPoint(seeds[i]), ffield, 0.7) ){
                    seed=seeds[i];
                    break;
                }
            }

            if( seed==-1 )
                return;
        }
        else{
            Atom origin=grid.getPoint(0);
            double resol=grid.getResol();
            seed=grid.getIndex( 
                    (int) Math.round((start.getX()-origin.getX())/resol),
                    (int) Math.round((start.getY()-origin.getY())/resol),
                    (int) Math.round((start.getZ()-origin.getZ())/resol)
                    );
        }

        ArrayDeque<Integer> fifo=new ArrayDeque<Integer>();
        fifo.add(seed);
        grid.setValue( seed, SOLVENT_VALUE);
        

        int [] neighbors=new int [GridProbe.NUM_NEIGHBORS];
        ArrayList<Atom> listNeigh=new ArrayList<Atom>();

        while( !fifo.isEmpty() ) {

            int curridx=fifo.poll();

            GridProbe.getNeighbors(curridx, grid, neighbors);

            // for all the grid points connected with current
            for( int i=0; i<neighbors.length; i++){

                if( neighbors[i]<0 )
                    continue;

                double val=grid.getValue(neighbors[i]);

                if( val==BORDER_VALUE || val==SOLVENT_VALUE )
                    continue;// if visited

                if( btree.isInsideVolume( grid.getPoint(neighbors[i]), ffield, 0.0) ){
                    //if exists an atom center nearer than atom radii + value
                    grid.setValue( neighbors[i], BORDER_VALUE);
                    grid.setInVolume( neighbors[i], true);
                }
                else{
                    //checks if a sphere with radii 1.4 is able to reach the
                    // neighbor position from current one (without collision)

                    Atom mid = Calc.add( grid.getPoint(neighbors[i]), grid.getPoint(curridx));
                    CalcGeom.product2(mid, 0.5);

                    int ninter = numInterSphereNeighbors( mid, btree, ffield, listNeigh );

                    if( ninter==0 ){
                        // if not collision
                        grid.setValue( neighbors[i], SOLVENT_VALUE);
                        fifo.add( neighbors[i] );
                    }
                    else{
                        listNeigh.clear();
                    }

                }

                
                
            }//end loop neighbors
            
        }

    }
    
//    /**
//     * scan the grid for interior region of the protein and
//     * increases the density of points between two borders
//     *
//     * @param btree
//     * @param ffield
//     * @param grid
//     * @param axis1 : axis order, if 0,1,2 scan z (fixed x,y)
//     * @param axis2
//     * @param axis3
//     */
//    private static void scanGridDensity( BSPTree btree, AtomGrid grid,
//            int axis1, int axis2, int axis3 ){
//
//        int npoint1, npoint2, npoint3;
//        int perm=0;
//
//
//        if( axis1==0 && axis2==1 ){ //x,y,z
//            npoint1=grid.getNpointx();
//            npoint2=grid.getNpointy();
//            npoint3=grid.getNpointz();
//            perm=0;
//        }
//        else if( axis1==1 && axis2==2){ //y,z,x
//            npoint1=grid.getNpointy();
//            npoint2=grid.getNpointz();
//            npoint3=grid.getNpointx();
//            perm=1;
//        }
//        else{//z,x,y
//            npoint1=grid.getNpointz();
//            npoint2=grid.getNpointx();
//            npoint3=grid.getNpointy();
//            perm=2;
//        }
//
//        ArrayList<Integer> borders=new ArrayList<Integer>();
//        int idx=0;
//
//        for(int i=0;i<npoint1;i++){
//            for(int j=0;j<npoint2;j++){
//
//                borders.clear();
//
//                for(int k=0;k<npoint3;k++){
//
//                    switch(perm){
//                        case 0: idx=grid.getIndex(i, j, k); break;
//                        case 1: idx=grid.getIndex(k, i, j); break;
//                        case 2: idx=grid.getIndex(j, k, i); break;
//                    }
//
//                    if( grid.getValue(idx)==BORDER_VALUE ){
//                        borders.add(k);
//                    }
//                }
//
//
//                //increases the density of points between two borders
//                while( borders.size()>1 ){
//                    int start=borders.get(0)+1;
//                    int end=borders.get(1);
//
//                    for(int k=start;k<end;k++){
//
//                        switch(perm){
//                            case 0: idx=grid.getIndex(i, j, k); break;
//                            case 1: idx=grid.getIndex(k, i, j); break;
//                            case 2: idx=grid.getIndex(j, k, i); break;
//                        }
//                        grid.setValue(idx, grid.getValue(idx)+1 );
//                    }
//
//                    borders.remove(0);
//                }
//            }
//        }
//
//    }

    /**
     * triangulate an octahedron
     *
     * @return : an array of 30 vectors
     */
    private static Atom [] triangOctahedron(){

        Atom [] l_rays=new Atom [30];

        //create the 6 vertices of the octahedron
        l_rays[0]=new AtomImpl();
        l_rays[0].setX(0); l_rays[0].setY(1); l_rays[0].setZ(0);

        l_rays[1]=new AtomImpl();
        l_rays[1].setX(1); l_rays[1].setY(0); l_rays[1].setZ(0);

        l_rays[2]=new AtomImpl();
        l_rays[2].setX(0); l_rays[2].setY(0); l_rays[2].setZ(-1);

        l_rays[3]=new AtomImpl();
        l_rays[3].setX(-1); l_rays[3].setY(0); l_rays[3].setZ(0);

        l_rays[4]=new AtomImpl();
        l_rays[4].setX(0); l_rays[4].setY(0); l_rays[4].setZ(1);

        l_rays[5]=new AtomImpl();
        l_rays[5].setX(0); l_rays[5].setY(-1); l_rays[5].setZ(0);

        Atom a=null;
        Atom b=null;
        Atom c=null;

        Atom ab=null;
        Atom bc=null;
        Atom ca=null;

        for(int i=0;i<4;i++){

            a=l_rays[0];
            b=l_rays[i+1];
            c=l_rays[((i+1)%4)+1];

            ab=Calc.add(a, b);
            CalcGeom.product2(ab, 0.5);

            bc=Calc.add(b, c);
            CalcGeom.product2(bc, 0.5);

            ca=Calc.add(c, a);
            CalcGeom.product2(ca, 0.5);

            l_rays[5+i*6+1]=Calc.unitVector( Calc.add(ab, bc) );
            l_rays[5+i*6+2]=Calc.unitVector( Calc.add(bc, ca) );
            l_rays[5+i*6+3]=Calc.unitVector( Calc.add(ca, ab) );

            a=l_rays[5];

            ab=Calc.add(a, b);
            CalcGeom.product2(ab, 0.5);

            bc=Calc.add(b, c);
            CalcGeom.product2(bc, 0.5);

            ca=Calc.add(c, a);
            CalcGeom.product2(ca, 0.5);

            l_rays[5+i*6+4]=Calc.unitVector( Calc.add(ab, bc) );
            l_rays[5+i*6+5]=Calc.unitVector( Calc.add(bc, ca) );
            l_rays[5+i*6+6]=Calc.unitVector( Calc.add(ca, ab) );

        }

        return l_rays;
    }
    
    /**
     *
     * @param btree
     * @param ffield
     * @param actSiteCenter
     * @param side_length
     * @param grid
     * @return
     */
    private static AtomGrid generateGridActSite(BSPTree btree, IForceField ffield,
            Atom actSiteCenter, double side_length, AtomGrid grid) {

        
        double gridResol=GRID_RESOL;

        int npointx=(int) Math.ceil( side_length / gridResol);
        int npointy=(int) Math.ceil( side_length / gridResol);
        int npointz=(int) Math.ceil( side_length / gridResol);

        //adjust grid resolution if the dimension is bigger
        while( npointx*npointy*npointz > grid.size() ){

            gridResol+=GRID_RESOL_INC;
            npointx=(int) Math.ceil( side_length / gridResol);
            npointy=(int) Math.ceil( side_length / gridResol);
            npointz=(int) Math.ceil( side_length / gridResol);
        }


        grid.clear();
        grid.setDimension(npointx,npointy,npointz);
        grid.setResol(gridResol);

        double minX=actSiteCenter.getX() - side_length*0.5;
        double minY=actSiteCenter.getY() - side_length*0.5;
        double minZ=actSiteCenter.getZ() - side_length*0.5;
        Atom point=new AtomImpl();

        // construct hte grid points
        for(int i=0;i<npointx;i++){
            for(int j=0;j<npointy;j++){
                for(int k=0;k<npointz;k++){

                    int idx=grid.getIndex(i, j, k);

                    point.setX( minX+i*gridResol);
                    point.setY( minY+j*gridResol);
                    point.setZ( minZ+k*gridResol);

                    grid.setPoint(idx, point, false);
                }
            }
        }

        //fill the grid starting in the center, borders stop the propagation
        fillExterior( grid, btree, ffield, actSiteCenter );

        int lim=npointx*npointy*npointz;

        //assign points not in solvent inside atoms volume
        for(int i=0;i<lim;i++){

            double val=grid.getValue(i);

            if( val!=BORDER_VALUE && val!=SOLVENT_VALUE ){
                grid.setValue( i, BORDER_VALUE);
                grid.setInVolume( i, true);
            }
        }

        return grid;
    }


    /**
     * calculates the pocket inside a box given a center and a side length
     * 
     * @param btree
     * @param actSiteCenter
     * @param side_length
     * @param grid
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static GridProbeCluster calcActiveSite( BSPTree btree,
            Atom actSiteCenter, double side_length, AtomGrid grid )
            throws StructureException {
        
        IForceField ffield = new TestForceField();

        //generate the grid
        generateGridActSite( btree, ffield, actSiteCenter, side_length, grid);

        // get the grid points that are outside the VdW volume
        int lim=grid.getNpointx()*grid.getNpointy()*grid.getNpointz();

        ArrayList<GridProbe> gridProbes=new ArrayList<GridProbe>();

        //gets the grid atoms of the solvent
        for(int i=0;i<lim;i++){

            Atom p=grid.getPoint(i);

            if( !grid.isInVolume(i) ){
                gridProbes.add( new GridProbe(i, p, 0) );
            }
        }


        // Calculation of buriedness values of grid probes 
        calcBuriedIndices( gridProbes, btree );

        //cluster the probes
        ArrayList<GridProbeCluster> clusters =
                calcClusters( gridProbes, grid, LOW_LIMIT, 30);

        if( clusters==null )
            return null;

        //get the cluster nearest to center
        GridProbeCluster ncluster=null;
        double mindist=java.lang.Double.MAX_VALUE;

        for( GridProbeCluster cluster : clusters ){
            double dist=CalcGeom.squareDistance(actSiteCenter, cluster.getCentroid() );

            if( dist < mindist ){
                mindist=dist;
                ncluster=cluster;
            }
        }

        return ncluster;
    }

    /**
     *
     * @param btree
     * @param actSiteCenter
     * @param radii
     * @param grid
     * @return
     * @throws org.biojava.bio.structure.StructureException
     */
    public static GridProbeCluster generateActiveSite( BSPTree btree,
            Atom actSiteCenter, double radii, AtomGrid grid )
            throws StructureException {

        IForceField ffield = new TestForceField();

        //generate the grid
        generateGridActSite( btree, ffield, actSiteCenter, radii*2.0, grid);

        // get the grid points that are outside the VdW volume
        int lim=grid.getNpointx()*grid.getNpointy()*grid.getNpointz();


        GridProbeCluster cluster=new GridProbeCluster();

        //gets the grid atoms of the solvent
        for(int i=0;i<lim;i++){

            Atom p=grid.getPoint(i);

            if( !grid.isInVolume(i) ){
                cluster.getProbes().add( new GridProbe(i, p, 0) );
            }
        }

        return cluster;
    }

    
}
