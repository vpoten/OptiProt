/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.molsurface;

import java.util.ArrayDeque;
import java.util.ArrayList;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.optiprot.maths.BSPTree;
import org.optiprot.maths.CalcGeom;
import org.optiprot.maths.CalcSphere;
import org.optiprot.potential.IForceField;
import org.optiprot.potential.TestForceField;

/**
 *
 * @author victor
 */
public class ReducedSurface {
    
    private ArrayList<RSVertex> vertices=new ArrayList<RSVertex>();
    private ArrayList<RSEdge> edges = new ArrayList<RSEdge>();
    private ArrayList<RSFace> faces = new ArrayList<RSFace>();

    static final private double PROBE_RADII=1.4;
    static final private double MAX_ARADII=2.0;

    static final private double _2_PI=2.0*Math.PI;


    public ReducedSurface() {
    }

    public void compute( Atom [] atoms , BSPTree btree) throws StructureException {

        RSFace face=initialFace( atoms, btree, PROBE_RADII );

        this.getFaces().add(face);
        this.getEdges().add( face.getEdge1() );
        this.getEdges().add( face.getEdge2() );
        this.getEdges().add( face.getEdge3() );
        this.getVertices().add( face.getVertex1() );
        this.getVertices().add( face.getVertex2() );
        this.getVertices().add( face.getVertex3() );

        rsComponent( face, btree );

        // TODO - connect other RSComponents
        
    }

    
    /**
     * 
     * @param face : initial face to extend
     * @param btree
     * @throws org.biojava.bio.structure.StructureException
     */
    private void rsComponent( RSFace face, BSPTree btree ) throws StructureException{

        ArrayDeque<RSEdge> fifo=new ArrayDeque<RSEdge>();
        ArrayList<Atom> list=new ArrayList<Atom>();
        double [] distances=new double [3];
        IForceField ffield = new TestForceField();

        fifo.add( face.getEdge1() );
        fifo.add( face.getEdge2() );
        fifo.add( face.getEdge3() );

        RSEdge currEdge=null;

        while( !fifo.isEmpty() ){

            currEdge=fifo.poll();

            // treatment of the edge

            list.clear();

            //get  the atoms closests to the edge
            getAtomsCloseEdge( currEdge.getVertex1(), currEdge.getVertex2(),
                    btree, PROBE_RADII, list );


            //calcs the center of the arc trajectory of probe ball
            Atom ab=Calc.substract( currEdge.getVertex2(), currEdge.getVertex1());

            double d = Calc.skalarProduct(
                    Calc.unitVector(
                    Calc.substract( currEdge.getFace1().getCenter(), currEdge.getVertex1())
                    ),
                    ab);

            Atom centArc = Calc.add( currEdge.getVertex1(),
                    CalcGeom.product(Calc.unitVector(ab), d) );

            Atom mc=Calc.unitVector(Calc.substract( currEdge.getFace1().getCenter(), centArc ));


            //get the vertex of the face that not is in the edge
            Atom c1=currEdge.getFace1().getOtherVertex(currEdge);
            Atom mc1=Calc.substract(c1, centArc);
            Atom rotAx=Calc.vectorProduct( mc1, mc);


            Atom newAt=null;
            Atom centerProbe=null;
            Atom orgProbe=null;
            double ang=_2_PI;


            for(Atom at : list ){

                if( at.equals(c1) )
                    continue;

                distances[0]=Calc.getDistance( currEdge.getVertex1(), currEdge.getVertex2());
                distances[1]=Calc.getDistance( currEdge.getVertex2(), at);
                distances[2]=Calc.getDistance( at, currEdge.getVertex1());

                Atom [] centers=null;

                //calc the centers of the probe sphere tangent to the three atoms
                centers = CalcSphere.trilateration(
                        currEdge.getVertex1(), PROBE_RADII+ffield.getVDWRadius(currEdge.getVertex1()),
                        currEdge.getVertex2(), PROBE_RADII+ffield.getVDWRadius(currEdge.getVertex2()),
                        at, PROBE_RADII+ffield.getVDWRadius(at), distances );

                Atom v1=Calc.unitVector( Calc.substract(centers[1], centArc) );
                Atom v2=Calc.substract(centers[2], centArc);

                double ang1=Calc.skalarProduct( v1, mc);
                double ang2=Calc.skalarProduct( v2, mc);

                if( Calc.skalarProduct( Calc.vectorProduct(mc,v1), rotAx)<0.0 )
                    ang1=_2_PI-ang1;

                if( Calc.skalarProduct( Calc.vectorProduct(mc,v2), rotAx)<0.0 )
                    ang2=_2_PI-ang2;

                //get the center with smallest arc

                if( ang1<ang ){
                    ang=ang1;
                    centerProbe=centers[1];
                    orgProbe=centers[0];
                    newAt=at;
                }

                if( ang2<ang ){
                    ang=ang2;
                    centerProbe=centers[2];
                    orgProbe=centers[0];
                    newAt=at;
                }

            }

            if( newAt!=null ){
                //creates a new face and adds the new edges to the fifo

                Atom dir=Calc.substract(centerProbe, orgProbe);

                face=new RSFace();

                face.setCenter(centerProbe);

                face.setVertex1( currEdge.getVertex1() );
                face.setVertex2( currEdge.getVertex2() );
                face.setVertex3( new RSVertex(newAt,ffield.getVDWRadius(newAt)) );

                Atom normal=CalcGeom.normalPlane( face.getVertex1(),
                        face.getVertex2(), face.getVertex3() );

                //calcs normal to the outside
                if( Calc.skalarProduct(normal, dir)<0.0 )
                    CalcGeom.product2(normal, -1.0);

                face.setNormal( normal );

                currEdge.setFace2( face );

                RSEdge edge2=new RSEdge();
                edge2.setFace1(face);
                edge2.setVertex1( face.getVertex2() );
                edge2.setVertex2( face.getVertex3() );

                RSEdge edge3=new RSEdge();
                edge3.setFace1(face);
                edge3.setVertex1( face.getVertex3() );
                edge3.setVertex2( face.getVertex1() );

                face.setEdge1( currEdge );
                face.setEdge2(edge2);
                face.setEdge3(edge3);

                // adds the new edges to the fifo
                fifo.add( face.getEdge2() );
                fifo.add( face.getEdge3() );

                //add the new elements
                this.getFaces().add(face);
                this.getEdges().add( face.getEdge2() );
                this.getEdges().add( face.getEdge3() );
                this.getVertices().add( face.getVertex3() );

            }//end if

        }

    }


    /**
     * @return the vertices
     */
    public ArrayList<RSVertex> getVertices() {
        return vertices;
    }

    /**
     * @return the edges
     */
    public ArrayList<RSEdge> getEdges() {
        return edges;
    }

    /**
     * @return the faces
     */
    public ArrayList<RSFace> getFaces() {
        return faces;
    }

    private RSFace initialFace( Atom [] atoms, BSPTree btree, double probe )
            throws StructureException{

        int axis=0;
        Atom leftAtom=atoms[0];

        IForceField ffield = new TestForceField();

        ArrayList<Atom> list=new ArrayList<Atom>();

        double [] distances=new double [3];

        //probe the 3 axis
        for( axis=0;axis<3;axis++){

            //get the left-most atom
            for(Atom at : atoms ){
                if( at.getCoords()[axis]<leftAtom.getCoords()[axis] ){
                    leftAtom=at;
                }
            }

            //get the nearests atoms from leftAtom
            btree.neighbours(leftAtom, MAX_ARADII*2.0, list);

            if( list.isEmpty() )
                continue;

            Atom atom2=null;

            //get the left most atom of the neighbors
            for(Atom at : list){

                if( at.equals(leftAtom) )
                    continue;

                if( atom2==null )
                    atom2=(RSVertex) at;

                if( at.getCoords()[axis]<atom2.getCoords()[axis] ){
                    atom2=(RSVertex) at;
                }
            }

            list.clear();

            getAtomsCloseEdge( leftAtom, atom2, btree, probe, list );

            //get the third atom of the face
            for(Atom atom3 : list){

                if( atom3.equals(leftAtom) || atom3.equals(atom2) )
                    continue;

                distances[0]=Calc.getDistance(leftAtom, atom2);
                distances[1]=Calc.getDistance(atom2, atom3);
                distances[2]=Calc.getDistance(atom3, leftAtom);

                Atom [] centers=null;

                //calc the centers of the probe sphere tangent to the three atoms
                centers = CalcSphere.trilateration(
                        leftAtom, probe+ffield.getVDWRadius(leftAtom),
                        atom2, probe+ffield.getVDWRadius(atom2),
                        atom3, probe+ffield.getVDWRadius(atom3), distances );

                boolean inter1 = false;
                boolean inter2 = false;

                // check if is a collision free probe
                for(Atom at : list){

                    if( at.equals(leftAtom) || at.equals(atom2) ||
                            at.equals(atom3) )
                        continue;

                    if( CalcSphere.isIntersection(at, ffield.getVDWRadius(at), centers[1], probe) )
                        inter1=true;

                    if( CalcSphere.isIntersection(at, ffield.getVDWRadius(at), centers[2], probe) )
                       inter2=true;

                    if( inter1 && inter2 )
                        break;
                }

                Atom centerProbe=null;

                if( !inter1 && !inter2 ){

                    centerProbe=centers[1];

                    if( centers[2].getCoords()[axis] <
                            centers[1].getCoords()[axis] ){
                        centerProbe=centers[2];
                    }

                }
                else if( !inter1 ){
                    centerProbe=centers[1];
                }
                else if( !inter2 ){
                    centerProbe=centers[2];
                }


                if( centerProbe!=null ){
                    //constructs and return the face

                    Atom dir=Calc.substract(centerProbe, centers[0]);

                    RSFace face=new RSFace();
                    
                    face.setCenter(centerProbe);

                    face.setVertex1( new RSVertex(leftAtom,ffield.getVDWRadius(leftAtom)) );
                    face.setVertex2( new RSVertex(atom2,ffield.getVDWRadius(atom2)) );
                    face.setVertex3( new RSVertex(atom3,ffield.getVDWRadius(atom3)) );

                    Atom normal=CalcGeom.normalPlane(leftAtom, atom2, atom3);

                    //calcs normal to the outside
                    if( Calc.skalarProduct(normal, dir)<0.0 )
                        CalcGeom.product2(normal, -1.0);

                    face.setNormal( normal );

                    RSEdge edge1=new RSEdge();
                    edge1.setFace1(face);
                    edge1.setVertex1( face.getVertex1() );
                    edge1.setVertex2( face.getVertex2() );

                    RSEdge edge2=new RSEdge();
                    edge2.setFace1(face);
                    edge2.setVertex1( face.getVertex2() );
                    edge2.setVertex2( face.getVertex3() );

                    RSEdge edge3=new RSEdge();
                    edge3.setFace1(face);
                    edge3.setVertex1( face.getVertex3() );
                    edge3.setVertex2( face.getVertex1() );

                    face.setEdge1(edge1);
                    face.setEdge2(edge2);
                    face.setEdge3(edge3);

                    return face;

                }//
            }
            
        }

        return null;
    }

    /**
     * get the atoms close enough to the edge at1-at2
     * 
     * @param at1
     * @param at2
     * @param btree
     * @param probe
     * @param list
     */
    private void getAtomsCloseEdge(Atom at1, Atom at2, BSPTree btree,
            double probe, ArrayList<Atom> list ){

        Atom midpoint = Calc.add(at1, at2);
        CalcGeom.product2( midpoint, 0.5);

        //get atoms close enough to the edge at1-at2
        btree.neighbours( midpoint, MAX_ARADII+probe*2.0, list);
    }
    
}
