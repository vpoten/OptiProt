/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io.mol2;

import org.jumpmind.symmetric.csv.CsvReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author victor
 */
public class Elements {

    private static ArrayList<ElementEntry> elements=null;
    private static HashMap<String, ElementEntry> table=null;

    static{
        try {
            loadElements("org/optiprot/data/elements.csv");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Elements.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Elements.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static void loadElements( String filename )
            throws FileNotFoundException, IOException{

        elements=new ArrayList<ElementEntry>();
        table=new HashMap<String, ElementEntry>();

        InputStream fileStream=Elements.class.getClassLoader().getResourceAsStream(filename);

        CsvReader reader=new CsvReader(fileStream, Charset.forName("UTF8"));

        reader.readHeaders();

        while (reader.readRecord() ){

            int Num=Integer.parseInt(reader.get(0));
            String Symb=reader.get(1);
            double ARENeg=Double.parseDouble(reader.get(2));
            double RCov=Double.parseDouble(reader.get(3));
            double RBO=Double.parseDouble(reader.get(4));
            double RVdW=Double.parseDouble(reader.get(5));
            int MaxBnd=Integer.parseInt(reader.get(6));
            double Mass=Double.parseDouble(reader.get(7));
            double ElNeg=Double.parseDouble(reader.get(8));
            double Ionization=Double.parseDouble(reader.get(9));
            double ElAffinity=Double.parseDouble(reader.get(10));
            double Red=Double.parseDouble(reader.get(11));
            double Green=Double.parseDouble(reader.get(12));
            double Blue=Double.parseDouble(reader.get(13));
            String Name = reader.get(14);

            ElementEntry ent=new ElementEntry( Num,  Symb,  ARENeg,  RCov,
             RBO,  RVdW,  MaxBnd, Mass,  ElNeg,
             Ionization, ElAffinity,
             Red,  Green, Blue,  Name);

            elements.add( ent );
            table.put(Symb.toUpperCase(), ent);

        }

        reader.close();
        
    }

    static public ElementEntry getElement( String symb ){
        return table.get(symb.toUpperCase());
    }

    static public ElementEntry getElement( int num ){
        return elements.get(num);
    }
}
