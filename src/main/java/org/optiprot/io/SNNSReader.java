/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.optiprot.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.StringTokenizer;
import org.optiprot.neural.NNPattern;

/**
 * SNNS pattern (neural network) file reader/writer
 * 
 * @author victor
 */
public class SNNSReader {

    //SNNS file fields
    final static private String HEADER = "SNNS pattern definition file V3.2";
    final static private String DATE = "generated at ";
    final static private String NPATS = "No. of patterns : ";
    final static private String NINPUTS = "No. of input units : ";
    final static private String NOUTS = "No. of output units : ";
    final static private String COMMENT = "#";
    

    
    static public ArrayList<NNPattern> parseSNNSFile( String path )
            throws FileNotFoundException, IOException {

        BufferedReader r = new BufferedReader( new FileReader(path) );

        //read SNNS file header

        String l_line=r.readLine();

        if( !l_line.startsWith(HEADER) )
            throw new IOException("Bad SNNS header");

        l_line=r.readLine();

        if( !l_line.startsWith(DATE) )
            throw new IOException("Bad SNNS header");

        l_line=r.readLine();
        l_line=r.readLine();
        l_line=r.readLine();

        if( !l_line.startsWith(NPATS) )
            throw new IOException("Bad SNNS header");

        int npatterns = Integer.parseInt( l_line.substring(NPATS.length()) );

        l_line=r.readLine();

        if( !l_line.startsWith(NINPUTS) )
            throw new IOException("Bad SNNS header");

        int ninputs = Integer.parseInt( l_line.substring(NINPUTS.length()) );

        l_line=r.readLine();

        int noutputs = 0;

        if( l_line.startsWith(NOUTS) ){
            noutputs = Integer.parseInt( l_line.substring(NOUTS.length()) );
            l_line=r.readLine();
        }

        ArrayList<NNPattern> list=new ArrayList<NNPattern>();

        double [] buffer = new double [ninputs+noutputs];
        int pos=0;

        String [] tokens = new String [30];
        

        //read the rest of the file
        while( l_line!=null )
        {
            l_line=l_line.trim();

            if( l_line.isEmpty() ){
                //nothing to do
            }
            else if( l_line.startsWith(COMMENT) ){
               // a comment
            }
            else{
                int ntokens = getTokens( l_line, " \t", tokens.length, tokens);

                for(int i=0; i<ntokens; i++)
                    buffer[pos++] = Double.parseDouble( tokens[i] );

                if( pos==buffer.length ){
                    pos=0;
                    list.add( new NNPattern( buffer, ninputs, noutputs) );
                }
            }

            l_line=r.readLine();
        }

        r.close();

        return list;
    }


    static public void writeSNNSFile( String path, ArrayList<NNPattern> list )
            throws FileNotFoundException, IOException {

        BufferedWriter w = new BufferedWriter( new FileWriter(path) );

        //write header
        w.write( HEADER+"\n" );

        Calendar cal=Calendar.getInstance();
        String date=String.format("%1$tY-%1$tm-%1$td %1$tH:%1$tM:%1$tS", cal);
        w.write( DATE+date+"\n" );
        w.newLine();
        w.newLine();
        w.write( NPATS+list.size()+"\n" );
        w.write( NINPUTS+list.get(0).getNInputs()+"\n" );
        w.write( NOUTS+list.get(0).getNOutputs()+"\n" );
        w.newLine();

        //write patterns
        for( NNPattern pat : list ){

            String line="";

            double [] vector=null;

            vector=pat.getInputs();

            //write inputs
            for( int i=0; i<vector.length; i++){
                line+=vector[i]+" ";
            }

            w.newLine();

            line="";
            vector=pat.getOutputs();

            //write outputs
            for( int i=0; i<vector.length; i++){
                line+=vector[i]+" ";
            }

            w.newLine();
            w.newLine();
        }

        w.close();
    }


    /**
     *
     * @param line : line to parse
     * @param delimit : delimiters
     * @param num : number of tokens to obtain
     * @param tokens : array of tokens I/O
     * @return number of readed tokens
     */
    static private int getTokens( String line, String delimit, int num, String [] tokens ){

        StringTokenizer l_tokenizer = new StringTokenizer(line, delimit);

        int i=0;

        while (l_tokenizer.hasMoreTokens())
        {
            tokens[i++] = l_tokenizer.nextToken();

            if( i==(num-1) )
                break;
        }

        return i;
    }

}
