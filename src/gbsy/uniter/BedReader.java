/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy.uniter;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 *
 * @author koen
 */

public class BedReader{

    private final BufferedReader bedBufferedReader;

    public BedReader(String bedFile) throws FileNotFoundException{
        this.bedBufferedReader = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(new File(bedFile)))));
    }

    public BedRead getNextRead(){
        try {
            String line = this.bedBufferedReader.readLine();
            if (line == null){
                return null;
            }
            String[] linesplit = line.split("\t");
            BedRead bedRead = new BedRead(linesplit[0], linesplit[1], linesplit[2]);
            return bedRead;
        } catch (IOException ex) {
            return null;
        }
    }

    public void close() throws IOException{
        this.bedBufferedReader.close();
    }

}
    