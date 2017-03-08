/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy.uniter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

/**
 *
 * @author koen
 */
public class MyWriter{
        
    private Writer logwriter;

    public MyWriter(String outputDir, String fileName) throws IOException{
        if (outputDir != null){
            File file = new File(outputDir + System.getProperty("file.separator") + fileName);
            file.createNewFile();
            this.logwriter = new BufferedWriter(new FileWriter(file));
        }else{
            this.logwriter = null;
        }
    }

    public void addToFile(String line) throws IOException{
        if (this.logwriter != null){
            this.logwriter.write(line);
        }else{
            System.out.print("" + line);
        }
    }

    public void closeWriter() throws IOException{
        if (this.logwriter != null){
            this.logwriter.close();
            this.logwriter = null;
        }
    }

    @Override
    public void finalize() throws Throwable{
        if (this.logwriter != null){
            this.logwriter.close();
            this.logwriter = null;
        }
        super.finalize();
    }
}
