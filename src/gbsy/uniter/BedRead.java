/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy.uniter;

/**
 *
 * @author koen
 */
    
public class BedRead{

    private final String chromosome;
    private final String beginpos;
    private final String endpos;

    public BedRead(String chromosome, String beginPos, String endPos){
        this.chromosome = chromosome;
        this.beginpos = beginPos;
        this.endpos = endPos;
    }

    public String getChromosome(){
        return this.chromosome;
    }

    public String getBeginPos(){
        return this.beginpos;
    }

    public String getEndPos(){
        return this.endpos;
    }
}
