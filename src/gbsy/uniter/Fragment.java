/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy.uniter;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 *
 * @author koen
 */

public class Fragment implements Comparable<Fragment>{

    private HashMap<String, Integer> sample_count_map;
    private final String chromosome;
    private final int begin;
    private final int end;
    private final int enzymeOverlap;

    public Fragment(String[] samples, String chromosome, int begin, int end, int enzymeOverlap){
        sample_count_map = new HashMap<>();
        for (String s : samples){
            sample_count_map.put(s, 0);
        }
        this.chromosome = chromosome;
        this.begin = begin;
        this.end = end;
        this.enzymeOverlap = enzymeOverlap;
    }

    public String[] getSamples(){
        String[] samples = new String[this.sample_count_map.keySet().size()];
        int index = 0;
        for (String s : this.sample_count_map.keySet()){
            samples[index] = s;
            index++;
        }
        return samples;
    }

    public String getChromosome(){
        return this.chromosome;
    }

    public int getBegin(){
        return this.begin;
    }

    public int getEnd(){
        return this.end;
    }

    public int getLength(){
        return this.end - this.begin;
    }

    public int getEnzymeOverlap(){
        return this.enzymeOverlap;
    }

    public void addSampleCount(String sample, int count){
        this.sample_count_map.put(sample, count + this.sample_count_map.get(sample));
    }

    public int getCounts(String sample){
        return this.sample_count_map.get(sample);
    }

    public double getCounts(String sample, Map<String, Double> sampleNormalizationFactorMap){
        return this.sample_count_map.get(sample) * sampleNormalizationFactorMap.get(sample);
    }

    public double getSizeFactor(String sample){
        double product = 1.0;
        int sample_number = 0;
        for (String s : this.sample_count_map.keySet()){
            if (this.sample_count_map.get(s) != 0){
                product = product * this.sample_count_map.get(s);
                sample_number++;
            }
        }
        double mean = Math.pow((product * 1.0), (1/(sample_number * 1.0)));
        double sizefactor = (this.sample_count_map.get(sample) * 1.0)/(mean *1.0);
        return sizefactor;
    }

    public double getSizeFactor(String sample, HashMap<String,Double> sample_library_size_map){
        double product = 1.0;
        int sample_number = 0;
        for (String s : this.sample_count_map.keySet()){
            if (this.sample_count_map.get(s) != 0){
                product = product * (this.sample_count_map.get(s) * sample_library_size_map.get(s));
                sample_number++;
            }
        }
        double mean = Math.pow((product * 1.0), (1/(sample_number * 1.0)));
        double sizefactor = (this.sample_count_map.get(sample) * 1.0)/(mean *1.0);
        return sizefactor;
    }

    public int getOccuranceSamples(){
        int sample_occurance = 0;
        for (String sample : this.sample_count_map.keySet()){
            if (this.sample_count_map.get(sample) != 0){
                sample_occurance++;
            }
        }
        return sample_occurance;
    }

    public HashSet<String> getOccuranceSamplesList(){
        HashSet<String> samples = new HashSet<>();
        for (String sample : this.sample_count_map.keySet()){
            if (this.sample_count_map.get(sample) != 0){
                samples.add(sample);
            }
        }
        return samples;
    }

    public int getOccuranceTotal(){
        int total_occurance = 0;
        for (String sample : this.sample_count_map.keySet()){
            total_occurance += this.sample_count_map.get(sample);
        }
        return total_occurance;
    }

    public boolean filter(int minSamples, double minOccurancesPerSample){
        boolean filter = true;
        if (this.getOccuranceSamples() < minSamples){
            filter = false;
        }else{
            int countMinOccuranceSamples = 0;
            for (String sample : this.sample_count_map.keySet()){
                if(this.sample_count_map.get(sample) >= minOccurancesPerSample){
                    countMinOccuranceSamples++;
                }
            }
            if (countMinOccuranceSamples < minSamples){
                filter = false;
            }
        }
        return filter;
    }

    public boolean delete(Map<String, Double> sampleNormalizationFactorMap, double cutoff){
        boolean delete = true;
        for (String s : sampleNormalizationFactorMap.keySet()){
            if (this.sample_count_map.containsKey(s)){
                if (this.sample_count_map.get(s) * sampleNormalizationFactorMap.get(s) > cutoff){
                    delete = false;
                }
            }
        }
        return delete;
    }

//        @Override
    public int compareTo2(Fragment o) {
        if (this.chromosome.compareTo(o.getChromosome()) != 0){
            return this.chromosome.compareTo(o.getChromosome());
        }else{
            if (this.getBegin() == o.getBegin() && this.getEnd() == o.getEnd()){
                return 0;
            }else{
                if (this.getBegin() + this.getEnzymeOverlap() > o.getEnd() - o.getEnzymeOverlap()){
                    //o is before this
                    return -1;
                }
                if (this.getEnd() - this.getEnzymeOverlap() < o.getBegin() + o.getEnzymeOverlap()){
                    //o is after this
                    return 1;
                }
                if (this.getBegin() + this.getEnzymeOverlap() < o.getBegin() - this.getEnzymeOverlap()){
                    //o is overlapping, but this was starting first
                    return -1;
                }
                if (this.getBegin() + this.getEnzymeOverlap() > o.getBegin() - o.getEnzymeOverlap()){
                    //o is overlapping, and o was starting first
                    return 1;
                }
                if (this.getBegin() + this.getEnzymeOverlap() == o.getBegin() + o.getEnzymeOverlap()){
                    //both starting on the same place
                    if (this.getEnd() - this.getEnzymeOverlap() < o.getEnd() - o.getEnzymeOverlap()){
                        //this ends first
                        return 1;
                    }else{
                        //other ends first
                        return -1;
                    }
                }
            }
        }
        return 1;
    }

    @Override
    public int compareTo(Fragment o) {
        if (this.chromosome.compareTo(o.getChromosome()) != 0){
            return this.chromosome.compareTo(o.getChromosome());
        }else{
            if (this.getBegin() == o.getBegin() && this.getEnd() == o.getEnd()){
                return 0;
            }else{
                if (this.getBegin() > o.getEnd()){
                    //o is before this
                    return 1;
                }
                if (this.getEnd() < o.getBegin()){
                    //o is after this
                    return -1;
                }
                //Overlapping
                if (this.getBegin() < o.getBegin()){
                    //o is overlapping, but this was starting first
                    return -1;
                }
                if (this.getBegin() == o.getBegin()){
                    //same start
                    if (this.getEnd() < o.getEnd()){
                        return 1;
                    }else{
                        return -1;
                    }
                }
            }
        }
        return 1;
    }

    public boolean canCombine(Fragment f, int enzymeOverlap){
        if (!this.getChromosome().equals(f.getChromosome())){
            return false;
        }
        if ((this.getBegin() == f.getBegin() || this.getBegin() + enzymeOverlap == f.getBegin() || this.getBegin() - enzymeOverlap == f.getBegin()) && 
                (this.getEnd() == f.getEnd() || this.getEnd() + enzymeOverlap == f.getEnd() || this.getEnd() - enzymeOverlap == f.getEnd())){
            return true;
        }
        return false;
    }
}
