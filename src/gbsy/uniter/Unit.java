/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy.uniter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;

/**
 *
 * @author koen
 */

public class Unit implements Comparable<Unit>{

    private HashSet<Fragment> fragment_set;
    private final String chromosome;
    private int begin;
    private int end;
    private final int enzymeOverlap;

    public Unit(Fragment f){
        this.fragment_set = new HashSet<>();
        this.fragment_set.add(f);
        this.chromosome = f.getChromosome();
        this.begin = f.getBegin();
        this.end = f.getEnd();
        this.enzymeOverlap = f.getEnzymeOverlap();
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

    public int getEnzymeOverlap(){
        return this.enzymeOverlap;
    }

    public ArrayList<Fragment> getFragments(){
        ArrayList<Fragment> fragmentList = new ArrayList<>(this.fragment_set);
        Collections.sort(fragmentList, new FragmentComparator());
        return fragmentList;
    }

    boolean canAddFragment(Fragment f){
        if (!f.getChromosome().equals(this.chromosome)){
            //f is on other chromosome
            return false;
        }
        if(f.getEnd() <= this.getBegin() + this.enzymeOverlap){
            //f ends before this begins
            return false;
        }
        if(f.getBegin() + this.enzymeOverlap > this.getEnd()){
            //f begins after this ends
            return false;
        }
        return true;
    }

    public boolean isSimpleUnit(){
        return (this.fragment_set.size() == 1);
    }

    public int getNumberOfFragments(){
        return this.fragment_set.size();
    }

    public int getUnitDepth(){
        int depth = 0;
        for (Fragment f : this.fragment_set){
            depth += f.getOccuranceTotal();
        }
        return depth;
    }

    public int getUnitDepth(String sample){
        int depth = 0;
        for (Fragment f : this.fragment_set){
            depth += f.getCounts(sample);
        }
        return depth;
    }

    public int getNumberOfSamples(){
        HashSet<String> sampleset = new HashSet<String>();
        for (Fragment f : this.fragment_set){
            sampleset.addAll(f.getOccuranceSamplesList());                
        }            
        return sampleset.size();
    }

    public void addFragment(Fragment f){
        this.fragment_set.add(f);
        if (this.begin > f.getBegin()){
            this.begin = f.getBegin();
        }
        if (this.end < f.getEnd()){
            this.end = f.getEnd();
        }
    }

    public double getSizeFactor(String sample){
        double product = 1.0;
        int sample_number = 0;
        for (String s : this.fragment_set.iterator().next().getSamples()){
            if (this.getUnitDepth(s) != 0){
                product = product * this.getUnitDepth(s);
                sample_number++;
            }
        }
        double mean = Math.pow((product * 1.0), (1/(sample_number * 1.0)));
        double sizefactor = (this.getUnitDepth(sample) * 1.0)/(mean *1.0);
        return sizefactor;
    }

    public boolean delete(Map<String, Double> sampleNormalizationFactorMap, double cutoff){
        HashSet<Fragment> todelete = new HashSet<>();
        for (Fragment f : this.fragment_set){
            if (f.delete(sampleNormalizationFactorMap, cutoff)){
                todelete.add(f);
            }
        }
        for (Fragment f : todelete){
            this.fragment_set.remove(f);
        }
        return this.fragment_set.isEmpty();
    }

    public void combineFragments(int enzymeOverlap){
        boolean redo = false;
        do{
            ArrayList<Fragment> fragmentList = new ArrayList<>(this.fragment_set);
            Collections.sort(fragmentList, new FragmentComparator());
            redo = false;
            for(int i = fragmentList.size() -1; ! redo && i > 0; i--){
                if (fragmentList.get(i).canCombine(fragmentList.get(i - 1), enzymeOverlap)){
                    //can combine => so combine
                    redo = true;
                    int begin = fragmentList.get(i).getBegin();
                    if (fragmentList.get(i - 1).getBegin() < begin){
                        begin = fragmentList.get(i -1).getBegin();
                    }
                    int end = fragmentList.get(i).getEnd();
                    if (fragmentList.get(i - 1).getEnd()< end){
                        end = fragmentList.get(i -1).getEnd();
                    }
                    Fragment f = new Fragment(fragmentList.get(i -1).getSamples(), fragmentList.get(i -1).getChromosome(), begin, end, enzymeOverlap);
                    for (String sample : fragmentList.get(i).getSamples()){
                        f.addSampleCount(sample, fragmentList.get(i).getCounts(sample));
                    }
                    for (String sample : fragmentList.get(i -1).getSamples()){
                        f.addSampleCount(sample, fragmentList.get(i-1).getCounts(sample));
                    }
                    this.fragment_set.remove(fragmentList.get(i -1));
                    this.fragment_set.remove(fragmentList.get(i));
                    this.fragment_set.add(f);                    
                }
            }
        }while(redo);
    }


    @Override
    public int compareTo(Unit o) {
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

}
    
