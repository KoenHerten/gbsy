/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy.uniter;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the Genomics Core Leuven
 */
public class Uniter {
    
    public static void main(String[] args){
        String input = "-b /home/koen/Downloads/bed_test/lib1_DLT_1.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_DLT_3.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_GRSP_2.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_PK_1.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_STY_1.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_STY_3.bed.LG10.bed "
                + "-b /home/koen/Downloads/bed_test/lib1_DLT_2.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_GRSP_1.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_GRSP_3.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_PK_3.bed.LG10.bed  "
                + "-b /home/koen/Downloads/bed_test/lib1_STY_2.bed.LG10.bed "
                + "-o /home/koen/Downloads/bed_test/result "
                + "-mo 5 ";
        args = input.split(" ");
        Uniter bedToCount = new Uniter(args);
    }
    
    public Uniter(String[] args){
        this.outputDir = null;
        this.minimum_occurance = 0;
        this.minimum_samples = 1;
        this.enzymeOverlap = 0;
        this.bedFiles = new ArrayList<>();
        this.chromosome_location_sample_countMap = new HashMap<>();
        this.chromosome_sample_countMap = new HashMap<>();
        this.filePrefix = "";
        double minimum_samples_percentage = 0;
        if (args.length == 0){
            args = new String[1];
            args[0] = "-h";
        }
        for (int i=0; i < args.length; i++){
            if (args[i].equals("-h") || args[i].equals("-help")){
                System.out.println("Transforms bed files to Genomic Count files");
                System.out.println("\t-b\tadd a bed file");
                System.out.println("\t-o\tchange the output from command line to files (given directory)");
                System.out.println("\t-mo\tminimum occurances of a location in the locations file (standard 0)");
                System.out.println("\t-ms\tminimum number of samples for a location in the locations file (standard 0)");
                System.out.println("\t-msp\tminimum number of samples for a location in the locations file as percentage (standard 0, max 1)"
                        + "\t\tThe highest value of -ms and -msp is taken (-msp expressed in number of samples)");
                System.out.println("\t-eo\tthe length of the overlap of the enzyme (will be used to find repeat enzyme sites) "
                        + "(standard 0, should be the overlap of the enzyme)");
                System.out.println("\t-w\tweigthed (true) or absolute counts (false) (standard false)");
                System.out.println("\t-r\tround: set values lower then the minimum to 0 (true) or keep these values (false) (standard false)");
                System.out.println("\t-prefix\tthe name of the comparison (added in front of the outputfiles)");
                System.exit(0);
            }
            if (args[i].equals("-b")){
                //bed file
                this.bedFiles.add(args[i+1]);
                i++;
            }else if(args[i].equals("-o")){
                //output
                this.outputDir = args[i+1];
                i++;
            }else if(args[i].equals("-mo")){
                this.minimum_occurance = Double.parseDouble(args[i+1]);
                if (this.minimum_occurance < 0){this.minimum_occurance = 0;}
                i++;
            }else if(args[i].equals("-ms")){
                this.minimum_samples = Integer.parseInt(args[i+1]);
                if(this.minimum_samples < 0){this.minimum_samples = 0;}
                i++;
            }else if(args[i].equals("-msp")){
                minimum_samples_percentage = Double.parseDouble(args[i+1]);
                if (minimum_samples_percentage < 0){minimum_samples_percentage = 0;}
                if (minimum_samples_percentage > 1){minimum_samples_percentage = 1;}
                i++;
            }else if(args[i].equals("-w")){
                this.weighted_counts = true;
                if (args[i+1].toLowerCase().equals("false")){
                    this.weighted_counts = false;
                }
                i++;
            }else if(args[i].equals("-r")){
                this.round_to_zero = true;
                if (args[i+1].toLowerCase().equals("false")){
                    this.round_to_zero = false;
                }
                i++;
            }else if(args[i].equals("-eo")){
                this.enzymeOverlap = Integer.parseInt(args[i+1]);
                i++;
            }else if(args[i].equals("-prefix")){
                this.filePrefix = args[i+1] + ".";
                i++;
            }
        }
        
        for (String bedFile : this.bedFiles){
            System.out.println("Now reading: " + bedFile);
            this.addBed(bedFile);
        }
        this.samples = new String[this.bedFiles.size()];
        int count_of_bedFiles = 0;
        //create an array of all samplenames (for order)
        for (String file : this.bedFiles){
            String[] bedfilesplit = file.split(System.getProperty("file.separator"));
            String sampleName = bedfilesplit[bedfilesplit.length - 1];
            this.samples[count_of_bedFiles] = sampleName;
            count_of_bedFiles++;
        }
        if(this.minimum_samples > count_of_bedFiles){
            this.minimum_samples = count_of_bedFiles;
        }
        if(this.minimum_samples < (count_of_bedFiles * minimum_samples_percentage)){
            this.minimum_samples = (int) (minimum_samples_percentage * count_of_bedFiles);
        }
        System.out.println("Start new approache test");
        System.out.println("Output:\t" + this.outputDir);
        System.out.println("Minimum Occurances:\t" + this.minimum_occurance);
        System.out.println("Minimum Samples:\t" + this.minimum_samples);
        System.out.println("Enzyme Overlap:\t" + this.enzymeOverlap);
        this.unitList = new ArrayList<>();
        this.unitGeneration();
        this.showFragmentAndUnitStats();
        this.printLengthMap();
        this.printUnits(this.unitList, "");
        //this.printSimpleUnits();
        ArrayList<Unit> simpleunits = this.filterOnSimpleUnits();
        this.printUnits(simpleunits, ".simple");
        ArrayList<Unit> doubleunits = this.filterOnSimpleAndDoubleUnits();
        this.printUnits(doubleunits, ".double");
        //this.normalization();
    }
    
    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        long factor = (long) Math.pow(10, places);
        value = value * factor;
        long tmp = Math.round(value);
        return (double) tmp / factor;
    }
    
    private String outputDir;
    private ArrayList<String> bedFiles;
    private HashMap<String, HashMap<String, HashMap<String, Integer>>> chromosome_location_sample_countMap;
    private HashMap<String, HashMap<String, Integer>> chromosome_sample_countMap;
    private double minimum_occurance;
    private int minimum_samples;
    private boolean weighted_counts;
    private boolean round_to_zero;
    private String[] samples;
    private int enzymeOverlap;
    private String filePrefix;
    private ArrayList<Unit> unitList;
     
    
    private void addBed(String bedFile){
        try {
            String[] bedfilesplit = bedFile.split("/");
            String sampleName = bedfilesplit[bedfilesplit.length - 1];
            BedReader bedReader = new BedReader(bedFile);
            BedRead bedRead;
            while ((bedRead = bedReader.getNextRead()) != null){
                if (this.chromosome_sample_countMap.containsKey(bedRead.getChromosome())){
                    HashMap<String, Integer> sample_countMap = this.chromosome_sample_countMap.get(bedRead.getChromosome());
                    if (sample_countMap.containsKey(sampleName)){
                        sample_countMap.put(sampleName, sample_countMap.get(sampleName) + 1);
                    }else{
                        sample_countMap.put(sampleName, 1);
                    }
                }else{
                    HashMap<String, Integer> sample_countMap = new HashMap<>();
                    sample_countMap.put(sampleName, 1);
                    this.chromosome_sample_countMap.put(bedRead.getChromosome(), sample_countMap);
                }
                
                String start_end_location = bedRead.getBeginPos() + "_" + bedRead.getEndPos();
                if (! this.chromosome_location_sample_countMap.containsKey(bedRead.getChromosome())){
                    this.chromosome_location_sample_countMap.put(bedRead.getChromosome(), new HashMap<String, HashMap<String, Integer>>());                    
                }
                HashMap<String, HashMap<String, Integer>> location_sample_countMap = this.chromosome_location_sample_countMap.get(bedRead.getChromosome());
                if (! location_sample_countMap.containsKey(start_end_location)){
                    location_sample_countMap.put(start_end_location, new HashMap<String, Integer>());
                }
                HashMap<String, Integer> sample_countMap = location_sample_countMap.get(start_end_location);
                if (! sample_countMap.containsKey(sampleName)){
                            sample_countMap.put(sampleName, 0);
                }
                sample_countMap.put(sampleName, sample_countMap.get(sampleName) + 1);
                
                
                
//                if (this.chromosome_location_sample_countMap.containsKey(bedRead.getChromosome())){
//                    HashMap<String, HashMap<String, Integer>> location_sample_countMap = this.chromosome_location_sample_countMap.get(bedRead.getChromosome());
//                    if (location_sample_countMap.containsKey(start_end_location)){
//                        HashMap<String, Integer> sample_countMap = location_sample_countMap.get(start_end_location);
//                        if (sample_countMap.containsKey(sampleName)){
//                            sample_countMap.put(sampleName, sample_countMap.get(sampleName) + 1);
//                        }else{
//                            sample_countMap.put(sampleName, 1);
//                        }
//                    }else{
//                        HashMap<String, Integer> sample_countMap = new HashMap<>();
//                        sample_countMap.put(sampleName, 1);
//                        location_sample_countMap.put(start_end_location, sample_countMap);
//                    }
//                }else{
//                    HashMap<String, HashMap<String, Integer>> location_sample_countMap = new HashMap<>();
//                    HashMap<String, Integer> sample_countMap = new HashMap<>();
//                    sample_countMap.put(sampleName, 1);
//                    location_sample_countMap.put(start_end_location, sample_countMap);
//                    this.chromosome_location_sample_countMap.put(bedRead.getChromosome(), location_sample_countMap);
//                }
            }
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Uniter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void unitGeneration(){
    
        ArrayList<Fragment> fragmentList = new ArrayList<>();
        
        int removedFragments = 0;
        for (String chromosome : this.chromosome_location_sample_countMap.keySet()){
            for (String location : this.chromosome_location_sample_countMap.get(chromosome).keySet()){
                int begin = Integer.parseInt(location.split("_")[0]);
                int end = Integer.parseInt(location.split("_")[1]);
                Fragment f = new Fragment(this.samples, chromosome, begin, end, this.enzymeOverlap);
                for (String sample : this.samples){
                    if (this.chromosome_location_sample_countMap.get(chromosome).get(location).containsKey(sample)){
                        f.addSampleCount(sample, this.chromosome_location_sample_countMap.get(chromosome).get(location).get(sample));
                    }else{
                        f.addSampleCount(sample, 0);
                    }
                }
                if (f.filter(this.minimum_samples, this.minimum_occurance)){
                    fragmentList.add(f);
                }else{
                    //do not add
                    removedFragments++;
                }
            }
        }
        System.out.println("\t\tNumber of Fragments removed:\t" + removedFragments + "");
        
        Collections.sort(fragmentList, new FragmentComparator());
        this.unitList = new ArrayList<>();
        
        Unit unit = null;
        for (Fragment f : fragmentList){
            if (unit == null){
                unit = new Unit(f);
                this.unitList.add(unit);
            }else{
                if (unit.canAddFragment(f)){
                    unit.addFragment(f);
                }else{
                    unit = new Unit(f);
                    this.unitList.add(unit);
                }
            }
        }
        
        
        System.out.println("First  \tUnits: " + this.unitList.size());
        
        System.out.println("Number of Fragments: " + fragmentList.size());
        System.out.println("Number of Units:     " + this.unitList.size());
        Unit testUnit = null;
        int numberOfSimpleUnits = 0;
        int unit2 = 0;
        int unit3 = 0;
        int unit4 = 0;
        int unitM = 0;
        for (Unit u : this.unitList){
            if (u.isSimpleUnit()){
                numberOfSimpleUnits++;
            }else{
                if (u.getNumberOfFragments() == 2){
                    unit2++;
                }
                if (u.getNumberOfFragments() == 3){
                    unit3++;
                }
                if (u.getNumberOfFragments() == 4){
                    unit4++;
                }
                if (u.getNumberOfFragments() > 4){
                    unitM++;
                    testUnit = u;
                }
            }
        }
        System.out.println("Number of Simple Units: " + numberOfSimpleUnits);
        System.out.println("Number of Units Sites 2: " + unit2);
        System.out.println("Number of Units Sites 3: " + unit3);
        System.out.println("Number of Units Sites 4: " + unit4);
        System.out.println("Number of Units Sites M: " + unitM);
    }
    
    private void normalization(){
        ArrayList<Fragment> fragmentList = new ArrayList<>();
        
        for(Unit u : this.unitList){
            fragmentList.addAll(u.getFragments());
        }
        //create a map with the total number of reads per sample
        HashMap<String, Integer> sample_totalMap = new HashMap<>();

        for (String sample : this.samples){
            sample_totalMap.put(sample, 0);
        }
        for (Unit u : this.unitList){
            for (Fragment f : u.getFragments()){
                for(String sample : this.samples){
                    sample_totalMap.put(sample, sample_totalMap.get(sample) + f.getCounts(sample));
                }
            }

        }
            
        //find the highest number of total reads
        int highest_read_number = 0;
        int lowest_read_number = -1;
        for (String s : this.samples){
            int read_number = sample_totalMap.get(s);
            if (highest_read_number <= read_number){
                highest_read_number = read_number;
            }
            if (lowest_read_number == -1 || lowest_read_number > read_number){
                lowest_read_number = read_number;
            }
        }
            
        //create a map with the factor for normalisation
        HashMap<String, Double> sample_factor_high_map = new HashMap<>();
        HashMap<String, Double> sample_factor_low_map = new HashMap<>();
        double highest_factor_normal = 1;
        double lowest_factor_normal = 1;
        for (String s : this.samples){
            sample_factor_high_map.put(s, ((highest_read_number * 1.0) / (sample_totalMap.get(s) * 1.0)));
            sample_factor_low_map.put(s, ((lowest_read_number * 1.0) / (sample_totalMap.get(s) * 1.0)));
            if (highest_factor_normal < ((highest_read_number * 1.0) / (sample_totalMap.get(s) * 1.0))){
                highest_factor_normal = ((highest_read_number * 1.0) / (sample_totalMap.get(s) * 1.0));
            }
            if (lowest_factor_normal > ((lowest_read_number * 1.0) / (sample_totalMap.get(s) * 1.0))){
                lowest_factor_normal = ((lowest_read_number * 1.0) / (sample_totalMap.get(s) * 1.0));
            }
        }
            
        System.out.println("Quartile size factor");
        HashMap<String, Double> sample_factor_fragNorm_map = new HashMap<>();
        HashMap<String, Double> sample_factor_fragQuart_map = new HashMap<>();
        for (String s : this.samples){
            ArrayList<Double> fragNorm_list = new ArrayList<>();
            for (Fragment f : fragmentList){
                if (f.getSizeFactor(s) != 0)
                    fragNorm_list.add(f.getSizeFactor(s, sample_factor_high_map));
            }
            Collections.sort(fragNorm_list);
            if (fragNorm_list.size() % 2 == 0){
                sample_factor_fragNorm_map.put(s, (
                        fragNorm_list.get(fragNorm_list.size() /2) + fragNorm_list.get(fragNorm_list.size() /2 -1)/2.0));
            }else{
                sample_factor_fragNorm_map.put(s, fragNorm_list.get(fragNorm_list.size() /2));
            }
            double quartalcount = 0;
            long start = Math.round(fragNorm_list.size()*(0.25));
            long end = Math.round(fragNorm_list.size()*(0.75));
            for (int i = (int) start; i<end; i++){
                quartalcount+=fragNorm_list.get(i);
            }
            System.out.println(""+s+"\t"+quartalcount/((end-start)*1.0));
            sample_factor_fragQuart_map.put(s, quartalcount/((end-start)*1.0));
        }
        HashMap<String, Double> sample_factor_unitNorm_map = new HashMap<>();
        for (String s : this.samples){
            ArrayList<Double> unitNorm_list = new ArrayList<>();
            for (Unit u : this.unitList){
                if (u.getSizeFactor(s) != 0)
                    unitNorm_list.add(u.getSizeFactor(s));
            }
            Collections.sort(unitNorm_list);
            if (unitNorm_list.size() % 2 == 0){
                sample_factor_unitNorm_map.put(s, (
                        unitNorm_list.get(unitNorm_list.size() /2) + unitNorm_list.get(unitNorm_list.size() /2 -1)/2.0));
            }else{
                sample_factor_unitNorm_map.put(s, unitNorm_list.get(unitNorm_list.size() /2));
            }
        }
        
        System.out.println("Quartile library normalization");
        HashMap<String, Double> sample_factor_quart_lib_norm_map = new HashMap<>();
        Double highest = 0.0;
        for (String s : this.samples){
            ArrayList<Integer> count_list = new ArrayList();
           for (Unit u : this.unitList){
               count_list.add(u.getUnitDepth(s));
           } 
           Collections.sort(count_list);
            double quartalcount = 0;
            long start = Math.round(count_list.size()*(0.25));
            long end = Math.round(count_list.size()*(0.75));
            for (int i = (int) start; i<end; i++){
                quartalcount+=count_list.get(i);
            }
            System.out.println(""+s+"\t"+quartalcount/((end-start)*1.0));
            sample_factor_quart_lib_norm_map.put(s, quartalcount/((end-start)*1.0));
            if(highest < quartalcount/((end-start)*1.0)) highest = quartalcount/((end-start)*1.0);
        }
        for (String s : this.samples){
            System.out.println(""+s+"\t"+sample_factor_quart_lib_norm_map.get(s)/highest);
        }
        
        
            
        System.out.println("Sample\tHigh\tLow\tFragNorm\tUnitNorm\tFragQuart");
        for (String s : this.samples){
            System.out.println(s + "\t" + sample_factor_high_map.get(s) + "\t" + sample_factor_low_map.get(s) + "\t" + sample_factor_fragNorm_map.get(s)
                    + "\t" + sample_factor_unitNorm_map.get(s) + "\t" + sample_factor_fragQuart_map.get(s));
        }
            
        
            
            
    }
    
    private ArrayList<Unit> filterOnSimpleUnits(){
        ArrayList<Unit> simpleUnitList = new ArrayList<>();
            for (Unit u : this.unitList){
                if(u.getNumberOfFragments()==1){
                    simpleUnitList.add(u);
                }
            }
        Collections.sort(simpleUnitList, new UnitComparator());
        return simpleUnitList;
    }
    
    
    private ArrayList<Unit> filterOnSimpleAndDoubleUnits(){
        ArrayList<Unit> doubleUnitList = new ArrayList<>();
            for (Unit u : this.unitList){
                if(u.getNumberOfFragments()<=2){
                    doubleUnitList.add(u);
                }
            }
        Collections.sort(doubleUnitList, new UnitComparator());
        return doubleUnitList;
    }
    
    
    private void showFragmentAndUnitStats(){
        
        ArrayList<Fragment> fragmentList = new ArrayList<>();
        
        for(Unit u : this.unitList){
            fragmentList.addAll(u.getFragments());
        }
        
        HashMap<String, Integer> frag_sample_count = new HashMap<>();
        HashMap<String, Integer> unit_sample_count = new HashMap<>();
        HashMap<String, Integer> frag_sample_uniq_count = new HashMap<>();
        HashMap<String, Integer> unit_sample_uniq_count = new HashMap<>();
        for (String s : this.samples){
            frag_sample_count.put(s, 0);
            unit_sample_count.put(s, 0);
            frag_sample_uniq_count.put(s, 0);
            unit_sample_uniq_count.put(s, 0);
        }
        for(Fragment f : fragmentList){
            for (String s : this.samples){
                if (f.getCounts(s) > 0){
                    frag_sample_count.put(s, frag_sample_count.get(s) + 1);
                    if(f.getOccuranceSamples() == 1){
                        frag_sample_uniq_count.put(s, frag_sample_uniq_count.get(s) + 1);
                    }
                }
            }
        }
        for (Unit u : this.unitList){
            for (String s : this.samples){
                if(u.getUnitDepth(s) > 0){
                    unit_sample_count.put(s, unit_sample_count.get(s) + 1);
                    HashSet<String> sampleList = new HashSet();
                    for (Fragment f : u.getFragments()){
                        for(String samples : f.getSamples()){
                            sampleList.add(samples);
                        }
                    }
                    if (u.getNumberOfSamples() == 1){
                        unit_sample_uniq_count.put(s, unit_sample_uniq_count.get(s) + 1);
                    }
                }
            }
        }
        System.out.println("sample\tfragment count\tunit count\tuniq fragment count\tuniq unit count");
        for (String s : samples){
            System.out.println("" + s + "\t" + frag_sample_count.get(s) + "\t" + unit_sample_count.get(s)
                                    + "\t" + frag_sample_uniq_count.get(s) + "\t" + unit_sample_uniq_count.get(s));
        }

    }
    
    private void printUnits(ArrayList<Unit> units, String fix){
        Collections.sort(units, new UnitComparator());
            
        try{
            MyWriter locationWriter = new MyWriter(this.outputDir, this.filePrefix + "location.variation" + fix + ".mo" + this.minimum_occurance + ".ms" 
                    + this.minimum_samples + ".eo" + this.enzymeOverlap + ".vgbsf");
            MyWriter bedWriter = new MyWriter(this.outputDir, this.filePrefix + "interesting.locations" + fix + ".mo" + this.minimum_occurance + ".ms" 
                    + this.minimum_samples + ".eo" + this.enzymeOverlap + ".bed");
            String line = "Chr\tBegin\tEnd\tLociVariants";
            String bedLine;
            for (String sample : this.samples){
                line += "\t" + sample;
            }
            locationWriter.addToFile(line + "\n");
            for (Unit u : units){
                line = u.getChromosome() + "\t" + u.getBegin() + "\t" + u.getEnd() + "\t";
                bedLine = u.getChromosome() + "\t" + u.getBegin() + "\t" + u.getEnd();
                int count = 0;
                for (Fragment f : u.getFragments()){
                    if (count != 0){
                        line += ",";
                    }
                    count++;
                    line += f.getChromosome() + ":" + f.getBegin() + "-" + f.getEnd();
                }
                for (String sample : this.samples){
                    line += "\t";
                    int count2 = 0;
                    for (Fragment f : u.getFragments()){
                        if (count2 != 0){
                            line += ",";
                        }
                        count2++;
                        line += f.getCounts(sample);
                    }
                }
                locationWriter.addToFile(line + "\n");
                bedWriter.addToFile(bedLine + "\n");
            }
            locationWriter.closeWriter();
            bedWriter.closeWriter();
        } catch (IOException ex) {
            Logger.getLogger(Uniter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void printSimpleUnits(){
        Collections.sort(this.unitList, new UnitComparator());
            
        try{
            MyWriter locationWriter = new MyWriter(this.outputDir, this.filePrefix + "location.variation.simple.mo" + this.minimum_occurance + ".ms" 
                    + this.minimum_samples + ".eo" + this.enzymeOverlap + ".vgbsf");
            MyWriter bedWriter = new MyWriter(this.outputDir, this.filePrefix + "interesting.locations.simple.mo" + this.minimum_occurance + ".ms" 
                    + this.minimum_samples + ".eo" + this.enzymeOverlap + ".bed");
            String line = "Chr\tBegin\tEnd\tLociVariants";
            String bedLine;
            for (String sample : this.samples){
                line += "\t" + sample;
            }
            locationWriter.addToFile(line + "\n");
            for (Unit u : this.unitList){
                if(u.getNumberOfFragments()==1){
                    line = u.getChromosome() + "\t" + u.getBegin() + "\t" + u.getEnd() + "\t";
                    bedLine = u.getChromosome() + "\t" + u.getBegin() + "\t" + u.getEnd();
                    int count = 0;
                    for (Fragment f : u.getFragments()){
                        if (count != 0){
                            line += ",";
                        }
                        count++;
                        line += f.getChromosome() + ":" + f.getBegin() + "-" + f.getEnd();
                    }
                    for (String sample : this.samples){
                        line += "\t";
                        int count2 = 0;
                        for (Fragment f : u.getFragments()){
                            if (count2 != 0){
                                line += ",";
                            }
                            count2++;
                            line += f.getCounts(sample);
                        }
                    }
                    locationWriter.addToFile(line + "\n");
                    bedWriter.addToFile(bedLine + "\n");
                }
            }
            locationWriter.closeWriter();
            bedWriter.closeWriter();
        } catch (IOException ex) {
            Logger.getLogger(Uniter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void printLengthMap(){
        HashMap<Integer, HashMap<String, Integer>> length_sample_count_map = new HashMap<>();

        for (Unit u : unitList){
            for (Fragment f : u.getFragments()){
                if (! length_sample_count_map.containsKey(f.getLength())){
                    length_sample_count_map.put(f.getLength(), new HashMap<String, Integer>());
                }
                for (String s : f.getSamples()){
                    if (length_sample_count_map.get(f.getLength()).containsKey(s)){
                        length_sample_count_map.get(f.getLength()).put(s, length_sample_count_map.get(f.getLength()).get(s) + f.getCounts(s));
                    }else{
                        length_sample_count_map.get(f.getLength()).put(s, f.getCounts(s));
                    }
                }
            }
        }
        try{
            MyWriter lengthWriter = new MyWriter(this.outputDir, this.filePrefix + "length_counts.txt");
            int maxLength = 0;
            for (int l : length_sample_count_map.keySet()){
                if (l > maxLength){
                    maxLength = l;
                }
            }
            String line = "Length";
            for (String s : this.samples){
                line += "\t" + s;
            }
            lengthWriter.addToFile(line + "\n");
            for (int i =0; i <= maxLength; i++){
                line = "" + i;
                for (String s : this.samples){
                    int sample_count = 0;
                    if (length_sample_count_map.containsKey(i)){
                        if (length_sample_count_map.get(i).containsKey(s)){
                            sample_count = length_sample_count_map.get(i).get(s);
                        }
                    }
                    line += "\t" + sample_count;
                }
                lengthWriter.addToFile(line + "\n");
            }
            lengthWriter.closeWriter();
        } catch (IOException ex) {
            Logger.getLogger(Uniter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    
    
    
    
}
