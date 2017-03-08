/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy;

import gbsy.uniter.Uniter;

/**
 *
 * @author koen
 */
public class GBSY {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        String input = "uniter "
//                + "-b /home/koen/Downloads/bed_test/lib1_DLT_1.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_DLT_3.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_GRSP_2.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_PK_1.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_STY_1.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_STY_3.bed.LG10.bed "
//                + "-b /home/koen/Downloads/bed_test/lib1_DLT_2.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_GRSP_1.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_GRSP_3.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_PK_3.bed.LG10.bed  "
//                + "-b /home/koen/Downloads/bed_test/lib1_STY_2.bed.LG10.bed "
//                + "-o /home/koen/Downloads/bed_test/result "
//                + "-mo 5 "
//                + "-msp 0.5 "
//                + "-prefix gbsy";
//        args = input.split(" ");
        // TODO code application logic here
        if(args.length==0){
            System.out.println("Subprograms are:");
            System.out.println("\tuniter");
        }else if(args[0].equals("uniter")){
            String[] newargs = new String[args.length-1];
            for(int i=1; i < args.length; i++){
                newargs[i-1] = args[i];
            }
            Uniter bedToCount = new Uniter(newargs);
        }
    }
    
}
