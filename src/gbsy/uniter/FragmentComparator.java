/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbsy.uniter;

import java.util.Comparator;

/**
 *
 * @author koen
 */
public class FragmentComparator implements Comparator<Fragment>{

    @Override
    public int compare(Fragment o1, Fragment o2) {
        return o1.compareTo(o2);
    }

}
    
