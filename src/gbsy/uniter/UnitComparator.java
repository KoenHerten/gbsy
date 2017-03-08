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
public class UnitComparator implements Comparator<Unit>{

    @Override
    public int compare(Unit o1, Unit o2) {
        return o1.compareTo(o2);
    }

}
