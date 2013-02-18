/*
 */

package pubsim.poly;

import pubsim.lattices.reduction.LLL;
import pubsim.lattices.reduction.LatticeReduction;

/**
 * @author Robby McKilliam
 */
public class Mbest extends Babai {

    /**
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public Mbest(int m, int n, int M) {
        this(m,n,M,new LLL());
    }
    
     /**
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public Mbest(int m, int n, int M, LatticeReduction lr) {
        super(m,n);
        npalgorithm = new pubsim.lattices.decoder.Mbest(lattice,M,lr);
    }
}