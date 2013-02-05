/*
 */

package pubsim.poly;

import Jama.Matrix;
import pubsim.lattices.VnmStarGlued;
import pubsim.lattices.decoder.Mbest;
import pubsim.lattices.reduction.LLL;
import pubsim.lattices.reduction.LatticeReduction;

/**
 *
 * @author Robby McKilliam
 */
public class MbestEstimator extends BabaiEstimator {

    /**
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public MbestEstimator(int m, int n, int M) {
        this(m,n,M,new LLL());
    }
    
     /**
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public MbestEstimator(int m, int n, int M, LatticeReduction lr) {
        super(m,n);
        npalgorithm = new Mbest(lattice,M,lr);
    }
}