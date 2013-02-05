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
        lattice = new VnmStarGlued(m, n - m - 1);
        this.m = m;
        npalgorithm = new Mbest(lattice,M,lr);
        ambiguityRemover = new AmbiguityRemover(m);

        ya = new double[n];
        p = new double[m+1];
        this.n = n;

        this.M = lattice.getMMatrix();
        Matrix Mt = this.M.transpose();
        K = Mt.times(this.M).inverse().times(Mt);
    }
}