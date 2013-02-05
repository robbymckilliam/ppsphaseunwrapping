/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import Jama.Matrix;
import pubsim.lattices.VnmStarGlued;
import pubsim.lattices.decoder.SphereDecoderSchnorrEuchner;

/**
 * Runs the VnmStarSampledEfficient nearest lattice point algorithm to
 * estimate m polynomial phase signal.
 * @author Robby McKilliam
 */
public class SphereDecoderEstimator extends BabaiEstimator {
    
    /** 
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public SphereDecoderEstimator(int m, int n) {
        lattice = new VnmStarGlued(m, n - m - 1);
        this.m = m;
        npalgorithm = new SphereDecoderSchnorrEuchner(lattice);
        ambiguityRemover = new AmbiguityRemover(m);

        ya = new double[n];
        p = new double[m+1];
        this.n = n;

        M = lattice.getMMatrix();
        Matrix Mt = M.transpose();
        K = Mt.times(M).inverse().times(Mt);
    }
   

}