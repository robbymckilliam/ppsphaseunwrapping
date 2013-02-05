/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

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
        super(m,n);
        npalgorithm = new SphereDecoderSchnorrEuchner(lattice);
    }
   

}