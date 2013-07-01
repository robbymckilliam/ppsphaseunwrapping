package pubsim.poly;

import pubsim.lattices.decoder.SphereDecoderSchnorrEuchner;

/**
 * Runs the VnmStarSampledEfficient nearest lattice point algorithm to
 * estimate m polynomial phase signal.
 * @author Robby McKilliam
 */
public class SphereDecoder extends Babai {
    
    /** 
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public SphereDecoder(int m, int n) {
        super(m,n);
        npalgorithm = new SphereDecoderSchnorrEuchner(lattice);
    }
   

}