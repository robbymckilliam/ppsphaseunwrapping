/*
 */
package pubsim.poly;

/**
 * Implements O'Shea's cubic phase function estimator
 * Peter O'Shea "A Fast Algorithm for Estimating the Parameters of a Quadratic FM Signal"
 * IEEE Trans. Signal Proc. Vol 53, Feb 2004.
 * @author Robby McKilliam
 */
public class CubicPhaseFunction extends AbstractPolynomialPhaseEstimator {

    public CubicPhaseFunction(int m, int n){
        super(m);
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        throw new UnsupportedOperationException("Not supported yet.");
    }    

        
}
