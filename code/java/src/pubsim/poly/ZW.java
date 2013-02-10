/*
 */
package pubsim.poly;

/**
 * An implementation of Zhou and Wang's Euclidean algorithm approach for increasing
 * the range of parameters for which the HAF/DPT applies.
 * @author Robby McKilliam
 */
public class ZW extends AbstractPolynomialPhaseEstimator {
    
    final protected HAF dpt1, dpt2;
    
    public ZW(int m, int n, int tau1, int tau2){
        super(m);
        dpt1 = new HAF(m, n,tau1);
        dpt2 = new HAF(m, n,tau2);
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}
