/*
 */
package pubsim.poly;

/**
 * Sets up polynomial phase estimators in a standard way
 * @author Robby McKilliam
 */
abstract public class AbstractPolynomialPhaseEstimator implements PolynomialPhaseEstimatorInterface {

    final protected AmbiguityRemoverRectangular ambiguityRemover;
    final int m;

    private AbstractPolynomialPhaseEstimator() {
        throw new UnsupportedOperationException("Should not be called!");
    }

    public AbstractPolynomialPhaseEstimator(int m) {
        ambiguityRemover = new AmbiguityRemoverRectangular(m);
        this.m = m;
    }

    @Override
    final public double[] error(double[] real, double[] imag, double[] truth) {
        double[] est = estimate(real, imag);
        double[] err = new double[est.length];

        for (int i = 0; i < err.length; i++) {
            err[i] = est[i] - truth[i];
        }
        err = ambiguityRemover.disambiguate(err);
        //for (int i = 0; i < err.length; i++) {
        //    err[i] *= err[i];
        //}
        return err;
    }
    
    @Override
    final public int getOrder() {
        return m;
    }
    
}
