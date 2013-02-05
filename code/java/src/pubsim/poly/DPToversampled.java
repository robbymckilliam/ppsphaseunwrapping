/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

/**
 * The is the DPT estimator that has been oversampled so that the
 * volume of parameters is equal to that of the identifiable region.
 * @author Robby McKilliam
 */
public class DPToversampled extends DPTEstimator {

    @Override
    public double[] error(double[] real, double[] imag, double[] truth) {

        //this computes the 'fair' error after oversampling
        double dptos = Math.pow( (n*1.0)/m , (m-1.0)/(m+1.0));

        double[] est = estimate(real, imag);
        double[] err = new double[est.length];

        for (int i = 0; i < err.length; i++) {
            double oscale = Math.pow(dptos,i);
            double oserr = oscale*(est[i] - truth[i]);
            double osr = 1.0;
            if( i > 1 ) osr = oscale*Math.pow(m/(n*1.0), i-1)/pubsim.Util.factorial(i);

            err[i] = osr*pubsim.Util.fracpart(oserr/osr);
        }

        //err = ambiguityRemover.disambiguate(err);
        for (int i = 0; i < err.length; i++) {
            err[i] *= err[i];
        }
        return err;
    }

}
