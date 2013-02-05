/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import pubsim.Util;

/**
 * Implementation of the High-order phase function estimator.
 * UNDER CONSTRUCTION!
 * @author robertm
 */
public class HPFEstimator implements PolynomialPhaseEstimator {

    public void setSize(int n) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getOrder() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double[] estimate(double[] real, double[] imag) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public double[] error(double[] real, double[] imag, double[] truth) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Returns the volume of the functional regions of
     * the HPF estimator.  This is just for comparison
     * with the identifiable region.
     * @return volume of functional region
     */
    public static double volumeOfFunctionalRegion(int n, int m){
        double prod = 1.0;
        for(int k = 2; k <= m; k++){
            prod *= 1.0/( Util.factorial(k) * Math.pow(n/2.0, k-1) );
        }
        return prod;
    }



}
