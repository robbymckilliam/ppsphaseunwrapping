/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import java.io.Serializable;

/**
 * 
 * @author Robby McKilliam
 */
public interface PolynomialPhaseEstimatorInterface extends Serializable{

    /** 
     * Return the order of this estimator.
     * This is the number of parameters.
     */
    public int getOrder();
    
    /** 
     * Run the estimator on received real and imaginary signal.
     * Returns an array of polynomial parameters.
     */
    public double[] estimate(double[] real, double[] imag);

    /**
     * Run the estimator and return the square error
     * wrapped modulo the the ambiguity region between
     * the estimate and the truth.
     */
    public double[] error(double[] real, double[] imag, double[] truth);
    
}
