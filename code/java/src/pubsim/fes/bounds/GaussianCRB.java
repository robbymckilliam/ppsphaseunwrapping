/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.fes.bounds;

/**
 * Standard Cramer-Rao bound for frequency estimation.
 * @author Robby McKilliam
 */
public class GaussianCRB {
    
    public static double getBound(double var, double n) {
        return getBound(var,1.0, n);
    }

    public static double getBound(double var, double amplitude, double n) {
        return 3.0*var/(Math.PI*Math.PI*amplitude*amplitude*n*(n*n-1));
    }
    
}
