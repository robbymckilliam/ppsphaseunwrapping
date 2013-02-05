/*
 * GlueAnstarEstimator.java
 *
 * Created on 12 August 2007, 12:28
 */

package pubsim.fes;

import pubsim.lattices.Vn1Star.Vn1StarGlued;

/**
 * Frequency estimator that uses Vn1 glue vector algorithm to solve the nearest
 * point problem for the frequency estimation lattice Vn1.  O(N^4log(N)).
 * @author Robby McKilliam
 */
public class GlueAnstarEstimator 
        extends LatticeEstimator
        implements FrequencyEstimator {
    
   public GlueAnstarEstimator(int N){
       super(N);
       lattice = new Vn1StarGlued(N-2);
    }
    
}
