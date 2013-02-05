/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.fes;

import pubsim.lattices.Vn1Star.Vn1StarZnLLS;

/**
 * O(N^3log(N)) version of the LSPU estimator.
 * @author Robby McKilliam
 */
public class ZnLLS extends LatticeEstimator implements FrequencyEstimator{

    public ZnLLS(int N){
        super(N);
        lattice = new Vn1StarZnLLS(N-2);
    }

}
