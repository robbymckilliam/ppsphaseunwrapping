/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly.bounds;

import Jama.Matrix;
import static pubsim.Util.powerSum;

/**
 * Interface for computing estimator bounds for polynomial phase parameters
 * @author Robby McKilliam
 */
public abstract class BoundCalculator {

    /**
     * Matrix describing covariances.  The
     * appears in every bound that I know of.
     */
    protected Matrix C;

    /**
     * Signal length.
     */
    int N;

    /**
     * polynomial order.
     */
    int m;

    protected BoundCalculator() {}

    /** 
     * @param N signal length
     * @param m polynomial order
     */
    public BoundCalculator(int N, int m){

        this.N = N;
        this.m = m;

        //if N is pretty big compared to m then the following works well
//        if( N > 20*m ){
            Matrix M = new Matrix(m+1, m+1);
            for(int i = 0; i <= m; i++){
                for(int j = 0; j <= m; j++){
                    M.set(i,j, 1.0/(i + j + 1));
                }
            }
            C = M.inverse();
//        }
        //else{
//            Matrix M = new Matrix(m+1, m+1);
//            for(int i = 0; i <= m; i++){
//                for(int j = 0; j <= m; j++){
//                    M.set(i,j, Math.pow(N, -(j+i))*powerSum(N,i+j));
//                }
//            }
//            C = M.inverse();
        //}

    }


}
