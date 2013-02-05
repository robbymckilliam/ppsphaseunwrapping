/*
 * KaysEstimator.java
 *
 * Created on 12 August 2007, 20:06
 */

package pubsim.fes;

/**
 * Implementation of Kay's phase unwrapping frequency estimator
 * @author Robby McKilliam
 */
public class KaysEstimator implements FrequencyEstimator{
    
    protected final int N;

    public KaysEstimator(int N){
        this.N = N;
    }
    
    /** Run the estimator on recieved data, @param y */
    public double estimateFreq(double[] real, double[] imag){  
        double ks = 0.0;
        for(int i=0; i<=N-2; i++){
            double re = real[i]*real[i+1] + imag[i]*imag[i+1];
            double im = real[i]*imag[i+1] - real[i+1]*imag[i];
            ks += (i + 1)*(N - 1 - i)*Math.atan2(im,re);
        }
        return ks*6.0/(2*Math.PI*N*(N*N - 1));
    }
    
}
