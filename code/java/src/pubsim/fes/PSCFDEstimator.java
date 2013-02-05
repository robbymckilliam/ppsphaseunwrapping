/*
 * PSCFDEstimator.java
 *
 * Created on 11 September 2007, 08:51
 */

package pubsim.fes;

/**
 * The Parabolic Smoothed Central Finite Difference Estimator
 * @author Robby McKilliam
 */
public class PSCFDEstimator implements FrequencyEstimator {
       
    protected final int N;

    public PSCFDEstimator(int N){
        this.N = N;
    }

    /** Run the estimator on received data, @param y */
    public double estimateFreq(double[] real, double[] imag){
        
        double ks = 0.0;
        double sumre = 0.0;
        double sumim = 0.0;
        for(int i=0; i<=N-2; i++){
            double re = real[i]*real[i+1] + imag[i]*imag[i+1];
            double im = real[i]*imag[i+1] - real[i+1]*imag[i];
            double mag = Math.sqrt((real[i]*real[i] + imag[i]*imag[i])*(real[i+1]*real[i+1] + imag[i+1]*imag[i+1]));
            sumre += 6.0*(i+1)*(N-1-i)/(N*(N*N-1))*re/mag;
            sumim += 6.0*(i+1)*(N-1-i)/(N*(N*N-1))*im/mag;
        }
        return Math.atan2(sumim,sumre)/(2*Math.PI);
    }
    
}
