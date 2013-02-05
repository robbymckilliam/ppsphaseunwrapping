/*
 * PeriodogramEstimator.java
 *
 * Created on 12 August 2007, 12:12
 */

package pubsim.fes;

/**
 * Periodogram Estimator for frequency. This uses a slow O(N^2) Fourier
 * transform.  Use PeriodogramFFTEstimator if you want the fast one.
 * @author Robby McKilliam
 */
public class PeriodogramEstimator implements FrequencyEstimator{
    
    protected final int N;
    protected final int num_samples;
    
    /**Max number of iterations for the Newton step */
    static final int MAX_ITER = 15;
    
    /**Step variable for the Newton step */
    static final double EPSILON = 1e-10;
    
    public PeriodogramEstimator(int N){
        num_samples = 100;
        this.N = N;
    }
    
    /** Constructor that sets the number of samples to be taken of
     * the periodogram.
     */
    public PeriodogramEstimator(int N, int samples){
        num_samples = samples;
        this.N = N;
    }
    
    public static double calculatePeriodogram(double[] real, double[] imag, double f) {
        double sumur = 0, sumui = 0;
        for (int i = 0; i < real.length; i++) {
            double cosf =  Math.cos(-2 * Math.PI * f * i);
            double sinf =  Math.sin(-2 * Math.PI * f * i);
            sumur += real[i]*cosf - imag[i]*sinf;
            sumui += imag[i]*cosf + real[i]*sinf;
        }
        return sumur * sumur + sumui * sumui;
    }
    
    /** Run the periodogram estimator on received data, @param y */
    @Override
    public double estimateFreq(double[] real, double[] imag){

        // Coarse search
        double maxp = 0;
        double fhat= 0.0;
        double fstep = 1.0/num_samples;
	for (double f = -0.5; f <= 0.5; f += fstep) {
	    double p = calculatePeriodogram(real, imag, f);
	    if (p > maxp) {
		maxp = p;
		fhat = f;
	    }
	}
        
        
        //Newton Raphson
        int numIter = 0;
        double f = fhat, lastf = f - 2 * EPSILON, lastp = 0;
        while (Math.abs(f - lastf) > EPSILON && numIter <= MAX_ITER) {
            double p = 0, pd = 0, pdd = 0;
	    double sumur = 0, sumui = 0, sumvr = 0, sumvi = 0,
	    sumwr = 0, sumwi = 0;
	    for (int i = 0; i < N; i++) {
                double cosf =  Math.cos(-2 * Math.PI * f * i);
                double sinf =  Math.sin(-2 * Math.PI * f * i);
		double ur = real[i]*cosf - imag[i]*sinf;
		double ui = imag[i]*cosf + real[i]*sinf;
		double vr = 2 * Math.PI * i * ui;
		double vi = -2 * Math.PI * i * ur;
		double wr = 2 * Math.PI * i * vi;
		double wi = -2 * Math.PI * i * vr;
		sumur += ur;
		sumui += ui;
		sumvr += vr;
		sumvi += vi;
		sumwr += wr;
		sumwi += wi;
	    }
	    p = sumur * sumur + sumui * sumui;
	    if (p < lastp)  //I am not sure this is necessary, Vaughan did it for period estimation.
		f = (f + lastf) / 2;
	    else {
		lastf = f;
		lastp = p;
		if (p > maxp) {
		    maxp = p;
		    fhat = f;
		}
		pd = 2 * (sumvr * sumur + sumvi * sumui);
		pdd = 2 * (sumvr * sumvr + sumwr * sumur
			   + sumvi * sumvi + sumwi * sumui);
		f += pd / Math.abs(pdd);
	    }
	    numIter++;
        }
        
        //System.out.println("f = " + fhat);
        
        return fhat;
    }
    
}
