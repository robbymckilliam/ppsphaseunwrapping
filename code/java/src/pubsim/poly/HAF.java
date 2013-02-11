/*
 */
package pubsim.poly;

import flanagan.complex.Complex;
import flanagan.math.FourierTransform;
import pubsim.Util;

/**
 * Implementation of the estimator based on iteratively maximising
 * the Discrete Polynomial Transform.  See Peleg and Friedlander,
 * Trans. Sig. Proc. Vol 43 August 1995.
 * @author Robby McKilliam
 */
public class HAF extends AbstractPolynomialPhaseEstimator {

    final protected Complex[] z;
    final protected double[] p;
    final protected int n;
    final protected int tau;
    protected int num_samples;
    protected FourierTransform fft;
    protected Complex[] sig;
    
    /**Max number of iterations for the Newton step */
    static final int MAX_ITER = 12;
    /**Step variable for the Newton step */
    static final double EPSILON = 1e-12;

    /** tau is set to default round(n/m) */
    public HAF(int m, int n) {
        this(m,n,(int)Math.round(((double) n)/m));
    }

    public HAF(int m, int n, int tau) {
        super(m);
        z = new Complex[n];
        p = new double[m+1];
        this.n = n;
        num_samples = 4 * n;
        this.tau = tau;
        int oversampled = 4;
        sig = new Complex[FourierTransform.nextPowerOfTwo(oversampled * n)];
        fft = new FourierTransform();
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(n != real.length) throw new RuntimeException("Data length does not equal " + n);
        
        for (int i = 0; i < n; i++) z[i] = new Complex(real[i], imag[i]);

        for (int i = m; i >= 0; i--) {
            p[i] = estimateM(z, i);
            for (int j = 0; j < z.length; j++) {
                double cs = Math.cos(-2.0 * Math.PI * p[i] * Math.pow(j + 1, i));
                double ss = Math.sin(-2.0 * Math.PI * p[i] * Math.pow(j + 1, i));
                z[j] = z[j].times(new Complex(cs, ss));
            }
        }
        
        return p;
    }

    /**
     * Estimate the parameter of order M from x.
     * This uses m coarse (discrete) search of the periodogram and
     * then Newton Raphson to climb to the peak.  This is equivalent
     * to the periodogram method of frequency estimation.
     * @param x the signal
     * @param M order of parameter to estimate
     * @return the estimated parameter
     */
    protected double estimateM(Complex[] x, int M) {
        if(M != 0) return estimateMUnormalised(x,M,tau) / Util.factorial(M) / Math.pow(tau, M - 1);
        else return estimateMUnormalised(x,M,tau);
    }
    protected double estimateMUnormalised(Complex[] x, int M, int tau) {

        //compute the PPT
        Complex[] d = PPT(M, x, tau);

        //this is the phase parameter, so just average
        if (M == 0) {
            Complex zsum = new Complex(0,0);
            for(int i = 0; i < x.length; i++) zsum = zsum.plus(x[i]);
            return zsum.argRad() / (2 * Math.PI);
        }

        for (int i = 0; i < d.length; i++) sig[i] = new Complex(d[i]);
        for (int i = d.length; i < sig.length; i++) sig[i] = new Complex(0.0, 0.0);

        fft.setData(sig);
        fft.transform();
        Complex[] ft = fft.getTransformedDataAsComplex();

        //note that the FFT is generally defined with exp(-jw) but
        //periodogram has exp(jw) so freq are -ve here.
        double maxp = 0;
        double fhat = 0.0;
        double f = 0.0;
        double fstep = 1.0 / ft.length;
        for (int i = 0; i < ft.length; i++) {
            double p = ft[i].squareAbs();
            if (p > maxp) {
                maxp = p;
                fhat = f;
            }
            f-=fstep;
        }
        fhat = refine(fhat, d, maxp);
        
        return fhat;
    }

    /**
     * Compute the polynomial phase transform of order m of z
     * with lag tau.
     */
    protected Complex[] PPT(int m, Complex[] y){
        return PPT(m,y,tau);
    }
    protected Complex[] PPT(int m, Complex[] y, int tau) {
        Complex[] trans = y;
        for (int i = 2; i <= m; i++) trans = PPT2(trans, tau);
        return trans;
    }

    /**
     * Compute the second order PPT.  This is used to
     * compute the PPT of higher orders.
     * @param y
     */
    protected Complex[] PPT2(Complex[] y) {
        return PPT2(y,tau);
    }
    protected Complex[] PPT2(Complex[] y, int tau) {
        Complex[] trans = new Complex[y.length - tau];
        for (int i = 0; i < y.length; i++) {
            if (i - tau >= 0) {
                trans[i - tau] = y[i].times(y[i - tau].conjugate());
            }
        //else {
        //    trans[i] = new Complex(0, 0);
        //}
        }
        return trans;
    }

    /** Refine coarse frequency estimate using Newton Raphson method */
    protected double refine(double fhat, Complex[] d, double maxp) {
        double f;
        int numIter = 0;
        f = fhat;
        double lastf = f - 2 * EPSILON;
        double lastp = 0;
        while (Math.abs(f - lastf) > EPSILON && numIter <= MAX_ITER) {
            double p = 0, pd = 0, pdd = 0;
            double sumur = 0, sumui = 0, sumvr = 0, sumvi = 0,
                    sumwr = 0, sumwi = 0;
            for (int i = 0; i < d.length; i++) {
                double cosf = Math.cos(-2 * Math.PI * f * i);
                double sinf = Math.sin(-2 * Math.PI * f * i);
                double ur = d[i].getReal() * cosf - d[i].getImag() * sinf;
                double ui = d[i].getImag() * cosf + d[i].getReal() * sinf;
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
            if (p < lastp) //I am not sure this is necessary, Vaughan did it for period estimation.
            {
                f = (f + lastf) / 2;
            } 
            else {
                lastf = f;
                lastp = p;
                if (p > maxp) {
                    maxp = p;
                    fhat = f;
                }
                pd = 2 * (sumvr * sumur + sumvi * sumui);
                pdd = 2 * (sumvr * sumvr + sumwr * sumur + sumvi * sumvi + sumwi * sumui);
                f += pd / Math.abs(pdd);
            }
            numIter++;
        }
        //System.out.println("fine fhat = " + fhat);
        fhat -= Math.round(fhat);
        return fhat;
    }
    
    /**
     * Returns the volume of the functional regions of
     * the DPT estimator.  This is just for comparison
     * with the identifiable region.
     * @return volume of functional region
     */
    public static double volumeOfFunctionalRegion(double tau, int m){
        double prod = 1.0;
        for(int k = 2; k <= m; k++) prod *= 1.0/( Util.factorial(k) * Math.pow(tau, k-1) );
        return prod;
    }

}
