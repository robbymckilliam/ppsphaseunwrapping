package pubsim.poly;

import flanagan.complex.Complex;
import pubsim.optimisation.Brent;
import pubsim.optimisation.SingleVariateFunction;
import static pubsim.VectorFunctions.prod;

/**
 * The product high order ambiguity function estimator. See paper: Sergio
 * Barbarossa, Anna Scaglione and Georgios B. Giannakis. "Product High-Order
 * Ambiguity Function for Multicomponent Polynomial-Phase Signal Modeling" IEEE
 * TRANSACTIONS ON SIGNAL PROCESSING, VOL. 46, NO. 3, MARCH 1998
 * 
 * This doesn't use the FFT, see below PHAF.PHAFfft for this.
 *
 * TO DO: Implement lower order parameters, only highest order parameter
 * estimated at the moment.
 *
 * @author Robby McKilliam
 */
public class PHAF extends AbstractPolynomialPhaseEstimator {

    /**
     * The HAF estimators we are going to use
     */
    final HAF[] hafs;
    /**
     * The number of products, i.e. individual HAFs
     */
    final int L;
    /**
     * Frequency scalers for rescaling HAF in frequency so that they align
     */
    final double[] freqscaler;
    /**
     * Oversampling factor for the objective function, n*OS sample will be
     * computed
     */
    final double OS;
    final protected int n;
    final Complex[] z;

    /**
     * OS defaults to 4
     */
    public PHAF(int m, int n, int tau[][]) {
        this(m, n, tau, 4.0);
    }

    /**
     * Construct PHAF estimator. Takes an array of arrays of lags, so each
     * tau[i] is is the set of lags for one of the HAF estimators. OS is how
     * hard the objective function will be oversampled, i.e. n*OS samples will
     * be computed.
     */
    public PHAF(int m, int n, int tau[][], double OS) {
        super(m);
        this.OS = OS;
        this.n = n;
        z = new Complex[n];
        L = tau.length;
        hafs = new HAF[L];
        for (int i = 0; i < L; i++) {
            hafs[i] = new HAF(m, n, tau[i]);
        }
        freqscaler = new double[L];
        double p1 = prod(tau[0]);
        for (int i = 0; i < L; i++) {
            freqscaler[i] = prod(tau[i]) / p1; //compute all the frequency scalers
        }
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        for (int i = 0; i < n; i++) {
            z[i] = new flanagan.complex.Complex(real[i], imag[i]);
        }
        double[] p = new double[m + 1];
        p[m] = estimateM(z, m); //only does the highest order parameter at the moment;
        return p;
    }

    protected double estimateM(Complex[] x, int M) {
        double fstep = 1.0 / (n * OS);
        double bestO = Double.NEGATIVE_INFINITY;
        double fhat = 0.0;
        for (double f = 0.0; f > -1.0; f -= fstep) {
            double thisO = calculateObjective(f, z, m);
            if (bestO < thisO) {
                fhat = f;
                bestO = thisO;
            }
        }
        //refine using Brent's method
        SingleVariateFunction func = new SingleVariateFunction() {
            public double value(double x) {
                return -calculateObjective(x, z, m);
            }
        };
        Brent opt = new Brent(func, fhat - fstep, fhat, fhat + fstep);
        fhat = opt.xmin();
        fhat -= Math.round(fhat);
        return fhat / hafs[0].normaliser(M);
    }

    /**
     * Returns the product of the HAF objective functions
     */
    public double calculateObjective(double f, Complex[] x, int M) {
        double prod = 1.0;
        for (int i = 0; i < L; i++) {
            prod = prod * hafs[i].calculateObjective(freqscaler[i] * f, x, M); //this is a herendously slow way to do this!
        }
        return prod;
    }

 /**
 * The product high order ambiguity function estimator. See paper: Sergio
 * Barbarossa, Anna Scaglione and Georgios B. Giannakis. "Product High-Order
 * Ambiguity Function for Multicomponent Polynomial-Phase Signal Modeling" IEEE
 * TRANSACTIONS ON SIGNAL PROCESSING, VOL. 46, NO. 3, MARCH 1998
 * 
 * Uses the FFT and is fast, but restricts the set of allowable lags.
 *
 * TO DO: Implement lower order parameters, only highest order parameter
 * estimated at the moment.
 * @author Robby McKilliam
 */
    public static class PHAFfft extends PHAF {

        public PHAFfft(int m, int n, int tau[][]) {
            super(m, n, tau);
            long prod1 = prod(hafs[0].gettau());
            for (int i = 1; i < L; i++) {
                if (prod(hafs[0].gettau()) != prod1) {
                    throw new RuntimeException("To compute PHAF with FFT the product of the lag sets must be the same");
                }
            }
        }

        @Override
        protected double estimateM(Complex[] x, int M) {
            Complex[] ft = computeFFTproduct(x, M);
            double fstep = 1.0 / ft.length;
            double maxp = 0;
            double fhat = 0.0;
            double f = 0.0;
            for (int i = 0; i < ft.length; i++) {
                double p = ft[i].squareAbs();
                if (p > maxp) {
                    maxp = p;
                    fhat = f;
                }
                f -= fstep;
            }
            //refine using Brent's method
            SingleVariateFunction func = new SingleVariateFunction() {
                public double value(double x) {
                    return -calculateObjective(x, z, m);
                }
            };
            Brent opt = new Brent(func, fhat - fstep, fhat, fhat + fstep);
            fhat = opt.xmin();
            fhat -= Math.round(fhat);
            return fhat / hafs[0].normaliser(M);
        }

        /**
         * Compute the product of HAFs using the FFT
         */
        public Complex[] computeFFTproduct(Complex[] x, int M) {
            Complex[] fftprod = hafs[0].FFTHAF(x, M);
            for (int i = 1; i < L; i++) {
                Complex[] fft = hafs[i].FFTHAF(x, M);
                for (int n = 0; n < fft.length; n++) {
                    fftprod[n] = fftprod[n].times(fft[n]);
                }
            }
            return fftprod;
        }
    }
}
