package pubsim.poly;

import flanagan.complex.Complex;
import pubsim.Util;
import pubsim.optimisation.Brent;
import pubsim.optimisation.SingleVariateFunction;

/**
 * The product high order ambiguity function estimator.  See paper:
 * Sergio Barbarossa, Anna Scaglione and Georgios B. Giannakis. 
 * "Product High-Order Ambiguity Function for Multicomponent Polynomial-Phase Signal Modeling"
 * IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 46, NO. 3, MARCH 1998
 * 
 * TO DO: Implement lower order parameters, only highest order parameter estimated at the moment.
 * @author Robby McKilliam
 */
public class PHAF extends AbstractPolynomialPhaseEstimator {
    
    /** The HAF estimators we are going to use */ 
    final HAF[] hafs;
    /** The number of products, i.e. individual HAFs */
    final int L; 
    /** Frequency scalers for rescaling HAF in frequency so that they align */
    final double[] freqscaler;
    /** Oversampling factor for the objective function, n*OS sample will be computed */
    final double OS;
    final protected int n;
    final Complex[] z;
    
    /** OS defaults to 4 */
    public PHAF(int m, int n, int tau[][]){
        this(m,n,tau,4.0);
    }
    
    /** 
     * Construct PHAF estimator.  Takes an array of arrays of lags, so each tau[i] is
     * is the set of lags for one of the HAF estimators.  OS is how hard the objective function
     * will be oversampled, i.e. n*OS samples will be computed.
     */
    public PHAF(int m, int n, int tau[][], double OS){
        super(m);
        this.OS = OS;
        this.n = n;
        z = new Complex[n];
        L = tau.length;
        hafs = new HAF[L];
        for(int i = 0; i < L; i++) hafs[i] = new HAF(m, n, tau[i]);
        freqscaler = new double[L];
        double p1 = pubsim.VectorFunctions.prod(tau[0]);
        for(int i = 0; i < L; i++) freqscaler[0] = pubsim.VectorFunctions.prod(tau[0]) / p1; //compute all the frequency scalers
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        for (int i = 0; i < n; i++) z[i] = new flanagan.complex.Complex(real[i], imag[i]);
        double[] p = new double[m+1];
        p[m] = estimateM(z, m); //only does the highest order parameter at the moment;
        return p;
    }
    
    protected double estimateM(Complex[] x, int M) {      
        double fstep = 1.0/(n*OS);
        double bestO = Double.NEGATIVE_INFINITY;
        double fhat = 0.0;
        for(double f = 0.0; f > -1.0; f-=fstep){ 
            double thisO = calculateObjective(f, z, m);
            if(bestO < thisO){
                fhat = f;
                bestO = thisO;
            }
        }
        //refine using Brent's method
        SingleVariateFunction func = new SingleVariateFunction() {
            public double value(double x) { return -calculateObjective(x,z,m); }
        };
        Brent opt = new Brent(func, fhat-fstep, fhat, fhat+fstep);
        fhat = opt.xmin();
        fhat -= Math.round(fhat);
        return fhat / hafs[0].normaliser(M);
    }
    
    /** Returns the product of the HAF objective functions */
    public double calculateObjective(double f, Complex[] x, int M) {
        double prod = 1.0;
        for(int i = 0; i < L; i++)
            prod = prod * hafs[i].calculateObjective(freqscaler[i]*f, x, M);
        return prod;       
    }
    
}
