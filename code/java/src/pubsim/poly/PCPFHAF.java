/*
 * @author Robby McKilliam
 */
package pubsim.poly;

import Jama.Matrix;
import pubsim.Complex;
import static pubsim.Util.factorial;
import pubsim.VectorFunctions;
import static pubsim.VectorFunctions.prod;
import pubsim.optimisation.Brent;
import pubsim.optimisation.SingleVariateFunction;

/**
 * The product CPF HAF estimator.  Maximises an objective function constructed by a
 * product of CPFHAF estimators in a similar way to how the PHAF works.  See Section 5 of
 * I. Djurovic and M. Simeunovic and S. Djukanovic and P. Wang 
 * "A Hybrid CPF-HAF Estimation of Polynomial-Phase Signals: Detailed Statistical Analysis"
 * IEEE Transactions on Signal Processing, 60:10, Oct. 2012.
 * TO DO: Estimate m-2 lower order parameters after dechirping
 * @author Robby McKilliam
 */
public class PCPFHAF extends AbstractPolynomialPhaseEstimator {
    
    /**
     * The HAF estimators we are going to use
     */
    final CPFHAF[] cpfhafs;
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
    final protected int N;
    final Complex[] z;
    final int[][] tau;
    
    /** Oversample by 8 by default */ 
    public PCPFHAF(int m, int N, int tau[][]) {
        this(m,N,tau,8);
    }
    
     /**
     * Construct PCPFPHAF estimator. Takes an array of arrays of lags, so each
     * tau[i] is is the set of lags for one of the CPFHAF estimators. OS is how
     * hard the objective function will be oversampled, i.e. n*OS samples will
     * be computed.
     */
    public PCPFHAF(int m, int N, int tau[][], double OS) {
        super(m);
        this.OS = OS;
        this.N = N;
        this.tau = tau;
        z = new Complex[N];
        L = tau.length;
        cpfhafs = new CPFHAF[L];
        for (int i = 0; i < L; i++) {
            cpfhafs[i] = new CPFHAF(m, N, tau[i]);
        }
        freqscaler = new double[L];
        double p1 = prod(tau[0]);
        for (int i = 0; i < L; i++) {
            freqscaler[i] = prod(tau[i]) / p1; //compute all the frequency scalers
        }
    }


    @Override
    public double[] estimate(double[] real, double[] imag) {
        for (int i = 0; i < N; i++)
            z[i] = new Complex(real[i], imag[i]);
        //maximise the CP in in the two places recomended   
        int n1 = 0;
        int n2 = (int)Math.round(0.11*N); //presumably O'Shea means to round here
        double w1 = maxOmega(n1);
        double w2 = maxOmega(n2); 
        double C2 = Math.pow(2,m-4) * VectorFunctions.prod(tau[0]) * factorial(m-1);
        double C3 = C2 * m / 3.0;
        double[] p = new double[m + 1];
        p[m] = (w2-w1)/(6*n2*C3);
        p[m-1] = w1/(2*C2);
        //TO DO:  Estimate remaining parameters after dechirping
        return cpfhafs[0].transformToStandardBasis(p);
    }
    
    protected double maxOmega(final int n) {
        //range for w. For some reason O'Shea doesn't specify this. I've guessed it from
        //the ambiguity requirements he states for the estimator.
        double samples = N*OS;
        double range = 2*(Math.PI/N + 3*Math.PI/2/N/N*Math.abs(n)); 
        double step = 2*range/samples;
        
        double CPbest = Double.NEGATIVE_INFINITY;
        double what = -range;
        for(double w = -range; w <= range; w+=step){
            double CPthis = calculateObjective(n,w);
            if(CPthis > CPbest){
                CPbest = CPthis;
                what = w;
            }
        }
        //refine using Brent's method
        SingleVariateFunction func = new SingleVariateFunction() {
            public double value(double x) {
                return -calculateObjective(n,x);
            }
        };
        Brent opt = new Brent(func, what - step, what, what + step);
        return opt.xmin();
    }
    
    /**
     * Returns the product of the CPFHAF objective functions
     */
    final protected double calculateObjective(int n, double f) {
        double prod = 1.0;
        for (int i = 0; i < L; i++) {
            prod = prod * cpfhafs[i].calculateObjective(n, freqscaler[i] * f, z); //this is a herendously slow way to do this!
        }
        return prod;
    }
    
}
