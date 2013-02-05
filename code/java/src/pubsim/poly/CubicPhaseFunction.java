/*
 */
package pubsim.poly;

import Jama.Matrix;
import pubsim.Complex;

/**
 * Implements O'Shea's cubic phase function estimator
 * Peter O'Shea "A Fast Algorithm for Estimating the Parameters of a Quadratic FM Signal"
 * IEEE Trans. Signal Proc. Vol 53, Feb 2004.
 * @author Robby McKilliam
 */
public class CubicPhaseFunction extends AbstractPolynomialPhaseEstimator {

    final int N; 
    final int samples;
    final protected double[] real;
    final protected double[] imag;
    final protected Matrix T; //transformation between polynomial bases
    
    /** Cubic phase function estimator of dimension N.  Will use sample the CP function
     * samples times to approximate maximum. */
    public CubicPhaseFunction(int N, int samples){
        super(3);
        if(N%2 == 0) throw new RuntimeException("n must be odd be the cubic phase function");
        this.N = N;
        this.samples = samples;
        real = new double[N];
        imag = new double[N];
        T = new Matrix(m+1,m+1);
        for(int t = 0; t <=m; t++) {
            for(int i = 0; i <= t; i++) {
                int k = -(N+1)/2;
                T.set(i,t, pubsim.Util.binom(t, i) * Math.pow(k, t-i) / 2 / Math.PI );          
            }
        }
    }
    
    /** Cubic phase function estimator of dimension N.  Will use sample the CP function
     * 4N times to approximate maximum. */
    public CubicPhaseFunction(int N){
        this(N,4*N);
    } 

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(real.length != N || imag.length != N)
            throw new ArrayIndexOutOfBoundsException();
        System.arraycopy(real, 0, this.real, 0, N); //copy input to local arrays
        System.arraycopy(imag, 0, this.imag, 0, N);
        
        //maximise the CP in in the two places recomended   
        int n1 = 0;
        int n2 = (int)Math.round(0.11*N); //presumably O'Shea means to round here
        double w1 = maxCP(n1);
        double w2 = maxCP(n2); 
        
        //now invert a 2x2 matrix to get quadratic and cubic estimates
        double det = 2*6*n2 - 2*6*n1;
        double a2 = (6*n2*w1 - 6*n1*w2)/det;
        double a3 = (-2*w1 + 2*w2)/det;
        
        //TO DO: use periodogram to get frequency and phase
        
        return transformToStandardBasis(new double[] {0.0,0.0,a2,a3});
    }  
    
    /** Re-centers indices to follow notation in O'Shea's paper */
    final protected Complex z(int n) {
        //index into array so that n appears centered
        int nn = n + (N-1)/2;
        if(nn < 0 || nn >= N) return Complex.zero;
        return new Complex(real[nn], imag[nn]);
    }
    
    /** The cubic phase function */
    final protected Complex CP(int n, double w){
        Complex sum = Complex.zero;
        for(int m = 0; m <= (N-1)/2; m++)
            sum = sum + (z(n+m) * z(n-m) * Complex.polar(1, -w*m*m));
        return sum;
    }
    
    /** 
     * Computes the maximum of CP over w for given n.  Uses 4*N samples.
     * O'Shea suggests that this can be done faster with a type of 'Fourier transform'.  He doesn't 
     * describe it or give a reference.  I have just implemented a direct sampling approach.
     */
    protected double maxCP(int n) {
        //range for w. For some reason O'Shea doesn't specify this. I've guessed it from
        //the ambiguity requirements he states for the estimator.
        double range = 2*(Math.PI/N + 3*Math.PI/2/N/N*Math.abs(n)); 
        double step = 2*range/samples;
        
        double CPbest = Double.NEGATIVE_INFINITY;
        double what = -range;
        for(double w = -range; w <= range; w+=step){
            double CPthis = CP(n,w).abs();
            if(CPthis > CPbest){
                CPbest = CPthis;
                what = w;
            }
            //System.out.print(w + ", " + CPthis + "; ");
        }
        //TO DO: REFINE what
        //System.out.println();
        return what;     
    }
    
    /** 
     * Transform parameters from O'Shea's origin centered basis to the 
     * standard polynomial basis.  The standard basis is the one used in my paper.
     */
    public double[] transformToStandardBasis(double[] p){
        return pubsim.VectorFunctions.matrixMultVector(T,p);
    }
 
    /** 
     * Transform parameters from the standard basis to O'Shea's origin centered 
     * polynomial basis.
     */
    public double[] transformToOriginCenterBasis(double[] p){
        return pubsim.VectorFunctions.matrixMultVector(T.inverse(),p);
    }
    
}
