/*
 */
package pubsim.poly;

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
    protected double[] real;
    protected double[] imag;
    
    /** Cubic phase function estimator of dimension N.  Will use sample the CP function
     * samples times to approximate maximum. */
    public CubicPhaseFunction(int N, int samples){
        super(3);
        if(N%2 == 0) throw new RuntimeException("n must be odd be the cubic phase function");
        this.N = N;
        this.samples = samples;
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
        this.real = real; this.imag = imag; //I need to work out a better ways to do this.
        
        //maximise the CP in in the two places recomended
        int n1 = 0;
        int n2 = (int)Math.round(0.11*N); //presumably OShea means to round here
        double w1 = maxCP(n1);
        double w2 = maxCP(n2); 
        
        //now invert a 2x2 matrix to get quadratic and cubic estimates
        double det = 2*6*n2 - 2*6*n1;
        double a2 = (6*n2*w1 - 6*n1*w2)/det;
        double a3 = (-2*w1 + 2*w2)/det;
        
        //TO DO: use periodogram to get frequency and phase
        
        return new double[] {0.0, 0.0, a2/2/Math.PI, a3/2/Math.PI};
    }  
    
    /** Re-centers indices to follow notation in O'Shea's paper */
    final protected Complex z(int n) {
        //index into array so that n appears centered
        int nn = n + (N-1)/2;
        if(nn < 0 || nn >= N) return Complex.zero;
        return new Complex(real[nn], imag[nn]);
    }
    
    /** The cubic phase function */
    protected Complex CP(int n, double w){
        Complex sum = Complex.zero;
        for(int m = 0; m <= (N-1)/2; m++)
            sum = sum + (z(n+m) * z(n-m) * (Complex.polar(1, -w*m*m)));
        return sum;
    }
    
    /** 
     * Computes the maximum of CP over w for given n.  Uses 4*N samples.
     * OShea suggests that this can be done fasterwith a type of 'Fourier transform'.  He doesn't 
     * describe it or give a reference.  I have just implemented a direct sampling approach.
     */
    protected double maxCP(int n) {
        //range for w. For some reason Oshea doesn't specify this. I've guess it from
        //the ambiguity requirements he states for the estimator.
        double range = 2*(Math.PI/N + 3*Math.PI/2/N/N*n); 
        double step = 2*range/samples;
        
        double CPbest = Double.NEGATIVE_INFINITY;
        double what = -range;
        for(double w = -range; w <= range; w+=step){
            double CPthis = CP(n,w).abs();
            if(CPthis > CPbest){
                CPbest = CPthis;
                what = w;
            }
        }
        //TO DO: REFINE what
        return what;     
    }
 
}
