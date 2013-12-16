/*
 */
package pubsim.poly;

import Jama.Matrix;
import pubsim.Complex;
import pubsim.optimisation.FunctionAndDerivatives;
import pubsim.optimisation.NewtonRaphson;

/**
 * Implements O'Shea's cubic phase function estimator
 * Peter O'Shea "A Fast Algorithm for Estimating the Parameters of a Quadratic FM Signal"
 * IEEE Trans. Signal Proc. Vol 53, Feb 2004.
 * @todo Dechirp and estimate phase and frequency, only quadratic and cubic term currently estimated.
 * @author Robby McKilliam
 */
public class CPF extends AbstractPolynomialPhaseEstimator {

    final int N; 
    final int samples;
    protected double[] real;
    protected double[] imag;
    final protected double[] wreal; //working memory
    final protected double[] wimag; //working memory
    final protected Complex[] y; //working memory
    final protected Matrix T, Tinv; //transformation between polynomial bases
    //final protected PeriodogramFFTEstimator fes;
 
    /** Cubic phase function estimator of dimension N.  Will use sample the CP function
     * samples times to approximate maximum. */
    public CPF(int N, int samples) {
        super(3);
        if(N%2 == 0) throw new RuntimeException("n must be odd be the cubic phase function");
        this.N = N;
        this.samples = samples;
        wreal = new double[N];
        wimag = new double[N];
        y = new Complex[N];
        T = constructOsheaBasisTransformaton(m,N);
        Tinv = T.inverse();
        //fes = new PeriodogramFFTEstimator(N);
    }
    
    /** Cubic phase function estimator of dimension N.  Will use sample the CP function
     * 10N times to approximate maximum. */
    public CPF(int N){
        this(N,10*N);
    } 

    @Override
    public double[] estimate(double[] real, double[] imag) {
        return transformToStandardBasis(estimateInOsheaBasis(real,imag));
    }  
    
    public double[] estimateInOsheaBasis(double[] real, double[] imag) {
        if(real.length != N || imag.length != N)
            throw new ArrayIndexOutOfBoundsException();
        this.real = real; this.imag = imag; //pointers to input signal
        
        //maximise the CP in in the two places recomended   
        int n1 = 0;
        int n2 = (int)Math.round(0.11*N); //presumably O'Shea means to round here
        double w1 = maxCP(n1);
        double w2 = maxCP(n2); 
        
        //now invert a 2x2 matrix to get quadratic and cubic estimates
        double det = 2*6*n2 - 2*6*n1;
        double a2 = (6*n2*w1 - 6*n1*w2)/det;
        double a3 = (-2*w1 + 2*w2)/det;
        
        //TO DO: Dechirp and estimate phase and frequency
//        dechirp(a2, a3); //subtract chirp signal
//        double freq = fes.estimateFreq(real, imag); //estimate frequency in standard basis
//        defreq(freq); //subtract frequency estimate
//        double phase = new SampleCircularMean().estimatePhase(y); //estimate phase in standard basis
//        
//        System.out.println(phase + ", " + freq);
//        
//        //get phase and freq in origin centered basis
//        double[] a = transformToOriginCenterBasis(new double[] {phase,freq,0,0}); 
  
        return new double[] {0.0,0.0,a2,a3};
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
    
    /** static verion of z for outsize use */
    public static Complex z(int n, Complex[] d){
        int N =d.length;
        int nn = n + (N-1)/2;
        if(nn < 0 || nn >= N) return Complex.zero;
        return d[nn];
    }
    
    /** static version of CP for outsize use */
    public static Complex CP(int n, double w, Complex[] d){
        int N = d.length;
        Complex sum = Complex.zero;
        for(int m = 0; m <= (N-1)/2; m++)
            sum = sum + (z(n+m,d) * z(n-m,d) * Complex.polar(1, -w*m*m));
        return sum;
    }
    
    /** 
     * Computes the maximum of CP over w for given n.  Uses 4*N samples.
     * O'Shea suggests that this can be done faster with a type of 'Fourier transform'.  He doesn't 
     * describe it or give a reference.  I have just implemented a direct sampling approach.
     */
    protected double maxCP(final int n) {
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
        //refine what
        FunctionAndDerivatives f = new FunctionAndDerivatives() {
            public double value(Matrix x) {
                return CP(n,x.get(0,0)).abs2();
            }
            public Matrix hessian(Matrix x) {
                Complex f = CP(n,x.get(0,0));
                Complex fd = CPdw(n,x.get(0,0));
                Complex fdd = CPdw2(n,x.get(0,0));
                double hes = 2*fd.abs() + 2*(fdd*f.conjugate()).re();
                return Matrix.identity(1,1).times(hes);
            }
            public Matrix gradient(Matrix x) {
                Complex fd = CPdw(n,x.get(0,0));
                Complex f = CP(n,x.get(0,0));
                double grad = 2 * (fd * f.conjugate()).re();
                return Matrix.identity(1,1).times(grad);
            }
        };
        NewtonRaphson optimiser = new NewtonRaphson(f);
        return optimiser.maximise(Matrix.identity(1,1).times(what)).get(0,0);
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
        return pubsim.VectorFunctions.matrixMultVector(Tinv,p);
    }

//    /** Remove quadratic and cubic parameters in the origin center basis */
//    protected void dechirp(double a2, double a3) {
//        //subtract estimated quadratic and cubic signal
//        for(int n = (N-1)/2; n <= (N-1)/2; n++) {
//            double phase = n*n*a2 + n*n*n*a3;
//            int t = n + (N-1)/2;
//            Complex temp = new Complex(real[t],imag[t]) * Complex.polar(1,-phase);
//            wreal[t] = temp.re();
//            wimag[t] = temp.im();
//        }
//    }
//
//    /** Remove frequency in standard basis */
//    protected void defreq(double freq) {
//        for(int t = 0; t < N; t++) 
//            y[t] = new Complex(wreal[t],wimag[t]) * Complex.polar(1, -2*Math.PI*freq*(t+1));
//    }
    
    /** First derivative of the cubic phase function with respect to w*/
    final protected Complex CPdw(int n, double w){
        Complex sum = Complex.zero;
        for(int m = 0; m <= (N-1)/2; m++)
            sum = sum - (z(n+m) * z(n-m) * Complex.polar(1, -w*m*m) * (new Complex(0,m*m)));
        return sum;
    }
    
    /** Second derivative of the cubic phase function with respect to w*/
    final protected Complex CPdw2(int n, double w){
        Complex sum = Complex.zero;
        for(int m = 0; m <= (N-1)/2; m++)
            sum = sum - (z(n+m) * z(n-m) * Complex.polar(1, -w*m*m) * (new Complex(m*m*m*m,0)));
        return sum;
    }

    public static Matrix constructOsheaBasisTransformaton(int m, int N) {
        Matrix T = new Matrix(m+1,m+1);
        for(int t = 0; t <=m; t++) {
            for(int i = 0; i <= t; i++) {
                int k = -(N+1)/2;
                T.set(i,t, pubsim.Util.binom(t, i) * Math.pow(k, t-i) / 2.0 / Math.PI );          
            }
        }
        return T;
    }
    
}
