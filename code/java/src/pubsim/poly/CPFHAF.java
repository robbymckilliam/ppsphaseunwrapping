/*
 */
package pubsim.poly;

import Jama.Matrix;
import pubsim.Complex;
import flanagan.math.FourierTransform;
import static pubsim.Util.factorial;
import pubsim.VectorFunctions;

/**
 * Combines the HAF and CPF, i.e. does phase differencing to get to a cubic phase signal
 * and then uses the CPF.  See:
 * I. Djurovic and M. Simeunovic and S. Djukanovic and P. Wang 
 * "A Hybrid CPF-HAF Estimation of Polynomial-Phase Signals: Detailed Statistical Analysis"
 * IEEE Transactions on Signal Processing, 60:10, Oct. 2012.
 * TO DO: Estimate m-2 lower order parameters after dechirping
 * @author Robby McKilliam
 */
public class CPFHAF extends AbstractPolynomialPhaseEstimator {
    
    protected final CPF cpf;
    protected final double[] realcpf;
    protected final double[] imagcpf;
    final protected Matrix T, Tinv; //transformation between polynomial bases
    final protected Complex[] z;
    final protected double[] p;
    final protected int n;
    final protected int[] tau;
     
    /** tau is set to default round(n/m) */
    public CPFHAF(int m, int n) {
        this(m,n,(n+1)/2/m); //tau default to half usual because of the symmetric definition of PPT here.
    }

    public CPFHAF(int m, int n, int tau) {
        this(m,n,VectorFunctions.filledArray(m-3, tau));
    }
    
    public CPFHAF(int m, int n, int[] tau) {
        super(m);
        if(m <= 3) throw new RuntimeException("m must be greater than 3");
        if(n%2 == 0) throw new RuntimeException("n must be odd be the cubic phase function");
        this.tau = tau;
        z = new Complex[n];
        p = new double[m+1];
        this.n = n;
        int cpfn = n - 2*VectorFunctions.sum(tau); //length adjusted for size of CPF after apply PPT m-3 times.
        cpf = new CPF(cpfn); 
        realcpf = new double[cpfn];
        imagcpf = new double[cpfn];
        T = CPF.constructOsheaBasisTransformaton(m,n);
        Tinv = T.inverse();
    }
    
    @Override
    public double[] estimate(double[] real, double[] imag) {
        for (int i = 0; i < n; i++) z[i] = new Complex(real[i], imag[i]); //complex version of recieved signal
        Complex[] d = CPFHAF.PPT(m-3, z,tau); //compute the PPT to get cubic phase signal
        for(int i = 0; i < d.length; i++){
            realcpf[i] = d[i].re();
            imagcpf[i] = d[i].im();
        }
        double [] pcpf = cpf.estimateInOsheaBasis(realcpf, imagcpf);
        double C2 = Math.pow(2,m-4) * VectorFunctions.prod(tau) * factorial(m-1);
        double C3 = C2 * m / 3.0;
        p[m] = pcpf[3]/C3;
        p[m-1] = pcpf[2]/C2;
        //TO DO:  Estimate remaining parameters after dechirping
        return transformToStandardBasis(p);
    }
    
    /** All the lags are the same */
    protected static Complex[] PPT(int m, Complex[] y, int[] tau) {
        Complex[] trans = y;
        for (int i = 0; i < m; i++) trans = CPFHAF.PPT2(trans, tau[i]);
        return trans;
    }
    
    /** Djorovic's symmetric version of PPT */
    protected static Complex[] PPT2(Complex[] y, int tau) {
        int N = y.length;
        Complex[] trans = new Complex[N - 2*tau];
        for (int i = tau; i < N-tau; i++) 
                trans[i-tau] = y[i+tau].times(y[i - tau].conjugate());
        return trans;
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
    
    /**
     * Calculates the CPF objective function after applying the PPT m-3 times.
     * Only does highest order parameters!
     */
     public double calculateObjective(int n, double omega, Complex[] x){
         Complex[] d = CPFHAF.PPT(m-3,x,tau); //compute the PPT to get cubic phase signal
         return CPF.CP(n, omega, d).abs();
     }
  
    
}
