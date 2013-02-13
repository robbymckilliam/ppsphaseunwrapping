/*
 */
package pubsim.poly;

import Jama.Matrix;
import flanagan.complex.Complex;
import flanagan.math.FourierTransform;
import static pubsim.Util.factorial;

/**
 * Combines the HAF and CPF, i.e. does phase differencing to get to a cubic phase signal
 * and then uses the CPF.  See:
 * I. Djurovic and M. Simeunovic and S. Djukanovic and P. Wang 
 * "A Hybrid CPF-HAF Estimation of Polynomial-Phase Signals: Detailed Statistical Analysis"
 * IEEE Transactions on Signal Processing, 60:10, Oct. 2012.
 * TO DO: Estimate m-2 lower order parameters after dechirping
 * @author Robby McKilliam
 */
public class CPFHAF extends HAF {
    
    protected final CPF cpf;
    protected final double[] realcpf;
    protected final double[] imagcpf;
    final protected Matrix T, Tinv; //transformation between polynomial bases
     
    /** tau is set to default round(n/m) */
    public CPFHAF(int m, int n) {
        this(m,n,(n+1)/2/m); //tau default to half usual because of the symmetric definition of PPT here.
    }

    public CPFHAF(int m, int n, int tau) {
        super(m,n,tau);
        if(n%2 == 0) throw new RuntimeException("n must be odd be the cubic phase function");
        int cpfn = n - 2*tau*(m-3); //length adjusted for size of CPF after apply PPT m-3 times.
        cpf = new CPF(cpfn); 
        realcpf = new double[cpfn];
        imagcpf = new double[cpfn];
        T = CPF.constructOsheaBasisTransformaton(m,n);
        Tinv = T.inverse();
    }
    
    @Override
    public double[] estimate(double[] real, double[] imag) {
        for (int i = 0; i < n; i++) z[i] = new Complex(real[i], imag[i]); //complex version of recieved signal
        Complex[] d = PPT(m-2, z, tau); //compute the PPT to get cubic phase signal
        for(int i = 0; i < d.length; i++){
            realcpf[i] = d[i].getReal();
            imagcpf[i] = d[i].getImag();
        }
        double [] pcpf = cpf.estimateInOsheaBasis(realcpf, imagcpf);
        double C2 = Math.pow(2,m-4) * Math.pow(tau, m-3) * factorial(m-1);
        double C3 = C2 * m / 3.0;
        p[m] = pcpf[3]/C3;
        p[m-1] = pcpf[2]/C2;
        //TO DO:  Estimate remaining parameters after dechirping
        return transformToStandardBasis(p);
    }
    
    /** Djorovic's symmetric version of PPT */
    @Override
    protected Complex[] PPT2(Complex[] y, int tau) {
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
    
}
