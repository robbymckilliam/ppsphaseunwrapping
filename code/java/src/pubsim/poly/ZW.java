/*
 */
package pubsim.poly;

import flanagan.complex.Complex;
import static pubsim.Util.extended_gcd;
import static pubsim.Util.gcd;
import static pubsim.Util.factorial;
import bignums.BigInteger;
import bignums.BigRational;
import bignums.RoundingMode;
import pubsim.Util;

/**
 * An implementation of Zhou and Wang's Euclidean algorithm approach for increasing
 * the range of parameters for which the HAF/DPT applies. 
 * @author Robby McKilliam
 */
public class ZW extends AbstractPolynomialPhaseEstimator {
    
    final protected HAF haf;
    final protected Complex[] z;
    final protected double[] p;
    final protected int tau1, tau2;
    final protected int n;
    
    public ZW(int m, int n, int tau1, int tau2){
        super(m);
        if( gcd(tau1,tau2) != 1 ) throw new RuntimeException("tau1 and tau2 must be relatively prime");
        this.n = n;
        z = new Complex[n];
        p = new double[m+1];
        this.tau1 = tau1;
        this.tau2 = tau2;
        haf = new HAF(m,n,tau1);
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(n != real.length) throw new RuntimeException("Data length does not equal " + n);
        
        for (int i = 0; i < n; i++) z[i] = new Complex(real[i], imag[i]);

        for (int i = m; i >= 0; i--) {
            if(i>1) p[i] = zwresolve(i);
            else p[i] = haf.estimateMUnormalised(z, i, tau1);
            for (int j = 0; j < z.length; j++) {
                double cs = Math.cos(-2.0 * Math.PI * p[i] * Math.pow(j + 1, i));
                double ss = Math.sin(-2.0 * Math.PI * p[i] * Math.pow(j + 1, i));
                z[j] = z[j].times(new Complex(cs, ss));
            }
        }
        
        return p;
    }

    /** Zhang and Wong's ambiguity resolver */
    protected double zwresolve(int P) {
        double f1 = haf.estimateMUnormalised(z, P, tau1);
        double f2 = haf.estimateMUnormalised(z, P, tau2);
        BigInteger a = new BigInteger(Integer.toString(tau1)).pow(P-1);
        BigInteger b = new BigInteger(Integer.toString(tau2)).pow(P-1);
        BigInteger[] t = extended_gcd(a,b);
        BigRational n1 = new BigRational(t[1]);
        BigRational n2 = new BigRational(t[2]);
        BigRational phat = new BigRational(f1,30) * n1 +  new BigRational(f2,30) * n2;
        return pubsim.Util.fracpart(phat.doubleValue())/factorial(P);
//          return haf1.estimateM(z, i);
    }
    
    public static BigRational fracpart(BigRational x) {
        return x - new BigRational(x.round());
    }
    
}
