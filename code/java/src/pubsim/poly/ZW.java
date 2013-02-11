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

/**
 * An implementation of Zhou and Wang's Euclidean algorithm approach for increasing
 * the range of parameters for which the HAF/DPT applies.
 * @author Robby McKilliam
 */
public class ZW extends AbstractPolynomialPhaseEstimator {
    
    final protected HAF haf1, haf2;
    final protected Complex[] z;
    final protected double[] p;
    final protected int tau1, tau2;
    final protected int n;
    
    public ZW(int m, int n, int tau1, int tau2){
        super(m);
        if( gcd(tau1,tau2) != 0 ) throw new RuntimeException("tau1 and tau2 must be relatively prime");
        this.n = n;
        z = new Complex[n];
        p = new double[m+1];
        this.tau1 = tau1;
        this.tau2 = tau2;
        haf1 = new HAF(m,n,tau1);
        haf2 = new HAF(m,n,tau2);
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(n != real.length) throw new RuntimeException("Data length does not equal " + n);
        
        for (int i = 0; i < n; i++) z[i] = new Complex(real[i], imag[i]);

        for (int i = m; i >= 0; i--) {
            p[i] = zwresolve(i);
            for (int j = 0; j < z.length; j++) {
                double cs = Math.cos(-2.0 * Math.PI * p[i] * Math.pow(j + 1, i));
                double ss = Math.sin(-2.0 * Math.PI * p[i] * Math.pow(j + 1, i));
                z[j] = z[j].times(new Complex(cs, ss));
            }
        }
        
        return p;
    }

    /** Zhang and Wong's ambiguity resolver */
    protected double zwresolve(int i) {
        BigRational f1 = new BigRational(haf1.estimateM(z, i),30);
        BigRational f2 = new BigRational(haf2.estimateM(z, i),30); 
        BigInteger a = new BigInteger(Integer.toString(tau1,2),2).pow(i-1);
        BigInteger b = new BigInteger(Integer.toString(tau2,2),2).pow(i-1);
        BigInteger M = ((f1 * new BigRational(b)) - (f2 * new BigRational(a))).round(); //Zhou and Wang's M
        BigInteger[] t = extended_gcd(a,b);
        assert(t[0].compareTo(BigInteger.ONE)==0); //assert GCD is one
        BigInteger n2 = t[1];
        BigInteger n1 = t[2].negate();
        BigInteger q = ((new BigRational(n1* M) - f1) * new BigRational(BigInteger.ONE,a)).round();
        BigInteger k1star = n1*M + q * a;
        BigRational phat = (new BigRational(k1star) + f1) / new BigRational(a * new BigInteger(Integer.toString(2)).pow(i-1) * new BigInteger(Long.toString(factorial(i))));
        return phat.doubleValue();
    }
    
}
