/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pubsim.poly;

import bignums.BigInteger;
import bignums.BigRational;
import static pubsim.Util.extended_gcd;
import static pubsim.Util.factorial;

/**
 * An implementation of Zhou and Wang's Euclidean algorithm approach for increasing
 * the range of parameters for which the HAF/DPT applies.  This implementation is slightly
 * different from Zhou and Wang's.
 * @author Robby McKilliam
 */
public class ZW extends HAFSimpleDiversity {
    
    public ZW(int m, int n, int tau1, int tau2) {
        super(m,n,tau1,tau2);
    }
    
    /** Zhang and Wong's ambiguity resolver */
    @Override
    protected double resolve(int P) {
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
    
}
