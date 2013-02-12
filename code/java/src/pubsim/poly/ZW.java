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
 * An implementation of Zhou and Wang's Euclidean algorithm approach for
 * increasing the range of parameters for which the HAF/DPT applies. This
 * implementation is slightly different from Zhou and Wang's.
 *
 * @author Robby McKilliam
 */
public class ZW extends HAFSimpleDiversity {

    public ZW(int m, int n, int tau1, int tau2) {
        super(m, n, tau1, tau2);
    }

    /**
     * Zhang and Wong's ambiguity resolver
     */
    @Override
    protected double resolve(int P) {
        double f1d = haf.estimateMUnormalised(z, P, tau1);
        double f2d = haf.estimateMUnormalised(z, P, tau2);
        BigRational f1 = new BigRational(f1d, 30);
        BigRational f2 = new BigRational(f2d, 30);
        BigInteger a = new BigInteger(Integer.toString(tau1)).pow(P - 1);
        BigInteger b = new BigInteger(Integer.toString(tau2)).pow(P - 1);
        BigInteger M = ((f1 * new BigRational(b)) - (f2 * new BigRational(a))).round();
        BigInteger[] t = extended_gcd(a, b);
        BigInteger n1 = t[2].negate();
        BigRational v = fracpart(new BigRational(n1*M) + f1, new BigRational(a));
        //System.out.println(new BigRational(a) + ", " + v.doubleValue() + ", " + (new BigRational(n1*M) + f1).doubleValue() + ", " + ((new BigRational(n1*M) + f1)/a).doubleValue() + ", " + ((new BigRational(n1*M) + f1)/a).round());
        BigRational phat = v / new BigRational(new BigInteger(Long.toString(factorial(P))) * a);
        return phat.doubleValue();
    }

    public static BigRational fracpart(BigRational x, BigRational a) {
        return x - (a * new BigRational((x/a).round()));
    }
}
