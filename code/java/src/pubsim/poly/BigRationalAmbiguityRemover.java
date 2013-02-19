/*
 */
package pubsim.poly;

import Jama.Matrix;
import bignums.BigRational;
import pubsim.VectorFunctions;
/**
 * Same as AmbiguityRemoverRectangular, but uses infinite precision rational numbers
 * to avoid numerical problems that can occur with polynomial phase signals of high order
 * INCOMPLETE!
 * @author Robby McKilliam
 */
public class BigRationalAmbiguityRemover extends AmbiguityRemover {
    
    BigRational[][] Mat; 
    
    public BigRationalAmbiguityRemover(int m){
        this.m = m;
        p = new double[m+1];
    }
    
    /**
     * Matrix containing coefficients of the integer valued polynomials
     * @return 
     */
    protected final BigRational[][] constructBasisMatrixBigRational() {
        BigRational[][] B = new BigRational[m+1][m+1];
        for (int j = 0; j <= m; j++) for (int i = 0; i <= m; i++) B[i][j] = BigRational.ZERO;
        BigRational[] c = {BigRational.ONE};
        for (int j = 0; j < m+1; j++) {
            for (int i = 0; i <= j; i++) {
                B[i][j] = c[i];
            }
            c = getNextColumnBigRational(c);
        }
        return B;
    }
    
    /**
     * Recursively generates the columns of the ambiguity lattice
     * generator matrix.  Essentially generates coefficients of the
     * integer valued polynomials.
     * @param c previous column.
     * @return next column
     */
    protected static BigRational[] getNextColumnBigRational(BigRational[] c) {
        BigRational[] r = new BigRational[c.length + 1];
        for(int i = 0; i < r.length; i++) r[i] = BigRational.ZERO;
        for (int i = 1; i < r.length; i++) r[i] = c[i - 1]; //shift c
        for (int i = 0; i < c.length; i++) r[i] = r[i] + (c[i] * new BigRational(c.length - 1,1));
        for (int i = 0; i < r.length; i++) r[i] = r[i] / new BigRational(c.length,1);
        return r;
    }
    
}
