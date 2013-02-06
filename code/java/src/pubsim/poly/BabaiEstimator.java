/*
 */
package pubsim.poly;

import Jama.Matrix;
import pubsim.VectorFunctions;
import pubsim.lattices.NearestPointAlgorithmInterface;
import pubsim.lattices.VnmStar;
import pubsim.lattices.decoder.Babai;

/**
 * Uses the Babai nearest plane algorithm
 * @author Robby
 */
public class BabaiEstimator extends AbstractPolynomialPhaseEstimator {

    final protected double[] ya,  p;
    final protected int n;
    final protected VnmStar lattice;
    protected NearestPointAlgorithmInterface npalgorithm;
    final protected Matrix M,  K;
    
    /** 
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public BabaiEstimator(int m, int n) {
        super(m);
        lattice = new VnmStar(m, n - m - 1);
        npalgorithm = new Babai(lattice);
        ya = new double[n];
        p = new double[m+1];
        this.n = n;
        M = lattice.getMMatrix();
        Matrix Mt = M.transpose();
        K = Mt.times(M).inverse().times(Mt);
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(n != real.length) throw new RuntimeException("Data length does not equal " + n);
        
        for (int i = 0; i < real.length; i++) {
            ya[i] = Math.atan2(imag[i], real[i]) / (2 * Math.PI);
        }
        npalgorithm.nearestPoint(ya);
        double[] u = npalgorithm.getIndex();

        double[] ymu = new double[ya.length];
        for (int i = 0; i < u.length; i++) {
            ymu[i] = ya[i] - u[i];
        }
        System.arraycopy(ya, u.length, ymu, u.length, ya.length - u.length);

        //compute the parameters
        VectorFunctions.matrixMultVector(K, ymu, p); 
        
        return ambiguityRemover.disambiguate(p);
    }

}
