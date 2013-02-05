/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pubsim.poly;

import Jama.Matrix;
import pubsim.VectorFunctions;
import pubsim.lattices.NearestPointAlgorithmInterface;
import pubsim.lattices.VnmStar;
import pubsim.lattices.VnmStarGlued;
import pubsim.lattices.decoder.Babai;

/**
 * Uses the Babai nearest plane algorithm
 * @author Robby
 */
public class BabaiEstimator implements PolynomialPhaseEstimator {

    protected double[] ya,  p;
    protected int n,  m;
    protected VnmStar lattice;
    protected NearestPointAlgorithmInterface npalgorithm;
    protected Matrix M,  K;
    protected AmbiguityRemover ambiguityRemover;
    
    protected BabaiEstimator() {}
    
    /** 
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public BabaiEstimator(int m, int n) {
        lattice = new VnmStarGlued(m, n - m - 1);
        this.m = m;
        npalgorithm = new Babai(lattice);
        ambiguityRemover = new AmbiguityRemover(m);

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


    /**
     * Run the estimator and return the square error wrapped modulo the the ambiguity
     * region between the estimate and the truth. i.e. dealias the estimate before
     * computing the square error.
     */
    public double[] error(double[] real, double[] imag, double[] truth){
        
        double[] est = estimate(real, imag);
        double[] err = new double[est.length];
        
        for (int i = 0; i < err.length; i++) {
            err[i] = est[i] - truth[i];
        }
        err = ambiguityRemover.disambiguate(err);
        for (int i = 0; i < err.length; i++) {
            err[i] *= err[i];
        }
        return err;
    }

    public int getOrder() {
        return m;
    }


}
