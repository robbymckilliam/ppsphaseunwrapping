/*
 */
package pubsim.poly;

import Jama.Matrix;
import pubsim.VectorFunctions;
import pubsim.lattices.NearestPointAlgorithmInterface;
import pubsim.lattices.Vnmstar.HilbertMatrix;
import pubsim.lattices.Vnmstar.VnmStar;
import pubsim.lattices.reduction.LLL;
import pubsim.lattices.reduction.LatticeReduction;

/**
 * Uses the Babai nearest plane algorithm
 * @author Robby McKilliam
 */
public class Babai extends AbstractPolynomialPhaseEstimator {

    final protected double[] ya,  p;
    final protected int n;
    final protected VnmStar lattice;
    protected NearestPointAlgorithmInterface npalgorithm;
    final protected Matrix M,  K;
    
    public Babai(int m, int n){
        this(m,n,new LLL());
    } 
    
    /** 
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public Babai(int m, int n, LatticeReduction lr) {
        super(m);
        lattice = new VnmStar(m, n - m - 1);
        npalgorithm = new pubsim.lattices.decoder.Babai(lattice, lr);
        ya = new double[n];
        p = new double[m+1];
        this.n = n;
        M = lattice.getMMatrix();
        Matrix Mt = M.transpose();
        //System.out.println(Mt.times(M).inverse().cond());
        K = new HilbertMatrix(m+1,n).KDouble();
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
