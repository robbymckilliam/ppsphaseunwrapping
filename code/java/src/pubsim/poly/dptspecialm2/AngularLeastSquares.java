/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly.dptspecialm2;

import Jama.Matrix;
import pubsim.lattices.Anstar.Anstar;
import pubsim.lattices.VnmStar;
import pubsim.Complex;
import pubsim.VectorFunctions;
import pubsim.lattices.Anstar.AnstarLinear;
import pubsim.poly.AmbiguityRemover;
import pubsim.poly.PolynomialPhaseEstimator;

/**
 * Returns the angular least squares estimator but only searches in the interval
 * on which the DPT estimator works.
 * @author Robby McKilliam
 */
public class AngularLeastSquares implements PolynomialPhaseEstimator{

    protected final int m = 2;
    protected int N;
    protected Matrix M,  K;

    private Anstar anstar;
    private double[] phi;
    private double[] z;
    private double[] pmw;
    private final AmbiguityRemover ambiguityRemover;

    public AngularLeastSquares(){
        ambiguityRemover = new AmbiguityRemover(m);
    }

    public void setSize(int n) {
        N = n;
        anstar = new AnstarLinear(N-1);
        M = VnmStar.getMMatrix(N-m-1, m);
        Matrix Mt = M.transpose();
        K = Mt.times(M).inverse().times(Mt);
        z = new double[N];
        pmw = new double[N];
        phi = new double[N];
    }

    public int getOrder() {
        return m;
    }

    public double[] estimate(double[] real, double[] imag) {
        if(real.length != N) setSize(real.length);

        final int numsamples = 4*N;

        for(int i = 0; i < N; i++){
            phi[i] = (new Complex(real[i], imag[i])).phase()/(2*Math.PI);
        }

        double p2start = -(2.0/N)/4.0;
        double p2end = (2.0/N)/4.0;
        double step = (p2end - p2start)/numsamples;

        double D = Double.POSITIVE_INFINITY;
        double[] hatp = new double[3];

        for(double p2 = p2start; p2 < p2end; p2 += step){
            for(double f = -0.5; f < 0.5; f+= 1.0/numsamples){

                for (int i = 0; i < N; i++) {
                    double phs = p2*(i+1)*(i+1) + f*(i+1);
                    z[i] = phi[i] - phs;
                }

                anstar.nearestPoint(z);
                double[] w = anstar.getIndex();

                for(int i = 0; i < N; i++) pmw[i] = phi[i] - w[i];

                double[] p = VectorFunctions.matrixMultVector(K, pmw);
                VnmStar.project(pmw, pmw, m);
                double dist = VectorFunctions.sum2(pmw);
                if(dist < D){
                    D = dist;
                    System.arraycopy(p, 0, hatp, 0, 3);
                }

            }
            //System.out.println(VectorFunctions.print(hatp));
        }

        return ambiguityRemover.disambiguate(hatp);
        
    }

    /**
     * Run the estimator and return the square error wrapped modulo the the ambiguity
     * region between the estimate and the truth.
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

}
