/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.fes;

import Jama.Matrix;
import pubsim.lattices.Vn1Star.Vn1Star;

/**
 * Abstract class for the frequency estimators that use the Vn1Star
 * lattice.
 * @author Robby McKilliam
 */
public abstract class LatticeEstimator implements FrequencyEstimator{

    protected final int N;
    protected Vn1Star lattice;
    protected final double[] y;
    //protected Matrix K;

    protected LatticeEstimator() { throw new UnsupportedOperationException(); }

    public LatticeEstimator(int N){
        this.N = N;
        y = new double[N];
    }

        /** Run the estimator on received data, @param y */
    @Override
    public double estimateFreq(double[] real, double[] imag){
        for(int i = 0; i < real.length; i++)
            y[i] = Math.atan2(imag[i],real[i])/(2*Math.PI);

        lattice.nearestPoint(y);

        //calculate f from the nearest point
        double f = 0;
        double N = this.N;
        double sumn = N*(N+1)/2;
        double sumn2 = N*(N+1)*(2*N+1)/6;
        double[] u = lattice.getIndex();
        for(int i = 0; i < this.N; i++)
            f += (N*(i+1) - sumn)*(y[i]-u[i]);

        f /= (sumn2*N - sumn*sumn);

        //System.out.println("f = " + f);

        return f - Math.round(f);

    }

}
