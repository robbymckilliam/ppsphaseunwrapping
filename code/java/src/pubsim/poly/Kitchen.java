/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import pubsim.Util;
import pubsim.VectorFunctions;
import pubsim.bearing.SampleCircularMean;

/**
 * Implementation of Kitchen's polynomial phase estimate.
 *
 * J. Kitchen, "A method for estimating the coefficients of a polynomial
 * phase signal", Signal Processing, vol 37, 1994.
 * 
 * This is essentially a generalisation of Kay's frequency estimator.
 * @author Robby McKilliam
 */
public class Kitchen extends AbstractPolynomialPhaseEstimator{

    final protected Double[] y;
    final protected double[] p;
    final protected int N;

    public Kitchen(int m, int N) {
        super(m);
        this.N = N;
        y = new Double[N];
        p = new double[m+1];
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(N != real.length) throw new RuntimeException("Data length does not equal " + N);

        for(int i = 0; i < N; i++)
            y[i] = Math.atan2(imag[i], real[i])/(2*Math.PI);

        for (int i = m; i >= 0; i--) {
            p[i] = estimateM(y, i);
            for (int j = 0; j < y.length; j++) {
                y[j] -= p[i]*Math.pow(j+1, i);
            }
        }
        
        return ambiguityRemover.disambiguate(p);
    }

    /** Estimate the Mth order parameter */
    protected double estimateM(Double[] y, int M){

        //if it's the last parameter just run sample circular mean
        //otherwise the estimator needs the phase parameter bounded,
        //which is really silly.
        if(M == 0){
            SampleCircularMean vm = new SampleCircularMean();
            return vm.estimateBearing(y);
        }

        double[] d = VectorFunctions.mthDifference(y, M);

        //compute the denominator
        double dprod = 1.0;
        for(int i = 0; i <= 2*M; i++){
            dprod *= d.length+i;
        }
        double Mfac = Util.factorial(M);
        double fprod = Util.factorial(2*M + 1)/( Mfac*Mfac*dprod );

        double aM = 0.0;
        for(int j = 0; j < d.length; j++){

            double prod = 1.0;
            for(int i = 1; i <= M; i++)
                prod *= (j + i)*(d.length - j - 1 + i);

            aM += Util.centeredFracPart(d[j])*prod;

        }
        return fprod * aM / Mfac;
    }

}
