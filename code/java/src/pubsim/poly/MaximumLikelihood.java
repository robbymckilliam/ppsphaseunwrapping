/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import Jama.Matrix;
import pubsim.lattices.util.PointInParallelepiped;
import pubsim.optimisation.AutoDerivativeFunction;
import pubsim.optimisation.FunctionAndDerivatives;
import pubsim.optimisation.NewtonRaphson;
import pubsim.Complex;
import pubsim.VectorFunctions;
import pubsim.fes.PeriodogramFFTEstimator;

/**
 * Implements a (approximate) maximum likelihood estimator for
 * polynomial phase signals.  This samples the identifiable region
 * performing Newton's method a whole lot of times.
 * @author Robby McKilliam
 */
public class MaximumLikelihood extends AbstractPolynomialPhaseEstimator{

    final protected int N;
    final protected PeriodogramFFTEstimator freqest;

    int samples[];
    final double[] realp, imagp;
    final Complex[] z;
    
    /**
     * @param m : polynomial order
     * @param samples : number of samples used per parameter in ML search.
     * Default samples = 100
     */
    public MaximumLikelihood(int m, int n){
        super(m);
        this.samples = new int[m-1];
        N = n;
        realp = new double[N];
        imagp = new double[N];
        z = new Complex[N];
        freqest = new PeriodogramFFTEstimator(N);
        for(int i = 0; i < samples.length; i++)
            samples[i] = (int)Math.round(4*Math.pow( N, i+2 ));
    }
    
    /**
     * @param samples : number of samples used per parameter in ML search.
     */
    public void setSamples(int[] samples){
        this.samples = samples;
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(N != real.length) throw new RuntimeException("Data length does not equal " + N);

        PolynomialPhaseLikelihood func
                = new PolynomialPhaseLikelihood(real, imag);
        NewtonRaphson newtonRaphson
                = new NewtonRaphson(func);

        PointInParallelepiped points
                = new PointInParallelepiped(
                    ambiguityRemover.getBasisMatrix().getMatrix(2, m, 2, m),
                                            samples);

        Matrix p = null;
        double D = Double.NEGATIVE_INFINITY;
        while(points.hasMoreElements()){
            //Matrix pt = newtonRaphson.maximise(points.nextElement());
            Matrix pt = points.nextElement();


            //compute the supposed frequency signal by remove higher order
            //parameters
            for (int j = 0; j < real.length; j++){
                z[j] = new Complex(real[j], imag[j]);
            }
            for (int j = 0; j < real.length; j++) {
                 for (int i = 0; i < m-1; i++) {
                    double cs = Math.cos(-2.0 * Math.PI * pt.get(i, 0) * Math.pow(j + 1, i+2));
                    double ss = Math.sin(-2.0 * Math.PI * pt.get(i, 0) * Math.pow(j + 1, i+2));
                    z[j] = z[j].times(new Complex(cs, ss));
                }
            }
            for (int j = 0; j < real.length; j++){
                realp[j] = z[j].re();
                imagp[j] = z[j].im();
            }

            //estimate the supposed frequency using the periodogram estimator
            double f = freqest.estimateFreq(realp, imagp);

            //now estimator the phase, just compute the complex mean
            for (int j = 0; j < real.length; j++) {
                double cs = Math.cos(-2.0 * Math.PI * f * (j+1));
                double ss = Math.sin(-2.0 * Math.PI * f * (j+1));
                z[j] = z[j].times(new Complex(cs, ss));
            }
            Complex c = new Complex(0,0);
            for (int j = 0; j < real.length; j++) c = c.add(z[j]);
            double pha = c.phase()/(2*Math.PI);
            //System.out.println(pha);

            //set the parameters vector
            Matrix params = new Matrix(m+1, 1);
            //System.out.println("params before = " + VectorFunctions.print(params));
            params.set(0,0, pha);
            params.set(1,0,f);
            for(int i = 2; i < m+1; i++) params.set(i,0, pt.get(i-2,0));

            //test it, if it's good save it.
            double dist = func.value(params);
            if(dist > D){
                D = dist;
                p = params.copy();
            }
        }
        //double[] parray = {0.1,0.1};
        //double dist = func.value(VectorFunctions.columnMatrix(parray));
        //p = VectorFunctions.columnMatrix(parray);

        //refine the best parameter using Newton's method
        try{
            p = newtonRaphson.maximise(p);
        }catch(Exception e){
            throw new ArithmeticException(e.getMessage());
        }
        //System.out.println(dist);
        //System.out.println(D);
        //System.out.println(VectorFunctions.print(p));
        return ambiguityRemover.disambiguate(VectorFunctions.unpackRowise(p));
    }

    public static class PolynomialPhaseLikelihoodAutoDerivative
            extends AutoDerivativeFunction {

        double[] yr, yi;
        int N;

        protected PolynomialPhaseLikelihoodAutoDerivative(){}

         /**
         * @param yr real part of signal
          * @param yi complex part of signal
         */
        public PolynomialPhaseLikelihoodAutoDerivative(double[] yr, double[] yi){
            super(1e-6);
            this.yr = yr;
            this.yi = yi;
            N = yr.length;
        }

        /**
         * x is m column vector containing the polynomial phase
         * parameter.  x = [p0, p1, p2, ... ]
         * @return Value of the likelihood function for these parameters
         */
        @Override
        public double value(Matrix x) {
            int M = x.getRowDimension();
            double val = 0.0;
            for(int n = 0; n < N; n++){
                double phase = 0.0;
                for(int m = 0; m < M; m++){
                    double p = x.get(m, 0);
                    phase += p * Math.pow(n+1, m);
                }
                double real = Math.cos(2*Math.PI*phase);
                double imag = Math.sin(2*Math.PI*phase);
                val += (yr[n] - real)*(yr[n] - real);
                val += (yi[n] - imag)*(yi[n] - imag);
            }
            return -val;
        }
    }

    public static class PolynomialPhaseLikelihood
            implements FunctionAndDerivatives{

        double[] yr, yi;
        double[] ymag, yphase;
        double[] spdiff, cpdiff;
        int N;

         /**
         * @param yr real part of signal
          * @param yi complex part of signal
         */
        public PolynomialPhaseLikelihood(double[] yr, double[] yi){
            N = yr.length;
            this.yr = yr;
            this.yi = yi;
            ymag = new double[N];
            yphase = new double[N];
            spdiff = new double[N];
            cpdiff = new double[N];
            for(int n = 0; n < N; n++){
                ymag[n] = Math.sqrt(yr[n]*yr[n] + yi[n]*yi[n]);
                yphase[n] = Math.atan2(yi[n], yr[n]);
            }
            N = yr.length;
        }

        @Override
        public double value(Matrix x) {
            int M = x.getRowDimension();
            double val = 0.0;
            for(int n = 0; n < N; n++){
                double phase = 0.0;
                for(int m = 0; m < M; m++){
                    double p = x.get(m, 0);
                    phase += p * Math.pow(n+1, m);
                }
                double real = Math.cos(2*Math.PI*phase);
                double imag = Math.sin(2*Math.PI*phase);
                val += (yr[n] - real)*(yr[n] - real);
                val += (yi[n] - imag)*(yi[n] - imag);
            }
            return -val;
        }

        @Override
        public Matrix hessian(Matrix x) {
            int M = x.getRowDimension();
            //precompute required sin values
            for(int n = 0; n < N; n++){
                double phase = 0.0;
                for(int m = 0; m < M; m++){
                    phase += x.get(m, 0)*Math.pow((n+1), m);
                }
                cpdiff[n] = Math.cos(2*Math.PI * phase - yphase[n]);
            }
            //compute Hessian elements
            Matrix H = new Matrix(M, M);
            for(int m = 0; m < M; m++){
                for(int k = 0; k < M; k++){
                    double grad2 = 0.0;
                    for(int n = 0; n < N; n++){
                        grad2 += ymag[n] * Math.pow((n+1), m) *
                                Math.pow((n+1), k) * cpdiff[n];
                    }
                    grad2 *= -8*Math.PI*Math.PI;
                    H.set(m,k, grad2);
                }
            }
            return H;
        }

        @Override
        public Matrix gradient(Matrix x) {
            int M = x.getRowDimension();
            //precompute required sin values
            for(int n = 0; n < N; n++){
                double phase = 0.0;
                for(int m = 0; m < M; m++){
                    phase += x.get(m, 0)*Math.pow((n+1), m);
                }
                spdiff[n] = Math.sin(2*Math.PI * phase - yphase[n]);
            }
            //compute gradients
            Matrix g = new Matrix(M, 1);
            for(int m = 0; m < M; m++){
                double grad = 0.0;
                for(int n = 0; n < N; n++){
                    grad += ymag[n] * Math.pow((n+1), m) * spdiff[n];
                }
                grad *= -4*Math.PI;
                g.set(m,0, grad);
            }
            return g;
        }
        
    }

}
