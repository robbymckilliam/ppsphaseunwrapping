/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly.dptspecialm2;

import Jama.Matrix;
import pubsim.Complex;
import pubsim.VectorFunctions;
import pubsim.optimisation.NewtonRaphson;

/**
 * Polynomial phase estimator that for m=2 that only searches the DPT region.
 * O(N^2logN) complexity.
 * @author Robby McKilliam
 */
public class MaximumLikelihood extends pubsim.poly.MaximumLikelihood {

    Complex[] z;
    double[] realp, imagp;

    public MaximumLikelihood(int n) {
        super(2,n);
    }

    
    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(N != real.length) throw new RuntimeException("Data length does not equal " + N);
        
        PolynomialPhaseLikelihood func
                = new PolynomialPhaseLikelihood(real, imag);
        NewtonRaphson newtonRaphson
                = new NewtonRaphson(func);

        int numsamples = 4*N;

        double p2start = -(2.0/N)/4.0;
        double p2end = (2.0/N)/4.0;
        double step = (p2end - p2start)/numsamples;

        Matrix p = null;
        double D = Double.NEGATIVE_INFINITY;

        for(double p2 = p2start; p2 < p2end; p2 += step){
            for (int j = 0; j < real.length; j++)
                z[j] = new Complex(real[j], imag[j]);

            for (int i = 0; i < N; i++) {
                double cs = Math.cos(-2.0 * Math.PI * p2 * (i+1)*(i+1) );
                double ss = Math.sin(-2.0 * Math.PI * p2 * (i+1)*(i+1) );
                z[i] = z[i].times(new Complex(cs, ss));
            }

            for (int j = 0; j < real.length; j++){
                realp[j] = z[j].re();
                imagp[j] = z[j].im();
            }

            //estimate the supposed frequency using the periodogram estimator
            double f = freqest.estimateFreq(realp, imagp);
            //System.out.println(f);

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
            Matrix params = new Matrix(getOrder()+1, 1);
            //System.out.println("params before = " + VectorFunctions.print(params));
            params.set(0,0, pha);
            params.set(1,0,f);
            params.set(2,0,p2);

            //test it, if it's good save it.
            double dist = func.value(params);
            //System.out.println("params = " + VectorFunctions.print(params));
            //System.out.println("dist = " + dist);
            if(dist > D){
                D = dist;
                p = params.copy();
                //System.out.println("p = " + VectorFunctions.print(p));
            }

        }

        //refine the best parameter using Newton's method
        try{
            p = newtonRaphson.maximise(p);
        }catch(Exception e){
            throw new ArithmeticException(e.getMessage());
        }
        //System.out.println(D);
        //System.out.println(VectorFunctions.print(p));
        return ambiguityRemover.disambiguate(VectorFunctions.unpackRowise(p));

    }

        
}
