/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.fes;

import flanagan.complex.Complex;
import flanagan.math.FourierTransform;


/**
 * Implementation of the Quinn and Fernandes frequency estimator
 *
 * @author Robby McKilliam
 */
public class QuinnFernades implements FrequencyEstimator{

    protected final int oversampled,  N;
    protected final FourierTransform fft;
    protected final Complex[] sig;

    /**eta term in Barry's paper.  This is related to the ARMA(1,1) model.*/
    protected Complex[] eta;

    /**Max number of iterations for the Newton step */
    protected static final int MAX_ITER = 15;

    /**Step variable for the Newton step */
    protected static final double EPSILON = 1e-10;

    /**
     * By default, oversample = 1. Quinn Fernandes does not need
     * zero padding.
     */
    public QuinnFernades(int N){
        oversampled = 1;
        this.N = N;
        sig = new Complex[FourierTransform.nextPowerOfTwo(oversampled * N)];
        fft = new FourierTransform();
        eta = new Complex[FourierTransform.nextPowerOfTwo(oversampled * N) + 1];
    }

    /**
     * Constructor that sets the number of samples to be taken of
     * the periodogram. You don't really need to oversample though.
     */
    public QuinnFernades(int N, int oversampled) {
        this.oversampled = oversampled;
        this.N = N;
        sig = new Complex[FourierTransform.nextPowerOfTwo(oversampled * N)];
        fft = new FourierTransform();
        eta = new Complex[FourierTransform.nextPowerOfTwo(oversampled * N) + 1];
    }

    @Override
    public final double estimateFreq(double[] real, double[] imag) {
        //construct zero padded complex signal
        for (int i = 0; i < N; i++) {
            sig[i] = new Complex(real[i], imag[i]);
        }
        for (int i = N; i < sig.length; i++) {
            sig[i] = new Complex(0.0, 0.0);
        }

        for (int i = N; i < sig.length; i++) {
            eta[i] = new Complex(0.0, 0.0);
        }

        fft.setData(sig);
        fft.transform();
        Complex[] ft = fft.getTransformedDataAsComplex();

        //FFT output is backwards, flip it
        double maxp = 0;
        double fhat = 0.0;
        double f = 0.0;
        double fstep = 1.0 / ft.length;
        for (int i = 0; i < ft.length; i++) {
            double p = ft[i].squareAbs();
            if (p > maxp) {
                maxp = p;
                fhat = f;
            }
            f-=fstep;
        }

        fhat -= Math.round(fhat);

        //Now implement QuinnFernandes iterations.
        int numIter = 0;
        double lastf = fhat - 2 * EPSILON;
        while(Math.abs(fhat - lastf) > EPSILON && numIter <= MAX_ITER){
            //System.out.println(fhat);
            lastf = fhat;
            Complex efhat = new Complex(Math.cos(2*Math.PI*fhat),
                                        Math.sin(2*Math.PI*fhat));
            //compute next eta
            eta[0] = new Complex(0,0);
            for(int t = 0; t < sig.length; t++){
                eta[t+1] = sig[t].plus( efhat.times(eta[t]) );
            }
            //compute next fhat
            Complex efhatc = efhat.conjugate();
            double num = 0.0, den = 0.0;
            for(int t = 0; t < sig.length; t++){
                num += efhatc.times(sig[t].times(eta[t].conjugate())).getImag();
                den += eta[t].squareAbs();
            }
            fhat += 2*num/den/2/Math.PI;

            numIter++;
        }

        return fhat - Math.round(fhat);
    }


}
