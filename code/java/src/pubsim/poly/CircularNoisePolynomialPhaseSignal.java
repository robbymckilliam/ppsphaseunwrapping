/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

/**
 * Like PolynommialPhaseSignal but the noise gets added to the phase.
 * @author Robby McKilliam
 */
public class CircularNoisePolynomialPhaseSignal extends PolynomialPhaseSignal{

    Double[] phasearray = new Double[0];

    public CircularNoisePolynomialPhaseSignal(int N){
        super(N);
    }
    
    /**
     * Returns the phase signal and places the real and imaginary parts
     * of the complex signal into real and imag arrays that can be retrieved
     * with getReal() and getImag()
     * @return
     */
    @Override
    public Double[] generateReceivedSignal() {
        if(phasearray.length != n) phasearray = new Double[n];

        for(int t = 0; t < n; t++){
            double phase = 0.0;
            for(int j = 0; j < params.length; j++)
                phase += Math.pow(t+1,j)*params[j];
            double pnoise = noise.getNoise();
            phasearray[t] = phase + pnoise;
            real[t] = Math.cos(2*Math.PI*(phase + pnoise));
            imag[t] = Math.sin(2*Math.PI*(phase + pnoise));
        }
        return phasearray;
    }

}
