/*
 * NoisyComplexSinusoid.java
 *
 * Created on 9 August 2007, 12:00
 */

package pubsim.fes;

import pubsim.SignalGenerator;
import pubsim.distributions.NoiseGenerator;

/**
 * Generator for a noisy sinusoid.
 * @author Robby McKilliam
 */
public class NoisyComplexSinusoid implements SignalGenerator<Double>{
    
    protected double f;
    protected double p;
    protected int N;
    
    protected double[] real;
    protected double[] imag;
    
    protected NoiseGenerator<Double> noise;

    /** default constructor: f = 0.1, rate = 1, N = 10, p = 0.0 */
    protected NoisyComplexSinusoid() {    }

    public NoisyComplexSinusoid(int N){
        setSize(N);
    }

    /** Set the number of samples */
    private void setSize(int N){
        f = 0.1;
        p = 0.0;
        this.N = N;
        real = new double[N];
        imag = new double[N];
    }

    
    /** Set the noise type for the signal */
    @Override
    public void setNoiseGenerator(NoiseGenerator noise){
        this.noise = noise;
    }
    
    @Override
    public NoiseGenerator getNoiseGenerator(){ return noise; }
    
    /** Set the frequency of the complex sinusoid in Hz */
    public void setFrequency(double f){
        this.f = f;
    }
    public double getFrequency(){ return f; }
    
   /** Set the phase of the sinusoid */
    public void setPhase(double p){
        this.p = p;
    }
    
    /** Generates the real and imaginary components and returns null.
     * Use getReal and getImag to get the individual components.
     */
    @Override
    public Double[] generateReceivedSignal(){
        for(int i = 0; i < N; i++){
            real[i] = Math.cos(i*2*Math.PI*f + p) + noise.getNoise();
            imag[i] = Math.sin(i*2*Math.PI*f + p) + noise.getNoise();
        }
        return null;
    }
    
    /** Return the noisy real component of the signal */
    public double[] getReal() { return real; }
    
    /** Return the noisy imaginary component of the signal */
    public double[] getImag() { return imag; }

    @Override
    public int getLength() {
        return N;
    }
    
}
