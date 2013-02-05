/*
 * GlueAnstarEstimatorTest.java
 * JUnit based test
 *
 * Created on 16 August 2007, 15:31
 */

package pubsim.fes;

import pubsim.fes.GlueAnstarEstimator;
import pubsim.fes.NoisyComplexSinusoid;
import junit.framework.*;

/**
 *
 * @author Robby
 */
public class GlueAnstarEstimatorTest extends TestCase {
    
    public GlueAnstarEstimatorTest(String testName) {
        super(testName);
    }

    /**
     * Test of estimateFreq method, of class simulator.fes.GlueAnstarEstimator.
     */
    public void testEstimateFreq() {
        System.out.println("estimateFreq");
        
        int iters = 100;
        double f = 0.2;
        int N = 21;
        GlueAnstarEstimator instance = new GlueAnstarEstimator(N);
        
        NoisyComplexSinusoid signal = new NoisyComplexSinusoid(N);
        signal.setFrequency(f);
        signal.setPhase(0.3);
        pubsim.distributions.GaussianNoise noise = new pubsim.distributions.GaussianNoise(0.0,0.001);
        signal.setNoiseGenerator(noise);
        
        signal.generateReceivedSignal();
        double result = instance.estimateFreq(signal.getReal(), signal.getImag());
        System.out.println("f = " + result);
        assertEquals(true, Math.abs(result - f)<0.001);
        
    }
    
}
