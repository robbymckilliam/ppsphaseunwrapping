/*
 * PeriodogramEstimatorTest.java
 * JUnit based test
 *
 * Created on 15 August 2007, 16:05
 */

package pubsim.fes;

import pubsim.fes.PeriodogramEstimator;
import pubsim.fes.NoisyComplexSinusoid;
import junit.framework.*;

/**
 *
 * @author Robby
 */
public class PeriodogramEstimatorTest extends TestCase {
    
    public PeriodogramEstimatorTest(String testName) {
        super(testName);
    }
    
    /**
     * Test of estimateFreq method, of class simulator.fes.PeriodogramEstimator.
     */
    public void testEstimateFreq() {
        System.out.println("estimateFreq");
        
        int iters = 100;
        double f = 0.1;
        PeriodogramEstimator instance = new PeriodogramEstimator(64);
        
        NoisyComplexSinusoid signal = new NoisyComplexSinusoid(64);
        signal.setFrequency(f);
        pubsim.distributions.GaussianNoise noise = new pubsim.distributions.GaussianNoise(0.0,0.001);
        signal.setNoiseGenerator(noise);
        
        for(int i=0; i < iters; i++ ){
            signal.generateReceivedSignal();
            double result = instance.estimateFreq(signal.getReal(), signal.getImag());
            assertEquals(true, Math.abs(result - f) < 0.02);
        }
        
    }
    
}
