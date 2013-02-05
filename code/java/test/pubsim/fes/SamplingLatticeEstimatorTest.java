/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.fes;

import junit.framework.TestCase;

/**
 *
 * @author Robby McKilliam
 */
public class SamplingLatticeEstimatorTest extends TestCase {
    
    public SamplingLatticeEstimatorTest(String testName) {
        super(testName);
    }            
    /**
     * Test of estimateFreq method, of class SamplingLatticeEstimator.
     */
    public void testEstimateFreq() {
        System.out.println("estimateFreq");
        
        int iters = 5;
        double f = 0.2;
        int N = 20;
        SamplingLatticeEstimator instance = new SamplingLatticeEstimator(N,100);
        
        NoisyComplexSinusoid signal = new NoisyComplexSinusoid(N);
        signal.setFrequency(f);
        pubsim.distributions.GaussianNoise noise = new pubsim.distributions.GaussianNoise(0.0,0.00001);
        signal.setNoiseGenerator(noise);
        
        signal.generateReceivedSignal();
        double result = instance.estimateFreq(signal.getReal(), signal.getImag());
        System.out.println("f = " + result);
        assertEquals(true, Math.abs(result - f)<0.001);
    }

}
