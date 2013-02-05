/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import pubsim.poly.PolynomialPhaseSignal;
import pubsim.distributions.GaussianNoise;
import junit.framework.TestCase;
import pubsim.VectorFunctions;

/**
 *
 * @author Robby
 */
public class PolynomialPhaseSignalTest extends TestCase {
    
    public PolynomialPhaseSignalTest(String testName) {
        super(testName);
    }            

    /**
     * Test of generateReceivedSignal method, of class PolynomialPhaseSignal.
     */
    public void testGenerateReceivedSignal() {
        System.out.println("generateReceivedSignal");
        
        int n = 4;
        double[] params = {1, 0.1, 0.2};
        
        PolynomialPhaseSignal instance = new PolynomialPhaseSignal(n);
        instance.setParameters(params);
        instance.setNoiseGenerator(new GaussianNoise(0, 0.0));
        
        //these are taken from matlab
        double[] expreal = {-0.3090, 1.0000, 0.8090, -0.8090};
        double[] expimag = {0.9511, -0.0000, 0.5878, -0.5878};
        
        instance.generateReceivedSignal();
        
        System.out.println("real = " + VectorFunctions.print(instance.getReal()));
        System.out.println("imag = " + VectorFunctions.print(instance.getImag()));
        
        assertEquals(VectorFunctions.distance_between(expreal, instance.getReal())<0.001, true);
        assertEquals(VectorFunctions.distance_between(expimag, instance.getImag())<0.001, true);
        
    }

        /**
     * Test of generateReceivedSignal method, of class PolynomialPhaseSignal.
     */
    public void testGenerateRandomParameters() {
        System.out.println("testGenerateRandomParameters");
        int m = 2;
        PolynomialPhaseSignal.RandomParameterGenerator pgen
                = new PolynomialPhaseSignal.RandomParameterGenerator(m);

        double[] p = pgen.generateParameters();
        assertEquals(3, p.length);
        System.out.println(VectorFunctions.print(p));

    }

}
