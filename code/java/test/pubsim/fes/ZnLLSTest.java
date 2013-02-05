/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pubsim.fes;

import org.junit.*;
import static org.junit.Assert.*;

/**
 *
 * @author mckillrg
 */
public class ZnLLSTest {
    
    public ZnLLSTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testSomeMethod() {
        System.out.println("test ZnLLS frequency estimator");
        
        int iters = 5;
        double f = 0.2;
        int N = 20;
        ZnLLS instance = new ZnLLS(N);
        
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
