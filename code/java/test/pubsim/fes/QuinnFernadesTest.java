/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.fes;

import static org.junit.Assert.assertEquals;
import org.junit.*;

/**
 *
 * @author Robby McKilliam
 */
public class QuinnFernadesTest {

    public QuinnFernadesTest() {
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

    /**
     * Test of estimateFreq method, of class QuinnFernades.
     */
    @Test
    public void testEstimateFreq() {
        System.out.println("estimateFreq");

        int iters = 10;
        double f = -0.41111;
        QuinnFernades instance = new QuinnFernades(64);

        NoisyComplexSinusoid signal = new NoisyComplexSinusoid(64);
        signal.setFrequency(f);
        pubsim.distributions.GaussianNoise noise = new pubsim.distributions.GaussianNoise(0.0,0.0001);
        signal.setNoiseGenerator(noise);

        for(int i=0; i < iters; i++ ){
            signal.generateReceivedSignal();
            double result = instance.estimateFreq(signal.getReal(), signal.getImag());
            System.out.println(result);
            assertEquals(result, f, 0.001);
        }

    }

}