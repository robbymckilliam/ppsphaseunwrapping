/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pubsim.poly;

import org.junit.*;
import static org.junit.Assert.*;
import pubsim.VectorFunctions;
import pubsim.distributions.circular.WrappedUniform;

/**
 *
 * @author mckillrg
 */
public class CircularNoisePolynomialPhaseSignalTest {

    public CircularNoisePolynomialPhaseSignalTest() {
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
     * Test of generateReceivedSignal method, of class
     * CircularNoisePolynomialPhaseSignal.
     */
    @Test
    public void testGenerateReceivedSignal() {
        System.out.println("generateReceivedSignal");
        int N = 5;
        CircularNoisePolynomialPhaseSignal instance = new CircularNoisePolynomialPhaseSignal(N);

        WrappedUniform noise = new WrappedUniform(0, 0.01);
        instance.setNoiseGenerator(noise);

        double[] p = {0.0, 0.0, 0.0};
        instance.setParameters(p);

        instance.generateReceivedSignal();
        System.out.println(VectorFunctions.print(instance.getReal()));
        System.out.println(VectorFunctions.print(instance.getImag()));
        
        instance.generateReceivedSignal();
        System.out.println(VectorFunctions.print(instance.getReal()));
        System.out.println(VectorFunctions.print(instance.getImag()));
    }
}
