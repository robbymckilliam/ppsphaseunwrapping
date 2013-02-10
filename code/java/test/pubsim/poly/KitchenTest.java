/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import pubsim.poly.Kitchen;
import pubsim.poly.PolynomialPhaseSignal;
import pubsim.distributions.GaussianNoise;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import pubsim.VectorFunctions;
import static org.junit.Assert.*;

/**
 *
 * @author Robby
 */
public class KitchenTest {

    public KitchenTest() {
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
     * Test of estimate method, of class Kitchen.
     */
    @Test
    public void testEstimate() {
        System.out.println("testEstimate");

        int n = 20;
        double[] params = {0.1, 0.2, 0.05};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));

        siggen.generateReceivedSignal();

        Kitchen inst = new Kitchen(m,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(p));

        assertTrue(VectorFunctions.distance_between(p, params) < 0.001);
    }

    /**
     * Test of estimateM method, of class Kitchen.
     */
    @Test
    public void testHighestOrderParameter() {
        System.out.println("testHighestOrderParameter");

        int n = 20;
        double[] params = {0.1, 0.2, 0.05};
        int a = params.length;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));

        siggen.generateReceivedSignal();

        Kitchen inst = new Kitchen(a,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(p));

        assertTrue(Math.abs(p[a-1] - params[a-1])< 0.0001);

    }

}