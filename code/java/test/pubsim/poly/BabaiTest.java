/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import pubsim.poly.PolynomialPhaseEstimatorInterface;
import pubsim.poly.BabaiTest;
import pubsim.poly.PolynomialPhaseSignal;
import Jama.Matrix;
import pubsim.distributions.GaussianNoise;
import pubsim.lattices.Vnmstar.VnmStar;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import pubsim.VectorFunctions;
import static org.junit.Assert.*;

/**
 *
 * @author Robby McKilliam
 */
public class BabaiTest {

    public BabaiTest() {
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
     * Test of estimate method, of class DPTEstimator.
     */
    @Test
    public void testEstimate() {
        System.out.println("testEstimate");

        int n = 24;
        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
        int m = params.length-1;

        Matrix M = VnmStar.getGeneratorMatrix(m, n-m-1);
        System.out.println(VectorFunctions.print(M));

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.00001));

        siggen.generateReceivedSignal();
        PolynomialPhaseEstimatorInterface inst = new Babai(m,n);
        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(p));

        assertTrue(VectorFunctions.distance_between(p, params) < 0.001);

    }

}