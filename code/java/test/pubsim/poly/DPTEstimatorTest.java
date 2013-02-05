/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package pubsim.poly;

import pubsim.poly.DPTEstimator;
import pubsim.poly.PolynomialPhaseSignal;
import pubsim.distributions.GaussianNoise;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import pubsim.Complex;
import pubsim.VectorFunctions;

/**
 *
 * @author Robby McKilliam
 */
public class DPTEstimatorTest {

    public DPTEstimatorTest() {
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
     * Test of PPT method, of class DPTEstimator.
     */
//    @Test
//    public void PPTHasFirstMElementsZero() {
//        int m = 4;
//        int n = 10;
//        Complex[] y = VectorFunctions.randomComplex(n);
//        DPTEstimator instance = new DPTEstimator(m);
//        instance.setSize(n);
//        Complex[] result = instance.PPT(m, y);
//        //System.out.print(VectorFunctions.print(result));
//        for(int i = 0; i <= 2; i++){
//            assertTrue(result[i].re() == 0.0);
//            assertTrue(result[i].im() == 0.0);
//        }
//    }

    /**
     * Test of PPT2 method, of class DPTEstimator.
     */
    @Test
    public void PPT2HasLastElementCorrect() {
        System.out.println("PPT2");
        int m = 4;
        int n = 10;
        Complex[] y = VectorFunctions.randomComplex(n);
        DPTEstimator instance = new DPTEstimator(m,n);
        flanagan.complex.Complex[] result = instance.PPT2(VectorFunctions.simComplexArrayToFlanComplexArray(y));
        //System.out.print(VectorFunctions.print(result));

        //test first element is zero
        //assertTrue(result[0].re() == 0.0);
        //assertTrue(result[0].im() == 0.0);

        //test last element
        Complex last = y[n-1].times(y[n - 1 - (int)Math.round(n/((double)m-1))].conjugate());
        assertEquals(last.im(), result[result.length-1].getImag(), 0.0000001);
        assertEquals(last.re(), result[result.length-1].getReal(), 0.0000001);
    }

    /**
     * Test of estimate method, of class DPTEstimator.
     */
    @Test
    public void testHighestOrderParameter() {
        System.out.println("testHighestOrderParameter");
        
        int n = 60;
        double[] params = {0.3, 0.1, 0.002};
        int a = params.length;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.00001));

        siggen.generateReceivedSignal();

        DPTEstimator inst = new DPTEstimator(params.length,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(p[a-1]);

        assertTrue(Math.abs(p[a-1] - params[a-1]) < 0.0001);

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

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.00001));

        siggen.generateReceivedSignal();

        DPTEstimator inst = new DPTEstimator(m,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(p));

        assertTrue(VectorFunctions.distance_between(p, params) < 0.001);

    }


}