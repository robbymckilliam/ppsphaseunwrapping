/*
 */

package pubsim.poly;

import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import pubsim.Complex;
import pubsim.VectorFunctions;
import pubsim.distributions.GaussianNoise;

/**
 *
 * @author Robby McKilliam
 */
public class HAFTest {

    public HAFTest() {
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
     * Test of PPT method, of class HAF.
     */
//    @Test
//    public void PPTHasFirstMElementsZero() {
//        int m = 4;
//        int n = 10;
//        Complex[] y = VectorFunctions.randomComplex(n);
//        HAF instance = new HAF(m);
//        instance.setSize(n);
//        Complex[] result = instance.PPT(m, y);
//        //System.out.print(VectorFunctions.print(result));
//        for(int i = 0; i <= 2; i++){
//            assertTrue(result[i].re() == 0.0);
//            assertTrue(result[i].im() == 0.0);
//        }
//    }

    /**
     * Test of PPT2 method, of class HAF.
     */
    @Test
    public void PPT2HasLastElementCorrect() {
        System.out.println("PPT2");
        int m = 4;
        int n = 10;
        Complex[] y = VectorFunctions.randomComplex(n);
        HAF instance = new HAF(m,n);
        int tau = instance.gettau()[0];
        flanagan.complex.Complex[] result = instance.PPT2(VectorFunctions.simComplexArrayToFlanComplexArray(y),tau);
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
     * Test of estimate method, of class HAF.
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

        HAF inst = new HAF(params.length,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(p[a-1]);

        assertTrue(Math.abs(p[a-1] - params[a-1]) < 0.0001);

    }

    /**
     * Test of estimate method, of class HAF.
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

        HAF inst = new HAF(m,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(p));

        assertTrue(VectorFunctions.distance_between(p, params) < 0.001);

    }
    
        /**
     * Test of estimate method, of class HAF.
     */
    @Test
    public void testEstimateMultiLag() {
        System.out.println("testEstimateMultiLag");

        int n = 24;
        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.00001));

        siggen.generateReceivedSignal();

        HAF inst = new HAF(m,n, new int[] {n/m+1, n/m-1});

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(p));

        assertTrue(VectorFunctions.distance_between(p, params) < 0.001);

    }

    /**
     * Test of estimate method, of class HAF.
     */
    @Test
    public void testEstimate5() {
        System.out.println("testEstimate order 5");

        int n = 34;
        double[] params = {0.11, 0.05002, 0.0205, 1e-4, 1e-6, 1e-7};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.0));

        siggen.generateReceivedSignal();

        HAF inst = new HAF(m,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(params));
        System.out.println(VectorFunctions.print(p));

        assertTrue(VectorFunctions.distance_between(p, params) < 0.001);

    }
    
    @Test
    public void testCalculcateObjective() {
        double tol = 1e-7;
        int n = 24;
        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.00001));

        siggen.generateReceivedSignal();
        flanagan.complex.Complex[] z = new flanagan.complex.Complex[n];
        for (int i = 0; i < n; i++) z[i] = new flanagan.complex.Complex(siggen.getReal()[i], siggen.getImag()[i]);
        
        HAF inst = new HAF(m,n);
       
        /** test highest order parameter m */
        flanagan.complex.Complex[] fft = inst.FFTHAF(z,m);
        double f = 0.0;
        double fstep = 1.0 / fft.length;
        for(int i = 0; i < fft.length; i++){
            //System.out.println(fft[i].squareAbs() + ", " + inst.calculateObjective(f, z, m));
            assertEquals(fft[i].squareAbs(), inst.calculateObjective(f, z, m), tol);
            f-=fstep;
        }
        
        /** test parameter m-1 */
        fft = inst.FFTHAF(z,m-1);
        f = 0.0;
        fstep = 1.0 / fft.length;
        for(int i = 0; i < fft.length; i++){
            //System.out.println(fft[i].squareAbs() + ", " + inst.calculateObjective(f, z, m));
            assertEquals(fft[i].squareAbs(), inst.calculateObjective(f, z, m-1), tol);
            f-=fstep;
        }
        
    }
    

}