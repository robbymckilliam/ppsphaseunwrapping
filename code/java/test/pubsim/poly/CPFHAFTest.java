/*
 */
package pubsim.poly;

import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import pubsim.VectorFunctions;
import pubsim.distributions.GaussianNoise;

/**
 * @author Robby McKilliam
 */
public class CPFHAFTest {
    
    public CPFHAFTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

//    /**
//     * Test of estimate method, of class HAF.
//     */
//    @Test
//    public void testSameAsCPFforOrder3() {
//        System.out.println("testEstimate");
//        double tol = 1e-8;
//        
//        int n = 25;
//        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
//        int m = params.length-1;
//
//        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
//        siggen.setParameters(params);
//        siggen.setNoiseGenerator(new GaussianNoise(0, 0.00001));
//
//        siggen.generateReceivedSignal();
//
//        CPFHAF inst = new CPFHAF(m,n);
//        CPF test = new CPF(n);
//
//        double[] pinst = inst.estimate(siggen.getReal(), siggen.getImag());
//        double[] ptest = inst.estimate(siggen.getReal(), siggen.getImag());
//
//        assertEquals(pinst[3], ptest[3],tol);
//        assertEquals(pinst[2], ptest[2],tol);
//
//    }
    
    /**
     * Test of estimate method, of class HAF.
     */
    @Test
    public void testEstimate4() {
        System.out.println("testEstimate order 4");
        double tol = 1e-4;  

        int n = 45;
        double[] params = {0.11, 0.05002, 0.0205, 1e-5, 1e-6};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));

        siggen.generateReceivedSignal();

        CPFHAF inst = new CPFHAF(m,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(params));
        System.out.println(VectorFunctions.print(p));

        assertEquals(params[m], p[m],tol);
        assertEquals(params[m-1], p[m-1],tol);
    }
    
    /**
     * Test of estimate method, of class HAF.
     */
    @Test
    public void testEstimate5() {
        System.out.println("testEstimate order 4");
        double tol = 1e-4;  

        int n = 49;
        double[] params = {0.11, 0.05002, 0.0205, 1e-5, 1e-6, 1e-7};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));

        siggen.generateReceivedSignal();

        CPFHAF inst = new CPFHAF(m,n);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(params));
        System.out.println(VectorFunctions.print(p));

        assertEquals(params[m], p[m],tol);
        assertEquals(params[m-1], p[m-1],tol);
    }
    
        /**
     * Test of estimate method, of class HAF.
     */
    @Test
    public void testEstimate5withdistinctlags() {
        System.out.println("testEstimate order 5 with distinct lags");
        double tol = 1e-4;  

        int n = 49;
        double[] params = {0.11, 0.05002, 0.0205, 1e-5, 1e-6, 1e-7};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));

        siggen.generateReceivedSignal();

        int[] tau = {3, 7};
        CPFHAF inst = new CPFHAF(m,n,tau);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(params));
        System.out.println(VectorFunctions.print(p));

        assertEquals(params[m], p[m],tol);
        assertEquals(params[m-1], p[m-1],tol);
    }

}
