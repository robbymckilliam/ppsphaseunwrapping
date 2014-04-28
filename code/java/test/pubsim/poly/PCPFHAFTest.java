/*
 * @author Robby McKilliam
 */
package pubsim.poly;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import pubsim.VectorFunctions;
import pubsim.distributions.GaussianNoise;

/**
 *
 * @author Robby McKilliam
 */
public class PCPFHAFTest {
    
    public PCPFHAFTest() {
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
      
        int[][] tau = { {3}, {5} };
        PCPFHAF inst = new PCPFHAF(m,n,tau);

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
        System.out.println("testEstimate order 5");
        double tol = 1e-4;  

        int n = 49;
        double[] params = {0.11, 0.05002, 0.0205, 1e-5, 1e-6, 1e-7};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));

        siggen.generateReceivedSignal();
      
        int[][] tau = { {3,5}, {4,4} };
        PCPFHAF inst = new PCPFHAF(m,n,tau);

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
    public void testEstimate5withNequal199() {
        System.out.println("testEstimate order 5 with N = 199");
        double tol = 1e-4;  

        int n = 199;
        double[] params = {0.25, 0.25, 6.0e-4, 1.0e-6, 1.32e-9, 1.32e-12};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));

        siggen.generateReceivedSignal();
      
        int[][] tau = { {20,20}, {15, 25}, {17, 23}, {13,29} };
        PCPFHAF inst = new PCPFHAF(m,n,tau);

        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());

        System.out.println(VectorFunctions.print(params));
        System.out.println(VectorFunctions.print(p));

        assertEquals(params[m], p[m],tol);
        assertEquals(params[m-1], p[m-1],tol);
    }

}