package pubsim.poly;

import flanagan.complex.Complex;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import pubsim.distributions.GaussianNoise;

/**
 *
 * @author harprobey
 */
public class PHAFTest {
    
    public PHAFTest() {
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
     * Test of estimate method, of class PHAF.
     */
    @Test
    public void testSameAsHAFWithOneLagSet() {
        double tol = 1e-4;
        System.out.println("same as HAF with one lag");
        int n = 24;
        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));
        siggen.generateReceivedSignal();
        
        HAF tester = new HAF(m,n);
        int[][] tau = new int[1][m-2];
        tau[0] = tester.gettau();
        PHAF inst = new PHAF(m, n, tau, 10);
        
        double testest = tester.estimate(siggen.getReal(), siggen.getImag())[m];
        double est = inst.estimate(siggen.getReal(), siggen.getImag())[m];
        
        System.out.println(est + ", " + testest);
        
        assertEquals(testest, est, tol);   
    }
    
    /**
     * Test of estimate method, of class PHAF.
     */
    @Test
    public void testEstimate2HAFS() {
        double tol = 1e-4;
        System.out.println("estimate with two lags");
        int n = 24;
        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));
        siggen.generateReceivedSignal();
        
        int[][] tau = new int[2][m-2];
        tau[0] = new int[] {n/m, n/m};
        tau[1] = new int[] {n/m-2, n/m+2};
        PHAF inst = new PHAF(m, n, tau, 10);
        
        double est = inst.estimate(siggen.getReal(), siggen.getImag())[m];
        
        System.out.println(est);
        
        assertEquals(params[m], est, tol);   
    }
    
        /**
     * Test of estimate method, of class PHAF.
     */
    @Test
    public void testEstimate3HAFS() {
        double tol = 1e-4;
        System.out.println("estimate with three lags");
        int n = 24;
        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
        int m = params.length-1;

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0.000001));
        siggen.generateReceivedSignal();
        
        int[][] tau = new int[3][m-2];
        tau[0] = new int[] {n/m, n/m};
        tau[1] = new int[] {n/m-2, n/m+2};
        tau[2] = new int[] {n/m-4, n/m+4};
        PHAF inst = new PHAF(m, n, tau, 10);
        
        double est = inst.estimate(siggen.getReal(), siggen.getImag())[m];
        
        System.out.println(est);
        
        assertEquals(params[m], est, tol);   
    }

}
