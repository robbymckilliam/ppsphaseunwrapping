/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pubsim.poly;

import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import pubsim.Complex;
import static pubsim.VectorFunctions.print;
import pubsim.distributions.GaussianNoise;

/**
 *
 * @author Robby McKilliam
 */
public class CubicPhaseFunctionTest {
    
    public CubicPhaseFunctionTest() {
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
     * Test of z method, of class CubicPhaseFunction.
     */
    @Test
    public void testZ() {
        System.out.println("test z");
        double tol = 1e-8;
        int n = 3;
        CubicPhaseFunction instance = new CubicPhaseFunction(n);
        double[] real = {1.0,2.0,3.0};
        double[] imag = {1.0,2.0,3.0};
        instance.estimate(real, imag);
        assertTrue( (instance.z(-2) - new Complex(0,0)).abs() < tol );
        assertTrue( (instance.z(-1) - new Complex(1,1)).abs() < tol );
        assertTrue( (instance.z(0) - new Complex(2,2)).abs() < tol );
        assertTrue( (instance.z(1) - new Complex(3,3)).abs() < tol );
        assertTrue( (instance.z(2) - new Complex(0,0)).abs() < tol );
    }
    
     /**
     * Test of estimate method, of class DPTEstimator.
     */
    @Test
    public void testCPmax() {
        System.out.println("testHighestOrderParameter");
        
        double tol = 1e-5;
        
        int N = 257;
        int m = 3;
        double[] oparams = {1.0, Math.PI/8, 0.005, 0.00001};
        CubicPhaseFunction inst = new CubicPhaseFunction(N,10*N);
        double[] params = inst.transformToStandardBasis(oparams);
        
        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(N);
        siggen.setParameters(params);
        siggen.setNoiseGenerator(new GaussianNoise(0, 0));
        siggen.generateReceivedSignal();
        inst.estimate(siggen.getReal(), siggen.getImag());
        
        //check that maxCP(0) is twice the quadratic parameter
        assertEquals(inst.maxCP(0)/2, oparams[2], tol); 

    }
    
//    /**
//     * Test estimate of quadratic and cubic parameters
//     */
//    @Test
//    public void testQuadraticAndCubic() {
//        System.out.println("test quadratic and cubic");
//        
//        double tol = 1e-5;
//        
//        int n = 515;
//        int m = 3;
//        double[] params = {1.0/2/Math.PI, Math.PI/8/2/Math.PI, 0.005/2/Math.PI, 0.00001/2/Math.PI};
//        int a = params.length;
//
//        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n);
//        siggen.setParameters(params);
//        siggen.setNoiseGenerator(new GaussianNoise(0, 0.0));
//
//        siggen.generateReceivedSignal();
//
//        CubicPhaseFunction inst = new CubicPhaseFunction(n);
//
//        double[] p = inst.estimate(siggen.getReal(), siggen.getImag());
//
//        System.out.println(print(p));
//
//        assertTrue(Math.abs(p[m] - params[m]) < tol);
//        assertTrue(Math.abs(p[m-1] - params[m-1]) < tol);
//
//    }

}