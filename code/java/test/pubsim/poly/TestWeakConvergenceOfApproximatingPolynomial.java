/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pubsim.poly;

import org.junit.*;
import static pubsim.Util.fracpart;
import pubsim.distributions.circular.CircularRandomVariable;
import pubsim.distributions.circular.CircularUniform;

/**
 *
 * @author Robby McKilliam
 */
public class TestWeakConvergenceOfApproximatingPolynomial {
    
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

    protected double polyn(double x, double N){
        if(Math.abs(x) < 1.0/2.0 - 1.0/2.0/N ) return 2.0;
        double y = fracpart(x - 0.5);
        return (4 - 15*N)/2.0 + 4*3*10*N*N*N*y*y - 6*5*12*N*N*N*N*N*y*y*y*y;
    }
    
    @Test
    public void test() {
        System.out.println("testPolyFunctionValue");
        
        long fromN = 100000;
        int numruns = 30;

        
        CircularRandomVariable rv = new CircularUniform();
        
        for(long N = fromN; N < fromN + numruns; N++){
            double sum = 0.0;
            for(int n = 0; n < N; n++) sum += polyn(rv.getNoise(),Math.round(Math.pow(N,1.0/2 + 0.01)));
            System.out.println("N = " + N + ", mean = " + (sum/N));
        }
                
        
    }
    
}
