/*
 */
package pubsim.poly;

import org.junit.*;
import static org.junit.Assert.*;
import pubsim.Util;
import pubsim.VectorFunctions;

/**
 *
 * @author Robby McKilliam
 */
public class AmbiguityRemoverRectangularTest {
    
    public AmbiguityRemoverRectangularTest() {
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
    
    @Test
    public void testInSideRectangle() {
        System.out.println("testInSideRectangle");
        int m = 5;
        AmbiguityRemoverRectangular instance = new AmbiguityRemoverRectangular(m);
        
        int iters = 20;
        for(int itr = 0; itr < iters; itr++){
            double[] p = VectorFunctions.randomGaussian(m+1, 0, 100);
            p = instance.disambiguate(p);
            //System.out.println(VectorFunctions.print(p));
            for(int i = 0; i < m+1; i++)
                assertTrue(Math.abs(p[i]) < 0.5/Util.factorial(i));
        }
    }
}
