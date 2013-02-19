/*
 */
package pubsim.poly;

import Jama.Matrix;
import bignums.BigRational;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author Robby McKilliam
 */
public class BigRationalAmbiguityRemoverTest {
    
    public BigRationalAmbiguityRemoverTest() {
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
    
    @Test
    public void testNextColumn() {
        System.out.println("testNextColumn");
        double tol = 1e-8;
        double[] testc = {1};
        BigRational[] c = {BigRational.ONE};
        for(int i = 0; i < c.length; i++) assertEquals(testc[i], c[i].doubleValue(), tol);
        testc = AmbiguityRemover.getNextColumn(testc);
        c = BigRationalAmbiguityRemover.getNextColumnBigRational(c);
        for(int i = 0; i < c.length; i++) assertEquals(testc[i], c[i].doubleValue(), tol);
        testc = AmbiguityRemover.getNextColumn(testc);
        c = BigRationalAmbiguityRemover.getNextColumnBigRational(c);
        for(int i = 0; i < c.length; i++) assertEquals(testc[i], c[i].doubleValue(), tol);
        testc = AmbiguityRemover.getNextColumn(testc);
        c = BigRationalAmbiguityRemover.getNextColumnBigRational(c);
        for(int i = 0; i < c.length; i++) assertEquals(testc[i], c[i].doubleValue(), tol);
        testc = AmbiguityRemover.getNextColumn(testc);
        c = BigRationalAmbiguityRemover.getNextColumnBigRational(c);
        for(int i = 0; i < c.length; i++) assertEquals(testc[i], c[i].doubleValue(), tol);
    }
    
    @Test
    public void testMatrixColumn() {
        System.out.println("test matrix");
        double tol = 1e-8;
        for(int m = 0; m <= 8; m++){
            BigRational[][] M = new BigRationalAmbiguityRemover(m).constructBasisMatrixBigRational();
            Matrix Mtest = new AmbiguityRemover(m).constructBasisMatrix();
            for(int i = 0; i < M.length; i++)
                for(int j = 0; j < M[i].length; j++)
                    assertEquals(M[i][j].doubleValue(), Mtest.get(i,j), tol);
        }
    }
}
