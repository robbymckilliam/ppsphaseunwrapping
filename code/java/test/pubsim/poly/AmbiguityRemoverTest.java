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
import pubsim.VectorFunctions;

/**
 *
 * @author Robby McKilliam
 */
public class AmbiguityRemoverTest {

    public AmbiguityRemoverTest() {
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
    public void throwsExceptionWhenIncorrectSize() {
        System.out.println("throwsExceptionWhenIncorrectSize");
        int m = 4;
        double[] p = {3, 2, 4, 5};
        AmbiguityRemover instance = new AmbiguityRemover(m);

        boolean caught = false;
        try{
            instance.disambiguate(p);
        } catch(RuntimeException e){
            caught = true;
        }
        assertTrue(caught);
    }

    @Test
    public void testDisambiguate() {
        System.out.println("testDisambiguate");
        double[] p1 = {1.3, -2.1};
        int m = p1.length - 1;
        AmbiguityRemover instance = new AmbiguityRemover(m);
        p1 = instance.disambiguate(p1);
        double[] y1 = {0.3, -0.1};
        assertTrue(VectorFunctions.distance_between2(p1, y1) < 0.000001);

        double[] p2 = {1.3, -2.1, 2.6};
        m = p2.length - 1;
        instance = new AmbiguityRemover(m);
        p2 = instance.disambiguate(p2);
        double[] y2 = {0.3, 0.4, 0.1};
        assertTrue(VectorFunctions.distance_between2(p2, y2) < 0.000001);
        
    }


        /**
     * Test of setSize method, of class AmbiguityRemover.
     */
    @Test
    public void testNextColumn() {
        System.out.println("testNextColumn");
        double[] c = VectorFunctions.eVector(0, 1);
        //System.out.println("c = " + VectorFunctions.print(c));
        c = AmbiguityRemover.getNextColumn(c);
        //System.out.println("c = " + VectorFunctions.print(c));
        double[] y1 = {0, 1};
        assertTrue(VectorFunctions.distance_between2(c, y1) < 0.000001);

        c = AmbiguityRemover.getNextColumn(c);
        //System.out.println("c = " + VectorFunctions.print(c));
        double[] y2 = {0, 1.0/2.0, 1.0/2.0};
        assertTrue(VectorFunctions.distance_between2(c, y2) < 0.000001);

        c = AmbiguityRemover.getNextColumn(c);
        //System.out.println("c = " + VectorFunctions.print(c));
        double[] y3 = {0, 1.0/3.0, 1.0/2.0, 1.0/6};
        assertTrue(VectorFunctions.distance_between2(c, y3) < 0.000001);

    }

}