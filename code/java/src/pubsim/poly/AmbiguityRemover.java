/*
 */
package pubsim.poly;

import Jama.Matrix;
import java.io.Serializable;
import pubsim.VectorFunctions;
import pubsim.lattices.Lattice;
import pubsim.lattices.NearestPointAlgorithmInterface;
import pubsim.lattices.decoder.SphereDecoder;

/**
 * This uses m nearest lattice point approach to remove the
 * ambiguities inherent in polynomial phase estimation.
 * @author Robby McKilliam
 */
public class AmbiguityRemover implements Serializable{

    protected int m;
    protected Matrix M;
    double[] p;
    NearestPointAlgorithmInterface np;

    protected AmbiguityRemover() {
    }

    /**
     * Public constructor.  Set the order of the polynomial phase signal.
     * @param m = the number of parameters ie. order of polynomial + 1
     */
    public AmbiguityRemover(int m) {
        this.m = m;
        p = new double[m+1];
        M = constructBasisMatrix();
        Lattice lattice = new Lattice(M);
        np = new SphereDecoder(lattice);
    }
    
    protected final Matrix constructBasisMatrix() {
        Matrix B = new Matrix(m+1, m+1);
        double[] c = VectorFunctions.eVector(0, 1);
        for (int j = 0; j < m+1; j++) {
            for (int i = 0; i <= j; i++) {
                B.set(i, j, c[i]);
            }
            c = getNextColumn(c);
        }
        return B;
    }

    /**
     * Recursively generates the columns of the ambiguity lattice
     * generator matrix.  Essentially generates coefficients of the
     * integer valued polynomials.
     * @param c previous column.
     * @return next column
     */
    protected static double[] getNextColumn(double[] c) {
        double[] r = new double[c.length + 1];

        //shift c
        for (int i = 1; i < r.length; i++) {
            r[i] = c[i - 1];
        }

        for (int i = 0; i < c.length; i++) {
            r[i] += c[i] * (c.length - 1);
        }

        for (int i = 0; i < r.length; i++) {
            r[i] /= c.length;
        }
        return r;
    }

    /**
     * This function removes the ambiguities from polynomial
     * phase signals.
     * @param p the parameter to remove ambiguity from
     * @return parameter in identifiable range.
     */
    public double[] disambiguate(double[] p) {
        if (m+1 != p.length) {
            throw new RuntimeException("Parameter vector p is not the correct size.");
        }
        np.nearestPoint(p);
        double[] np = this.np.getLatticePoint();

        for (int i = 0; i < p.length; i++) {
            this.p[i] = p[i] - np[i];
        }
        return this.p;
    }

    /** Return the basis matrix for the ambiguity lattice */
    public Matrix getBasisMatrix(){
        return M;
    }
    
    /** Return the volume of the identifiable region */
    public double getRegionVolume(){
        Matrix Mt = M.transpose();
        return Math.sqrt(Mt.times(M).det());
    }
}
