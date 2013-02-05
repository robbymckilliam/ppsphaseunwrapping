/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pubsim.poly;

import pubsim.lattices.GeneralLattice;
import pubsim.lattices.decoder.BabaiNoLLL;

/**
 *
 * @author Robby McKilliam
 */
public class AmbiguityRemoverRectangular extends AmbiguityRemover {

    public AmbiguityRemoverRectangular(int m) {
        this.m = m;
        p = new double[m+1];
        M = constructBasisMatrix();
        GeneralLattice lattice = new GeneralLattice(M);
        np = new BabaiNoLLL(lattice);
    }
}
