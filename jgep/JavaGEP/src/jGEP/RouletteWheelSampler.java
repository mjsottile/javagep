/** LANL:license
 * -------------------------------------------------------------------------
 * This SOFTWARE has been authored by an employee or employees of the
 * University of California, operator of the Los Alamos National Laboratory
 * under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.
 * The U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE.  The public may copy, distribute, prepare derivative works and
 * publicly display this SOFTWARE without charge, provided that this Notice
 * and any statement of authorship are reproduced on all copies.  Neither
 * the Government nor the University makes any warranty, express or implied,
 * or assumes any liability or responsibility for the use of this SOFTWARE.
 * If SOFTWARE is modified to produce derivative works, such modified
 * SOFTWARE should be clearly marked, so as not to confuse it with the
 * version available from LANL.
 * -------------------------------------------------------------------------
 * LANL:license
 * -------------------------------------------------------------------------
 */
package jGEP;

import java.util.*;

/**
 * Class that provides a Roulette Wheel Sampling method for the
 * selection step in GEP programs.
 *
 * @author   Matthew Sottile
 * @version  1.0
 */
public class RouletteWheelSampler implements Sampler {
    private Random r; // PRNG
    private double epsilon; // epsilon of tolerance for weight sum

    /**
     * Constructor.  Provide a pseudo-random number generator that is
     * or inherits from java.util.Random.
     *
     * @param  r    PRNG object.
     * @param  eps  Epsilon for tolerance for weight sum deviation from 1.0
     */
    public RouletteWheelSampler(Random r, double eps) {
        this.r = r;
        epsilon = eps;
    }

    /**
     * Sample from a population.  Given a set of weights (between 0 and 1),
     * return a list of equal length that contains indices of those
     * individuals selected.  For example, the weight array could be:
     *
     *      weights = 0.25 0.25 0.5
     *      return  = 1 3 3
     *      return  = 3 2 3
     *      etc...
     *
     * @param  weights  Array of individual weights
     * @return          Array of indices indicating sampled individuals.
     */
    public int[] sample(double weights[]) {
        int selected[] = new int[weights.length];
        
        //
        // first, verify that the weights make sense
        // (IE: sum up to 1.0)
        //
        double total = 0.0;
        for (int i = 0; i < weights.length; i++) {
            total += weights[i];
        }

        if (Math.abs(total-1.0) > epsilon) {
            System.err.println("Warning: weight sum outside tolerance ("+
                               total+")");
        }

        //
        // assuming weights are OK, now we sample.
        // O(n^2) algorithm, n=weights.length
        //
        int cur = 0;
        while (cur < weights.length) {
            double rnum = r.nextDouble();
            double accum = 0.0;

            for (int i = 0; i < weights.length; i++) {
                if (rnum == accum || rnum < accum+weights[i]) {
                    selected[cur] = i;
                    break;
                }
                accum += weights[i];
            }

            cur++;
        }
        
        // return array of selected indices.
        return selected;
    }
}
