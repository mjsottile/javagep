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
 * Class that provides a Stocastic Universal Sampling method for the
 * selection step in GEP programs.
 *
 * @author   Matthew Sottile
 * @version  1.0
 */
public class StochasticUniversalSampler implements Sampler {
    private Random r; // PRNG

    /**
     * Constructor.  Provide a pseudo-random number generator that is
     * or inherits from java.util.Random.
     *
     * @param  r  PRNG object.
     */
    public StochasticUniversalSampler(Random r) {
        this.r = r;
    }

    /**
     * Sample from a population.  Weights are IGNORED, but the array is used
     * to determine the size of the population.
     *
     * @param  weights  Array of individual weights
     * @return          Array of indices indicating sampled individuals.
     */
    public int[] sample(double weights[]) {
        int selected[] = new int[weights.length];
        
        //
        // sampling step
        //
        for (int i = 0; i < selected.length; i++) {
            double samp = r.nextDouble();
            
            selected[i] = 
                (int)java.lang.Math.floor(samp*(double)selected.length);
        }
        
        // return array of selected indices.
        return selected;
    }
}
