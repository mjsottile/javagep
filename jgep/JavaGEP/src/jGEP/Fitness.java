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

import java.util.Vector;
import java.util.Hashtable;

/**
 * This class is intended to serve as a fitness testing harness.  Given
 * a set of test cases and individuals, it will manage testing and result
 * tracking.
 *
 * @author    Matthew Sottile
 * @version   1.0
 */
public class Fitness {
    protected Vector testValues;
    private   double maxFitness;

    /**
     * Constructor
     *
     * @param tests Vector of test value sets.  Each element of the vector
     *              should be a hash table mapping a name to a value of the
     *              type expected by an expression node of the corresponding
     *              individuals.
     * @param max   Maximum possible fitness value.
     */
    public Fitness(Vector tests, double max) {
        testValues = tests;
        maxFitness = max;
    }

    /**
     * Maximum fitness.
     *
     * @return      Value of maximum fitness.
     */
    public double getMaxFitness() {
        return maxFitness;
    }

    /**
     * Provide the vector of test values that fitness will be
     * evaluated against.  For a given entry (a hashtable of terminal
     * symbols and their values), one entry of the "expected" output
     * is provided, under the key "Expected".  This is the only
     * reserved key in the hashtables, and MUST be present.
     * 
     * @param  tvs   Vector of test values.
     */
    public void setTestValues(Vector tvs) {
        testValues = tvs;
        for (int i = 0; i < tvs.size(); i++) {
            Hashtable vs = (Hashtable)tvs.elementAt(i);
            Double e = (Double)vs.get("Expected");
            if (e == null) {
                System.err.println("ERROR: hashtable missing `Expected' key.");
            }
        }
    }

    public double evaluate(Individual ind) {
        ExpressionNode roots[] = ind.express();
        double fVals[] = new double[testValues.size()];
        double fval = 0.0;

        // for now we test the first node since we haven't dealt with
        // connectives for the forest of trees.
        for (int i = 0; i < testValues.size(); i++) {
            Hashtable vals = (Hashtable)testValues.elementAt(i);
            
            try {
                fVals[i] = ((Double)roots[0].evaluate(vals)).doubleValue();
            } catch (Exception e) {
                fVals[i] = -1000000000.0;
            }
            
            double expected = ((Double)(vals.get("Expected"))).doubleValue();
            
            fval += (maxFitness - Math.abs(fVals[i] - expected));
        }
        
        return fval;
    }
}
