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

/**
 * Interface representing an expression node.
 *
 * @author   Matthew Sottile
 * @version  1.0
 */
public interface ExpressionNode {
    /**
     * The evaluate method takes a hashtable of values and returns an
     * object.  This treats the individual as a map from the tuple
     * space represented by the hashtable to a 1D space of objects.
     * For example, arithmetic expressions map from tuples of real numbers
     * to a single real value (R^n -> R).
     *
     * @param   values  The key/value pair of values being passed in.
     * @return          The object that this individual evaluates to
     *                  given the parameters in the hashtable.
     */
    public java.lang.Object evaluate(java.util.Hashtable values) 
        throws Exception;

    public String stringRepresentation();
}

