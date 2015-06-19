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
 * This class represents the characters that make up the genes
 * as functions and terminals.  Very basic, but allows us to
 * hide the terminals/functions arrays and do anything fancy
 * behind the scenes later that might be of interest.
 *
 * @author  Matthew Sottile
 * @version 1.0
 */
public class Genome {
    private char terminals[]; // set of terminal characters
    private char functions[]; // set of function characters
    private int  nt, nf;      // terminal and function counts
    private int  maxArity;    // maximum arity of any function
    private int  h, t;        // head and tail lengths
    
    /**
     * Constructor.
     *
     * @param  ts   Array of characters representing terminal symbols.
     * @param  fs   Array of characters representing functions.
     * @param  ma   Maximum arity of the functions.  For example, if '+'
     *              has the maximum arity, this is 2 since + is binary.
     * @param  hl   Head length in a gene. 
     */
    public Genome(char ts[], char fs[], int ma, int hl) {
	terminals = ts;
	functions = fs;
	nt = ts.length;
	nf = fs.length;
        maxArity = ma;
        this.h = hl;
        this.t = h*(maxArity-1) + 1;
    }

    /**
     * Length of a head in a gene.  The head can contain both terminal and
     * function characters.
     *
     * @return   The head length.
     */
    public int getHeadLength() {
        return h;
    }

    /**
     * Length of the tail in a gene.  The tail can only contain terminal
     * characters.  This is a requirement to ensure that any resulting
     * chromosome after a genetic operator is applied is a valid
     * chromosome.
     *
     * @return   The tail length.
     */
    public int getTailLength() {
        return t;
    }

    /**
     * Return the length of a gene (tail + head).
     *
     * @return   The length of a gene.
     */
    public int getGeneLength() {
        return t+h;
    }

    /**
     * Return the maximum arity of the function set in this genome.
     * For the basic case of arithmetic (+,-,*,/), we have a maximum
     * arity of 2.  
     *
     * @return   Maximum arity of the function set.
     */
    public int getMaxArity() {
        return maxArity;
    }
    
    /**
     * Return the size of the genome (the number of possible characters).
     * This is the sum of the number of terminals and the number of functions.
     *
     * @return  Genome size.
     */
    public int getSize() {
        return nt+nf;
    }

    /**
     * Return the number of terminal symbols available in the genome.
     *
     * @return   Number of terminals.
     */
    public int getNumTerminals() {
	return nt;
    }

    /**
     * Return the number of function symbols available in the genome.
     * 
     * @return   Number of functions.
     */
    public int getNumFunctions() {
	return nf;
    }

    /**
     * Return the terminal symbol at position n.  This simply a way of
     * uniquely identifying symbols consistently across all objects.
     * (Similar to saying 'b' is the 2nd (index 1) element of the english
     * alphabet).
     *
     * @param  n  The index into the terminal array.
     * @return    The character at index n.
     */
    public char getTerminal(int n) {
	return terminals[n];
    }

    /**
     * Return the function symbol at index n.  For details on what this 
     * means, look at getTerminal().
     * @see jGEP.Genome#getTerminal(int n) getTerminal(n)
     *     
     * @param  n  The index into the function array.
     * @return    The character at index n.
     */
    public char getFunction(int n) {
	return functions[n];
    }

    /**
     * Return a boolean true or false whether the character is a function
     * or not.
     *
     * @param   c  The character.
     * @return     Boolean if it is a function.
     */
    public boolean isFunction(char c) {
        for (int i = 0; i < nf; i++) {
            if (functions[i] == c) {
                return true;
            }
        }

        return false;
    }
}
