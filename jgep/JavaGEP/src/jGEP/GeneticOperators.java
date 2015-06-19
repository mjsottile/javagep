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
 * This class contains all operators on genes that can occur either
 * individually (such as mutation) or during reproduction (such as
 * recombination).  Users who wish to customize any of these can
 * extend this class and re-implement whichever methods they like.
 * A population will have an instance of this object, since it operates
 * on that specific population.  Probabilities of the operators occuring
 * also exist in this object instead of elsewhere.  This allows other 
 * operators to be created and their probabilities will accompany them.
 *
 * @author   Matthew Sottile
 * @version  1.0
 */
public class GeneticOperators {
    private Genome g;                         // Genome that this object
                                              // recognizes.
    private java.util.Random r;               // Pseudo-random number generator

    //
    // operator probabilities
    //
    private double pMutate;                      // mutation
    private double p1Point, p2Point, pGRecomb;   // recombination
    private double pGTrans, pISTrans, pRISTrans; // transposition

    /**
     * Constructor.
     *
     * @param   g   The genome that these operators act upon.
     * @param   r   A PRNG object representing a random variable that
     *              the operators can sample from.
     */
    public GeneticOperators(Genome g, java.util.Random r) {
        this.r = r;
        this.g = g;
    }

    /**
     * Return the crossover rate - this is the sum of the rates for all
     * types of recombination (gene, 1pt, 2pt).
     *
     * @return   The crossover rate for the genome.
     */
    public double getCrossoverRate() {
        return pGRecomb+p1Point+p2Point;
    }

    /**
     * Set the probability of gene recombination.
     * 
     * @param  p  The probability.
     */
    public void setPGRecomb(double p) {
        pGRecomb = p;
    }

    /**
     * Get the probability of gene recombination.
     *
     * @return   The probability.
     */
    public double getPGRecomb() {
        return pGRecomb;
    }

    /**
     * Set the probability of gene transposition.
     * 
     * @param  p  The probability.
     */
    public void setPGTrans(double p) {
        pGTrans = p;
    }

    /**
     * Set the probability of insertion sequence transposition.
     * 
     * @param  p  The probability.
     */
    public void setPISTrans(double p) {
        pISTrans = p;
    }

    /**
     * Set the probability of root insertion sequence transposition.
     * 
     * @param  p  The probability.
     */
    public void setPRISTrans(double p) {
        pRISTrans = p;
    }

    /**
     * Get the probability of gene transposition.
     *
     * @return   The probability.
     */
    public double getPGTrans() {
        return pGTrans;
    }

    /**
     * Get the probability of insertion sequence transposition.
     *
     * @return   The probability.
     */
    public double getPISTrans() {
        return pISTrans;
    }

    /**
     * Get the probability of root insertion sequence transposition.
     *
     * @return   The probability.
     */
    public double getPRISTrans() {
        return pRISTrans;
    }

    /**
     * Setter for the mutation probability.
     *
     * @param   p   The probability value.
     */
    public void setPMutate(double p) {
        //
        // if the probability is bogus, set it to 1.0.
        // later this might be an exception.
        //
        if (!(p >= 0.0 && p <= 1.0)) {
            pMutate = 1.0;
        } else {
            pMutate = p;
        }
    }

    /**
     * Getter for the mutation probability.
     *
     * @return   The mutation probability.
     */
    public double getPMutate() {
        return pMutate;
    }

    /**
     * Setter for the one-point recombination probability.
     *
     * @param   p   The one-point recombination probability.
     */
    public void setP1Point(double p) {
        //
        // if the probability is bogus, set it to 1.0.
        // later this might be an exception.
        //
        if (!(p >= 0.0 && p <= 1.0)) {
            p1Point = 1.0;
        } else {
            p1Point = p;
        }
    }

    /**
     * Getter for the two-point recombination probability.
     *
     * @return   The two-point recombination probability.
     */
    public double getP1Point() {
        return p1Point;
    }

    /**
     * Setter for the two-point recombination probability.
     *
     * @param   p   The two-point recombination probability.
     */
    public void setP2Point(double p) {
        //
        // if the probability is bogus, set it to 1.0.
        // later this might be an exception.
        //
        if (!(p >= 0.0 && p <= 1.0)) {
            p2Point = 1.0;
        } else {
            p2Point = p;
        }
    }

    /**
     * Getter for the two-point recombination probability.
     *
     * @return   The two-point recombination probability.
     */
    public double getP2Point() {
        return p2Point;
    }

    /**
     * Insertion-sequence transposition.
     *
     * @param  s      The string to perform IS transposition on.
     * @return        The string after transposition.
     */
    public String IStranspose(String s) {
        // step 1, pick the insertion sequence
        int isLength = r.nextInt(g.getHeadLength()-1) + 1;
        int isStart = r.nextInt(s.length()-isLength);
        String is_string = s.substring(isStart,isStart+isLength);

        // pick where the insertion sequence is going to go.  It MUST
        // go in the head of a gene, and must not start at the root element
        // of a gene.
        int targetGene = r.nextInt(s.length() / g.getGeneLength());
        int genePosition = r.nextInt(g.getHeadLength() - 1) + 1;

        // now pull the chromosome apart and paste the transposition element
        // in.
        String cHead = 
            s.substring(0,targetGene*g.getGeneLength()+genePosition);
        String cTail =
            s.substring(targetGene*g.getGeneLength()+genePosition+isLength);

        // return the chromosome reassembled with transposition complete
        //
        return cHead+is_string+cTail;
    }

    /** 
     * Root-insertion-sequence transposition.
     *
     * @param s     The string to perform RIS transposition on.
     * @return      The string after transposition.
     */
    public String RIStranspose(String s) {
        // pick a gene
        int gene = r.nextInt(s.length() / g.getGeneLength());
        String sgene = s.substring(gene*g.getGeneLength(),g.getGeneLength());

        // pick the RIS - find a function and pick sequence downstream from 
        // there.
        int pos = r.nextInt(g.getHeadLength());
        while ((pos < g.getHeadLength()) && !(g.isFunction(s.charAt(pos)))) {
            pos++;
        }

        // no sequence found - do nothing.
        if (pos == g.getHeadLength()) {
            return new String(s);
        }
        
        // sequence found - get it.
        int ris_length = r.nextInt(g.getGeneLength()-(pos-1))+1;
        String ris = sgene.substring(pos,pos+ris_length);
        String leftover = sgene.substring(ris.length());
        sgene = ris+leftover;

        String cHead, cTail;

        // chop the gene out of the chromosome.
        if (gene != 0) {
            cHead = s.substring(0,gene*g.getGeneLength());
            if ((gene+1)*g.getGeneLength() == s.length()) {
                cTail = "";
            } else {
                cTail = s.substring((gene+1)*g.getGeneLength());
            }
        } else {
            cHead = "";
            cTail = s.substring(g.getGeneLength());
        }

        // reassemble the chromosome

        return cHead+sgene+cTail;
    }

    /**
     * Gene transposition.
     *
     * @param  s   The string to perform gene transposition on.
     * @return     The string after transposition.
     */
    public String GeneTranspose(String s) {
        // find number of genes in chromosome
        int numGenes = s.length() / g.getGeneLength();

        // mono-genic is a NOP
        if (numGenes == 1) {
            return new String(s);
        }

        // since we have more than one gene, pick two to transpose
        int g1 = r.nextInt(numGenes);
        int g2 = r.nextInt(numGenes);

        // transposing same gene is a NOP
        if (g1 == g2) {
            return new String(s);
        }
        
        // at this point, transposing two different genes.
        String gs1, gs2;

        // extract the gene strings for each
        //
        if (g1 == numGenes-1) {
            gs1 = s.substring(g1 * g.getGeneLength());
        } else {
            gs1 = s.substring(g1 * g.getGeneLength(), g.getGeneLength());
        }

        if (g2 == numGenes-1) {
            gs2 = s.substring(g2 * g.getGeneLength());
        } else {
            gs2 = s.substring(g2 * g.getGeneLength(), g.getGeneLength());
        }

        // this part below is to extract the material around the genes from
        // the chromosome so we can put it back together with the genes
        // in their new locations.
        //
        String before, middle, after;
        int lo, hi;
        String slo, shi;

        if (g1 < g2) { lo = g1; hi = g2; } else { lo = g2; hi = g1; }

        if (lo == 0) {
            before = "";
        } else {
            before = s.substring(0,lo*g.getGeneLength());
        }

        if (hi == numGenes-1) {
            after = "";
        } else {
            after = s.substring((hi+1)*g.getGeneLength());
        }

        if (hi-lo == 1) {
            middle = "";
        } else {
            middle = s.substring((lo+1)*g.getGeneLength(),
                                 hi*g.getGeneLength());
        }

        if (g1 == lo) {
            slo = gs1; shi = gs2;
        } else {
            slo = gs2; shi = gs1;
        }

        // low gene moves to high position, high to low.  very basic.
        //
        return before+shi+middle+slo+after;
    }

    /**
     * Gene recombination.  Involves picking a gene in each chromosome
     * and swapping them.
     *
     * @param s    Array of two chromosome strings to perform gene
     *             recombination on.
     * @return     Array of two chromosome strings after recombination.
     */
    public String[] GeneRecombination(String s[]) {
        String rets[] = new String[2];

        // figure out how many genes are in a chromosome
        int numGenes = (int)s[0].length()/g.getGeneLength();
        
        if (!(numGenes > 1)) {
            // chromosomes are mono-genic.  nothing to do, so we just swap
            // the chromosomes - no real effect.
            rets[0] = s[1];
            rets[1] = s[0];
        } else {
            int g1 = r.nextInt(numGenes); // gene from chromosome 1
            int g2 = r.nextInt(numGenes); // gene from chromosome 2

            String gs1, gs2;

            // extract the genes to swap
            if (g1 != numGenes-1) { 
                gs1 = s[0].substring((g1*g.getGeneLength()),g.getGeneLength());
            } else {
                gs1 = s[0].substring((g1*g.getGeneLength()));
            }

            if (g2 != numGenes-1) {
                gs2 = s[1].substring((g2*g.getGeneLength()),g.getGeneLength());
            } else {
                gs2 = s[1].substring((g2*g.getGeneLength()));
            }

            // extract the genes before and after those moving so we can
            // paste the new chromosomes together.
            String h1 = "";
            if (g1 != 0) {
                h1 = s[0].substring(0,(g1*g.getGeneLength()));
            }
            String t1 = s[0].substring((g1+1)*g.getGeneLength());

            String h2 = "";
            if (g2 != 0) {
                h2 = s[1].substring(0,(g2*g.getGeneLength()));
            }
            String t2 = s[1].substring((g2+1)*g.getGeneLength());

            // paste the new chromosomes together
            rets[0] = h1+gs2+t1;
            rets[1] = h2+gs1+t2;
        }
            
        return rets;
    }

    /**
     * Given a string and the number of mutations to cause, return
     * a string with that number of mutations.  Note that an external
     * object will get the probability of mutation FROM this object, and then
     * use it to determine the numMutations parameter.
     * 
     * @param  s              The chromosome to mutate.
     * @param  numMutations   The number of mutations.
     * @return                The chromosome with mutations.
     */
    public String mutate(String s, int numMutations) {
        String src = new String(s);

        for (int i = 0; i < numMutations; i++) {
            int gene = r.nextInt(s.length() / g.getGeneLength());
            int pos = r.nextInt(g.getHeadLength());
            int v = r.nextInt(g.getSize());
            char schars[] = src.toCharArray();
            int offset = g.getGeneLength() * gene;

            if (v >= g.getNumFunctions()) {
                schars[pos+offset] = g.getTerminal(v-g.getNumFunctions());
            } else {
                schars[pos+offset] = g.getFunction(v);
            }

            src = new String(schars);
        }

        return src;
    }

    /**
     * Given a two-element array containing two string chromosomes, 
     * perform a one point recombination between them and return a new
     * two-element array with the new chromosomes.  This is called by
     * an external object that has used the probability stored in this
     * object to determine whether or not a one-point recombination
     * should occur.
     *
     * @param    s    The array of chromosomes.  Assume length is two.
     *                eventually, if this is not true we will throw an
     *                exception.
     * @return        The two new chromosomes after recombination.
     */
    public String[] OnePointRecombination(String s[]) {
        String front1, back1;
        String front2, back2;
        String outStrings[] = new String[2];

        outStrings[0] = new String();
        outStrings[1] = new String();

        // pick point to do recombination
        int pos = r.nextInt(s[0].length());

        front1 = s[0].substring(0,pos);
        back1  = s[0].substring(pos);

        front2 = s[1].substring(0,pos);
        back2  = s[1].substring(pos);

        outStrings[0] = front1+back2;
        outStrings[1] = front2+back1;

        return outStrings;
    }

    /**
     * Given a two-element array containing two string chromosomes, 
     * perform a two point recombination between them and return a new
     * two-element array with the new chromosomes.  This is called by
     * an external object that has used the probability stored in this
     * object to determine whether or not a two-point recombination
     * should occur.
     *
     * @param    s    The array of chromosomes.  Assume length is two.
     *                eventually, if this is not true we will throw an
     *                exception.
     * @return        The two new chromosomes after recombination.
     */
    public String[] TwoPointRecombination(String s[]) {
        int pos1 = r.nextInt(s[0].length());
        int pos2 = r.nextInt(s[0].length());
        int hi, lo;

        // crossover points
        if (pos1 > pos2) {
            hi = pos1; lo = pos2;
        } else {
            hi = pos2; lo = pos1;
        }

        String front1, mid1, back1;
        String front2, mid2, back2;
        String outGenes[] = new String[2];
        outGenes[0] = new String();
        outGenes[1] = new String();

        front1 = s[0].substring(0,lo);
        mid1   = s[0].substring(lo,hi);
        back1  = s[0].substring(hi);

        front2 = s[1].substring(0,lo);
        mid2   = s[1].substring(lo,hi);
        back2  = s[1].substring(hi);

        outGenes[0] = front1+mid2+back1;
        outGenes[1] = front2+mid1+back2;

        return outGenes;
    }
}
