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

import java.util.Random;

/**
 * Class representing a single individual.
 *
 * @author    Matthew Sottile
 * @version   1.0
 */
abstract public class Individual {
    protected String chromosome; /* the chromosome of this individual.
                                    this is one valid string in the
                                    space of valid strings (the genome). */
    protected Genome genome;     /* the genome is the space from which
                                    chromosomes are derived.  */
    protected int    genes;      /* number of genes in the chromosome */
    
    /**
     * Constructor.  An individual is created with a chromosome containing
     * their genes.
     *
     * @param  c   The chromosome.
     * @param  g   The genome
     * @param  gc  Number of genes in the chromosome
     */
    public Individual(String c, Genome g, int gc) throws Exception {
        //
        // check for errors on call
        //
        if (g == null) {
            throw new Exception("Cannot create individual with null Genome.");
        }
        if (c == null) {
            throw new Exception("Call other constructor.");
        }
        if (gc <= 0) {
            throw new Exception("Bogus gene count.");
        }
        
        //
        // sanity check
        //
        if (c.length() % g.getGeneLength() != 0) {
            throw new Exception("Chromosome wrong length.");
        }
        
        this.chromosome = new String(c);
        genome = g;
        genes = gc;
    }
    
    /**
     * Constructor.  This constructor creates an individual with no
     * genetic code yet.
     * 
     * @param  g   The genome
     * @param  gc  Number of genes in the chromosome
     */
    public Individual(Genome g, int gc) throws Exception {
        if (g == null) {
            throw new Exception("Cannot create individual with null genome.");
        }
        if (gc <= 0) {
            throw new Exception("Bogus gene count.");
        }
        
        chromosome = null;
        genome = g;
        genes = gc;
    }
    
    /**
     * Repicates this individual.  Essentially a clone.
     *
     * @return   A new individual cloning this one.
     */
    abstract public Individual replicate();

    /**
     * Set the chromosome for this individual to a random setting.
     */
    public void randomChromosome(Random r) {
        if (genome == null) {
            System.err.println("Cannot generate chromosome from null genome.");
            return;
        }
        
        if (r == null) {
            System.err.println("Pass in non-null PRNG.");
            return;
        }
        
        chromosome = new String("");
        int nf = genome.getNumFunctions();
        int nt = genome.getNumTerminals();
        
        for (int i = 0; i < genes; i++) {
            for (int j = 0; j < genome.getHeadLength(); j++) {
                // flip for function or terminal
                int f = r.nextInt(2);
                if (f == 0) {
                    chromosome += genome.getFunction(r.nextInt(nf));
                } else {
                    chromosome += genome.getTerminal(r.nextInt(nt));
                }
            } 
            for (int j = 0; j < genome.getTailLength(); j++) {
                chromosome += genome.getTerminal(r.nextInt(nt));
            }
        }
    }
    
    /**
     * Set the genome.
     *
     * @param  g   The genome object.
     */
    public void setGenome(Genome g) {
        genome = g;
    }

    /**
     * Set the chromosome.
     *
     * @param  c   The string containing the chromosome.
     */
    public void setChromosome(String c) {
        chromosome = new String(c);
    }

    /**
     * Return the chromosome.
     *
     * @return  The chromosome.
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * Expression is the process of translating from a flat chromosome
     * to the functional, structurally meaningful construct that it
     * represents.  A single chromosome can contain an arbitrary number
     * of genes, each expressed into a distinct structure.  No
     * explicit relation between genes and their expressed structures is
     * encoded in the chromosomes.
     *
     * @return   Array of ExpressionNode objects.  The individual is
     *           expected to provide a derived class from ExpressionNode
     *           that represents the structures specific to that class
     *           of individual.  This is easily achieved through the
     *           use of inner classes (which also provide a means to logically
     *           bind the expression node subclass to the individual 
     *           subclass.)
     */
    abstract public ExpressionNode[] express();
}
