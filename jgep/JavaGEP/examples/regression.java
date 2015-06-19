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
import jGEP.*;
import java.util.Vector;
import java.util.Hashtable;

public class regression {
    public static void main(String args[]) {
        // single variable for the expressions
        char ts[] = {'a'};
        char fs[] = {'+','-','*'};

        // create the genome
        Genome g = new Genome(ts,fs,2,15);
        
        // make a PRNG
        java.util.Random rand = new java.util.Random();

        // sampler
        RouletteWheelSampler samp = 
            new RouletteWheelSampler(rand);

        // create a new population
        Population p = new Population(samp, 100, g);

        // populate with random individuals
        for (int i = 0; i < p.getSize(); i++) {
            int pos = 0;
            char c[] = new char[g.getGeneLength()];

            while (pos < g.getGeneLength()) {
                int r;
                if (pos < g.getHeadLength()) {
                    r = rand.nextInt(g.getSize());
                } else {
                    r = rand.nextInt(g.getNumTerminals());
                }

                if (r < g.getNumTerminals()) {
                    c[pos] = g.getTerminal(r);
                } else {
                    c[pos] = g.getFunction(r-g.getNumTerminals());
                }
                pos++;
            }

            String chromosome = new String(c);
            ArithmeticIndividual ai = new ArithmeticIndividual(chromosome,g);
            System.err.println("Created individual: "+chromosome);

            try {
                p.addIndividual(ai);
            } catch (Exception e) {
                System.err.println("Exception seeding initial population.");
                System.err.println(e.toString());
            }
        }

        // create a genetic operator object, with a bunch of probabilities
        // setup.
        GeneticOperators gops = new GeneticOperators(g,rand);
        gops.setP1Point(0.2);
        gops.setP2Point(0.2);
        gops.setPGRecomb(0.3);
        gops.setPGTrans(0.01);
        gops.setPISTrans(0.01);
        gops.setPMutate(0.1);
        gops.setPRISTrans(0.01);

        //
        // now enter the loop of:
        //   1. express
        //   2. test
        //   3. fitness (end if ideal)
        //   4. selection
        //   5. operators
        //      a. replication
        //      b. mutation
        //      c. is transpose
        //      d. ris transpose
        //      e. gene transpose
        //      f. 1pt recomb
        //      g. 2pt recomb
        //      h. gene recomb
        //   6. goto 1
        boolean keepGoing = true;
        Individual bestIndividual = null;

        //
        // fitness test points and the expected results.
        //
        Hashtable testValues[] = new Hashtable[10];
        double expected[] = new double[10];
        for (int i = 0; i < 10; i++) {
            testValues[i] = new Hashtable();
            testValues[i].put("a",new Double((double)(0.0-5.0)+(double)i));
            expected[i] = ((0.0-5.0)+(double)i)*((0.0-5.0)+(double)i); // x^2
        }

        double score[] = new double[p.getSize()];

        while (keepGoing) {
            Vector individuals = p.getIndividuals();
            for (int i = 0; i < p.getSize(); i++) {
                Individual ind = (Individual)individuals.elementAt(i);

                // step 1. express
                ExpressionNode roots[] = ind.express();

                // step 2. test
                score[i] = -100000.0;
                try {
                    for (int j = 0; j < 10; j++) {
                        // assume only roots[0] is meaningful
                        double val = roots[0].evaluate(testValues[j]);
                        if (val > -100000.0 && val < 100000.0) {
                            score[i] += Math.abs(val - expected[j]);
                        }
                    }
                } catch (Exception e) {
                    System.err.println("EXCEPTION: "+e);
                }
            }

            if (keepGoing) {
                // step 3. fitness weighting

                int bestIdx = 0;
                int worstIdx = 0;
                for (int idx = 1; idx < p.getSize(); idx++) {
                    if (score[idx] > score[bestIdx]) {
                        worstIdx = idx;
                    }
                    if (score[idx] < score[bestIdx]) {
                        bestIdx = idx;
                    }
                }

                if (score[bestIdx] == 0.0) {
                    keepGoing = false;
                    Individual theBest = 
                        (Individual)individuals.elementAt(bestIdx);
                    System.out.println(theBest.getChromosome());
                }
                
                System.out.println("WORST="+score[worstIdx]+
                                   "  BEST="+score[bestIdx]);

                // step 4. selection
                p.select(score, bestIdx);
                
                // step 5. genetic operators
                
            }
        }

    }
}
