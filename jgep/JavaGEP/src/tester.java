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

import java.util.Hashtable;
import java.util.Vector;

public class tester {
    public static void reportError(String s) {
        System.err.println("ERROR: "+s);
    }
    
    //
    // test
    //
    public static void testloop() {
        /**
         * step 1 : create the genome with a set of terminals and
         *          functions.
         */
        char ts[] = {'a'};
        char fs[] = {'+','-','/','*'};
        Genome g = new Genome(ts,fs,2,50);

        /**
         * step 2 : create the population that progresses with a
         *          specific sampling technique
         */
        int popsize = 75;
        java.util.Random r = new java.util.Random();
        RouletteWheelSampler rws = new RouletteWheelSampler(r,0.00000001);
        Population p = new Population(rws, popsize, g);

        /**
         * step 3 : initialize the population with random individuals
         */
        for (int i = 0; i < popsize; i++) {
            ArithmeticIndividual ai = null;
            try {
                ai = new ArithmeticIndividual(g,1);
            } catch (Exception e) {
                reportError(e.toString());
            }
            if (ai != null) {
                ai.randomChromosome(r);
                try {
                    p.addIndividual(ai);
                } catch (Exception e) {
                    reportError(e.toString());
                }
            } else {
                reportError("null individual encountered.");
            }
        }

        //
        // set up the vector of test values
        //
        double ins[] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
        double outs[] = {0.841471,0.909297, 0.14112, -0.756802, 
                         -0.958924, -0.279415, 0.656987, 0.989358,
                         0.412118, -0.544021};
        Vector v = new Vector();
        for (int i = 0; i < ins.length; i++) {
            Hashtable ht = new Hashtable();
            ht.put("a",new Double(ins[i]));
            ht.put("Expected",new Double(outs[i]));
            v.addElement(ht);
        }

        //
        // fitness tester that we'll use later
        //
        Fitness fitness = new Fitness(v,100.0);

        //
        // set up the genetic operators - set the probabilities
        // for different operators.
        //
        GeneticOperators gops = new GeneticOperators(g,r);
        gops.setP1Point(0.3);
        gops.setP2Point(0.15);
        gops.setPGRecomb(0.0);
        gops.setPGTrans(0.0);
        gops.setPISTrans(0.05);
        gops.setPMutate(0.4);
        gops.setPRISTrans(0.03);        

        boolean keepGoing = true;
        
        while (keepGoing) {

            /**
             * step 4 : evaluate fitness
             */
            Vector individuals = p.getIndividuals();
            double fitnesses[] = new double[popsize];
            double bestFitness = -10000000.0; // really small
            int    bestIndex = 0;
            double weights[] = new double[popsize];
            for (int i = 0; i < popsize; i++) {
                weights[i] = (double)(1.0/(double)popsize);
            }
            
            for (int i = 0; i < popsize; i++) {
                Individual ind = (Individual)individuals.elementAt(i);
                fitnesses[i] = fitness.evaluate(ind);
                if (fitnesses[i] > bestFitness) {
                    bestFitness = fitnesses[i];
                    bestIndex = i;
                }
            }
            
            Individual ind = (Individual)individuals.elementAt(bestIndex);
            System.err.println("Best="+bestFitness);
            
            ExpressionNode roots[] = ind.express();            
            System.err.println("    ="+roots[0].stringRepresentation());
            
            /**
             * step 5 : selection and advancing the population
             */
            p.select(weights,bestIndex);
            individuals = p.getIndividuals();

            /**
             * step 6 : genetic operators
             */
            for (int i = 1; i < popsize; i++) {
                double draw;
                ind = (Individual)individuals.elementAt(i);

                //
                // mutation
                //
                draw = r.nextDouble();
                if (draw < gops.getPMutate()) {
                    String newChromosome = 
                        gops.mutate(ind.getChromosome(),1);
                    ind.setChromosome(newChromosome);
                }

                //
                // crossover
                //
                draw = r.nextDouble();
                if (draw < gops.getCrossoverRate()) {
                    String chromosomes[] = new String[2];
                    int which = r.nextInt(2);

                    //
                    // crossover requires another individual --
                    // pick another other than the one we're dealing with
                    // right now.  also leave index 0 alone - that's the
                    // best individual, and we want them to carry across
                    // generations in case none of the others are any better
                    //
                    int otherIdx = r.nextInt(popsize);
                    while (otherIdx == i || otherIdx == 0) {
                        otherIdx = r.nextInt(popsize);
                    }
                    
                    Individual otherind = 
                        (Individual)individuals.elementAt(otherIdx);

                    chromosomes[0] = ind.getChromosome();
                    chromosomes[1] = otherind.getChromosome();

                    switch (which) {
                    case 0: // 1pt
                        chromosomes = 
                            gops.OnePointRecombination(chromosomes);
                        break;
                    case 1: // 2pt
                        chromosomes = 
                            gops.TwoPointRecombination(chromosomes);
                        break;
                    case 2: // gene
                        System.err.println("Gene recombination not done.");
                        break;
                    default:
                        System.err.println("Not supposed to happen.");
                    }

                    ind.setChromosome(chromosomes[0]);
                    otherind.setChromosome(chromosomes[1]);
                }
            }
            
            /**
             * Goto step 4 if fitness not converged to acceptable region
             */
            if (bestFitness > 999.0) {
                keepGoing = false;
            }
        }
    }

    //
    // main
    //
    public static void main(String args[]) {
	testloop();
    }
}
