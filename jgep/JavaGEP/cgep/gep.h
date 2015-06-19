/**
 * gep.h : c structs for gene expression programming experiments
 *
 * author: matthew sottile (matt@cs.uoregon.edu)
 *
 * Licensed under the terms of the GNU Public Licence.  See LICENCE.GPL 
 * for details.
 */
#ifndef __GEP_H__
#define __GEP_H__

/**
 * Structure containing information about the genome from which individuals
 * are coded.  The genome contains terminals and functions (non-terminals),
 * chromosome and gene lengths, and arity limits on functions.
 */
typedef struct gep_genome {
  /**
   * Array of terminal characters.
   */
  char *terminals;

  /**
   * Array of function (non-terminal) characters.
   */
  char *functions;

  /**
   * Number of terminal characters.
   */
  unsigned int  num_terminals;

  /**
   * Number of function (non-terminal) characters.
   */
  unsigned int  num_functions;

  /**
   * Length of the gene head.  This portion of a gene can contain both
   * terminal and non-terminal characters.
   */
  unsigned int  head_length;

  /**
   * Maximum arity for functions in this genome.  Used for calculation of
   * gene tail length.
   */
  unsigned int  max_arity;

  /** 
   * The length of an individual is determined as:
   *
   *   length = (head_length + (head_length * (max_arity - 1)) + 1) * 
   *              genes_per_chromosome
   *
   * The tail contains only terminal characters and is the minimum length 
   * necessary to contain the leaf nodes of the expression tree formed within
   * the head in the worst case.  In most cases, the tail will contain 
   * non-coding regions.
   */
  unsigned int  individual_length;

  /**
   * Genes per chromosome.
   */
  unsigned int  genes_per_chromosome;
} gep_genome_t;

/**
 * Structure containing information about a population.  This includes:
 * 
 *  - The genome from which individuals get their genetic material from.
 *  - The probability of various genetic operators (mutation, crossover, transposition)
 *  - The size of the population.
 *  - The individuals and their corresponding fitness metrics.
 *  - Seed and state information for random number generation used for sampling and
 *    genetic operators.  This must be managed carefully in the event that this
 *    code is used in a parallel or multi-threaded environment.
 */
typedef struct gep_population {
  /**
   * Genome from which individuals get their genetic material.
   */
  gep_genome_t *genome;

  /**
   * Probabilities of various genetic operators.
   *
   *   x1   = 1 point crossover.
   *   x2   = 2 point crossover.
   *   m    = mutation.
   *   rist = root insertion sequence transposition
   *   ist  = insertion sequence transposition
   *   gt   = gene transposition
   */
  float        p_x1, p_x2, p_m, p_rist, p_ist, p_gt; 

  /**
   *  seed for random number generation.  
   */
  long         rseed;

  /**
   * number of individuals
   */
  int           num_individuals;

  /**
   * individuals.  The individuals array is the current population.
   * the selectionbuffer is used for building the next generation
   * in the selection process.  The pointers are then switched - simple
   * double buffering to prevent copy cost and excess temporary allocation.
   */
  char         **individuals;
  char         **selectionbuffer;

  /**
   * double fitness values 
   */
  double       *fitnesses;

  /**
   * best fitness
   */
  double        max_fitness;

  /* note: additional population-wide parameters can go here to assist in
     making this library thread safe.  This might involve locks around the
     random number generator or similar things to prevent number from
     being drawn twice on different processors or threads. */
} gep_population_t;

#ifdef __cplusplus
extern "C" {
#endif
  
  gep_genome_t *gep_create_genome(int nt, char *tc, int nf, 
                                  char *fc, int hl, int ma);

  gep_population_t *gep_create_population(gep_genome_t *genome, int ni,
                                          float px1, float px2, float pm,
                                          double mf);

  void gep_next_generation(gep_population_t *p);

  void gep_destroy_genome(gep_genome_t *g);

  void gep_destroy_population(gep_population_t *p);
  
#ifdef __cplusplus
}
#endif

#endif /* __GEP_H__ */
