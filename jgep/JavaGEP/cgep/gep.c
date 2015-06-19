/**
 * FILE: gep.c
 *
 * author: matthew sottile (matt@cs.uoregon.edu)
 *
 * Licensed under the terms of the GNU Public Licence.  See LICENCE.GPL 
 * for details.
 */
#include "gep.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* allow assertions to be disabled */
#ifdef NOASSERTS
#undef assert
#define assert(a) { } 
#endif /* NOASSERTS */

gep_genome_t *gep_create_genome(int nt, char *tc, int nf, char *fc, 
                                int hl, int ma) {
  gep_genome_t *g;
  
  /* sanity checks */
  assert(nt > 0);
  assert(nf > 0);
  assert(hl > 0);
  assert(ma > 0);

  /* make sure arrays are non-null */
  assert(tc != NULL);
  assert(fc != NULL);

  g = (gep_genome_t *)calloc(1, sizeof(gep_genome_t));
  assert(g != NULL);

  g->num_terminals = nt;
  g->num_functions = nf;
  g->head_length = hl;
  g->max_arity = ma;
  g->individual_length = (hl * (ma - 1)) + 1 + hl;

  g->terminals = (char *)malloc(sizeof(char) * nt);
  assert(g->terminals != NULL);
  memcpy(g->terminals, tc, nt*sizeof(char));

  g->functions = (char *)malloc(sizeof(char) * nf);
  assert(g->functions != NULL);
  memcpy(g->functions, fc, nf*sizeof(char));

  return g;
}

gep_population_t *gep_create_population(gep_genome_t *genome, int ni,
                                        float px1, float px2, float pm,
                                        double mf) {
  gep_population_t *p;
  int i, j, n;
  double rnum_f;

  /* need a genome! */
  assert(genome != NULL);

  /* more than 0 individuals */
  assert(ni > 0);

  /* make sure probabilities make sense */
  assert(px1 >= 0.0);
  assert(px1 <= 1.0);
  assert(px2 >= 0.0);
  assert(px2 <= 1.0);
  assert(pm >= 0.0);
  assert(pm <= 1.0);

  p = (gep_population_t *)malloc(sizeof(gep_population_t));
  assert(p != NULL);

  p->num_individuals = ni;
  p->p_x1 = px1;
  p->p_x2 = px2;
  p->p_m = pm;
  p->genome = genome;

  p->fitnesses = (double *)malloc(sizeof(double)*ni);
  assert(p->fitnesses != NULL);

  p->individuals = (char **)malloc(sizeof(char *)*ni);
  assert(p->individuals != NULL);

  p->selectionbuffer = (char **)malloc(sizeof(char *)*ni);
  assert(p->selectionbuffer != NULL);

  for (i = 0; i < ni; i++) {
    p->individuals[i] = (char *)malloc(sizeof(char)*genome->individual_length);
    assert(p->individuals[i] != NULL);
    p->selectionbuffer[i] = 
      (char *)malloc(sizeof(char)*genome->individual_length);
    assert(p->selectionbuffer[i] != NULL);
  }

  /* create individuals */
  for (i = 0; i < ni; i++) {
    for (j = 0; j < genome->individual_length; j++) {
      rnum_f = (double)random() / (double)RAND_MAX;

      /* treat head different from tail */
      if (j < genome->head_length) {
        n = (int)(floor((genome->num_terminals + genome->num_functions)*rnum_f));
        if (n >= genome->num_terminals) {
          (p->individuals[i])[j] = 
            genome->functions[n - genome->num_terminals];
        } else {
          (p->individuals[i])[j] = genome->terminals[n];
        }
      } else {
        (p->individuals[i])[j] = 
          genome->terminals[(int)(floor(genome->num_terminals * rnum_f))];
      }
    }
  }

  return p;
}

void gep_destroy_genome(gep_genome_t *g) {
  if (g == NULL) return;

  if (g->terminals != NULL) free(g->terminals);
  if (g->functions != NULL) free(g->functions);
  free(g);
}

void gep_destroy_population(gep_population_t *p) {
  int i;

  if (p == NULL) return;

  if (p->fitnesses != NULL) free(p->fitnesses);
  if (p->individuals != NULL) {
    for (i=0; i<p->num_individuals; i++)
      if (p->individuals[i] != NULL) free(p->individuals[i]);
    free(p->individuals);
  }

  free(p);
}

void gep_next_generation(gep_population_t *p) {
  double rv;
  int i,a,b,j,aa,bb,rpos;
  int *idx;
  char *tmp;

  if (p == NULL) return;

  idx = (int *)malloc(sizeof(int)*p->num_individuals);
  assert(idx != NULL);

  /* the process for a time step is as follows:

  - selection based on fitness here (currently in driver)
  - mutation
  - transposition
  - 1pt recomb, 2pt recomb, gene recomb
  */
  
  for (i = 0; i < p->num_individuals; i++) {
    rv = (double)rand() / (double)RAND_MAX;
    /* flip and see if a mutation occurs in this individual */
    if (rv < p->p_m) {

      /* pick position */
      a = (int)(((double)rand() / (double)RAND_MAX)*
                (double)p->genome->individual_length);

      /* and generate a random value */
      rv = (double)rand() / (double)RAND_MAX;

      if (a>=p->genome->head_length) { /* in tail -- only terminals */
        if (p->genome->num_terminals > 1) {
          b = (int)(rv * (double)(p->genome->num_terminals));
          p->individuals[i][a] = 
            p->genome->terminals[b];
        }
      } else { /* in head : either terminals or functions */
        b = (int)(rv * (double)(p->genome->num_terminals + 
                                p->genome->num_functions));

        if (b > p->genome->num_functions) {
          p->individuals[i][a] = 
            p->genome->terminals[b - p->genome->num_functions];
        } else {
          p->individuals[i][a] = 
            p->genome->functions[b];
        }
      }
    }
  }

  for (i = 0; i < p->num_individuals; i++) {
    idx[i] = i;
  }

  i = p->num_individuals;
  for (j=0;j<i;j++) {
    rv = (double)rand() / (double)RAND_MAX;
    if (rv < p->p_x1) {
      rv = (double)rand() / (double)RAND_MAX;
      a = (int)((double)i * rv);
      b = a;
      while (b == a) {
        rv = (double)rand() / (double)RAND_MAX;
        b = (int)((double)i * rv);
      }
      aa = idx[a];
      bb = idx[b];
      idx[a] = idx[i-1];
      idx[b] = idx[i-2];
      i-=2;

      rv = (double)rand() / (double)RAND_MAX;
      rpos = (int)(rv * (double)p->genome->individual_length);
      a = p->genome->individual_length - rpos;

      tmp = (char *)malloc(sizeof(char) * a);

      memcpy(tmp,(p->individuals[bb])+rpos,a);
      memcpy((p->individuals[bb])+rpos,(p->individuals[aa])+rpos,a);
      memcpy((p->individuals[aa])+rpos,tmp,a);
      free(tmp);
    }
  }

  free(idx);

  return;
}
