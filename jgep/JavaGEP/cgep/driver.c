/**
 * FILE: driver.c
 *
 * author: matthew sottile (matt@cs.uoregon.edu)
 *
 * Licensed under the terms of the GNU Public Licence.  See LICENCE.GPL 
 * for details.
 */
#include "gep.h"
#include "regress.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char terminals[] = {'a'};
char functions[] = {'*','-','+','/'};

double inputs[] = \
  {-2.0,-1.9,-1.8,-1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, \
   -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, \
   0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};

double outputs[] = \
  {0.110115, -0.378666, -0.0884169, 0.0496072, -0.317578,-0.756975,-0.853041, \
   -0.581914, -0.15522, 0.222213, 0.454649, 0.540203, 0.520615, 0.443212, \
   0.344088, 0.245474, 0.158992, 0.0898458, 0.0399881, 0.00999983, 0,   \
   0.00999983, 0.0399881, 0.0898458, 0.158992, 0.245474, 0.344088, 0.443212, \
   0.520615, 0.540203, 0.454649, 0.222213, -0.15522, -0.581914, -0.853041, \
   -0.756975, -0.317578, 0.0496072, -0.0884169, -0.378666, 0.110115};

typedef struct flist {
  double fv;
  int idx;
} flist_t;

int flcompare(const void *aa, const void *bb) {
  flist_t *a = (flist_t *)aa;
  flist_t *b = (flist_t *)bb;
  if (a->fv > b->fv) return -1;
  if (a->fv == b->fv) return 0;
  return 1;
}

#define POPSIZE 200
#define TESTCASES 41

int main(int argc, char **argv) {
  int generation = 0;
  gep_genome_t *genome;
  gep_population_t *pop;
  int i,j;
  et_node_t *r;
  double d;
  flist_t *fl;
  double fitcur, fitprev;
  int keepgoing = 1;
  int cutoff = (int)((double)POPSIZE * 0.5);
  char **tmp;

  if (argc > 1) {
    srand(atoi(argv[1]));
  }

  fl = (flist_t *)malloc(sizeof(flist_t)*POPSIZE);
  assert(fl != NULL);
  genome = gep_create_genome(1,terminals,4,functions,150,2);
  pop = gep_create_population(genome,POPSIZE,0.34,0.4,0.3,4000.0);

  fitprev = fitcur = -100000.0;

  while (keepgoing==1) {
    
    for(i=0;i<POPSIZE;i++) {
      r = r_express(pop->individuals[i],pop->genome->individual_length);
      pop->fitnesses[i] = pop->max_fitness;
      
      for (j=0;j<TESTCASES;j++) {
        d = eval_etnode(r,inputs[j]);
        pop->fitnesses[i] -= ((outputs[j] - d)*100)*((outputs[j] - d)*100);
      }
   
      destroy_et(r);
   
      fl[i].fv = pop->fitnesses[i];
      fl[i].idx = i;
    }
    
    qsort(fl,POPSIZE,sizeof(flist_t),flcompare);
    
    if (fl[0].fv == 0.0) {
      keepgoing = 0;
    }

    /* best always continues */
    memcpy(pop->selectionbuffer[0],
           pop->individuals[fl[0].idx],
           pop->genome->individual_length * sizeof(char));
    
    for (i = 1; i < POPSIZE; i++) {
      j = (int)(((double)rand() / (double)RAND_MAX)*(double)cutoff);
      memcpy(pop->selectionbuffer[i],
             pop->individuals[fl[j].idx],
             pop->genome->individual_length * sizeof(char));
    }
    
    tmp = pop->individuals;
    pop->individuals = pop->selectionbuffer;
    pop->selectionbuffer = tmp;

    fitcur = fl[0].fv;

    if (fitcur > fitprev) {
      printf("GENERATION: %d\n",generation); 
      printf("best individual so far (%f): \n",fitcur);
      for (i=0;i<pop->genome->individual_length;i++)
        printf("%c",pop->individuals[fl[0].idx][i]);
      printf("\n");

      fitprev = fitcur;
    }

    generation++;
    if (generation == 2000) keepgoing = 0;
    gep_next_generation(pop);
  }

  printf("best individual: \n");
  for (i=0;i<pop->genome->individual_length;i++)
    printf("%c",pop->individuals[fl[0].idx][i]);
  printf("\n");

  r = r_express(pop->individuals[fl[0].idx],
                pop->genome->individual_length);
  print_expression(r);
  printf("\n");
  destroy_et(r);


  gep_destroy_population(pop);
  gep_destroy_genome(genome);
  
  return(EXIT_SUCCESS);
}
