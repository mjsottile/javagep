/**
 * FILE: regress.h
 *
 * author: matthew sottile (matt@cs.uoregon.edu)
 *
 * Licensed under the terms of the GNU Public Licence.  See LICENCE.GPL 
 * for details.
 */
#ifndef __REGRESS_H__
#define __REGRESS_H__

typedef enum { TERMINAL, FUNCTION } node_t;

/**
 * expression tree node for basic arith expressions
 */
typedef struct et_node {
  char tag;
  node_t ty;
  struct et_node *left, *right;
} et_node_t;

#ifdef __cplusplus
extern "C" {
#endif

  void destroy_et(et_node_t *n);
  void print_et(et_node_t *n);
  void print_expression(et_node_t *n);
  double eval_etnode(et_node_t *n, double x);
  et_node_t *r_express(char *ind, int len);

#ifdef __cplusplus
}
#endif

#endif /* __REGRESS_H__ */
