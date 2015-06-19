/**
 * FILE: regress.c
 *
 * author: matthew sottile (matt@cs.uoregon.edu)
 *
 * Licensed under the terms of the GNU Public Licence.  See LICENCE.GPL 
 * for details.
 */
/**
 * basic regression example to test gep code.
 */
#include "gep.h"
#include "regress.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>

/**
 * print the expression tree, explicit indentation
 */
void __print_et(et_node_t *n, int indent) {
  int i;

  if (n == NULL) return;
  for (i=0;i<indent;i++) { printf("-"); }

  printf("%c\n",n->tag);
  __print_et(n->left,indent+1);
  __print_et(n->right,indent+1);
}

/**
 * print the expression tree, wrapper around actual function.
 */
void print_et(et_node_t *n) {
  __print_et(n,0);
}

void print_expression(et_node_t *n) {
  if (n==NULL) return;
  if (n->tag != 'a' && n->tag != '\0') {
    printf("(");
    print_expression(n->left);
    printf("%c ",n->tag);
    print_expression(n->right);
    printf(")");
  } else printf("a");
}

typedef enum { DIVBYZERO=1, NULLPTR=2 } eval_error_t;

/**
 * evaluate an expression tree containing binary functions and
 * one variable.
 */
double eval_etnode(et_node_t *n, double x) {
  double l,r;
  double retval;

  /* return x if a terminal */
  if (n->ty == TERMINAL) return x;

  /* if left is non-null, evaluate it. */
  if (n->left != NULL) {
    l = eval_etnode(n->left,x);
  } else {
    /*    longjmp(env,NULLPTR); */
  }

  if (n->right != NULL) {
    r = eval_etnode(n->right,x);
  } else {
    /* longjmp(env,NULLPTR); */
  }

  switch (n->tag) {
  case '+':
    retval = l+r;
    break;
  case '-':
    retval = l-r;
    break;
  case '*':
    retval = l*r;
    break;
  case '/':
    if (r == 0) {
      retval = -1000000.0;
    } else
      retval = l/r;
    break;
  case '^':
    if (r < 0) retval = -1000000.0;
    else
      retval = pow(l,r);
    break;
  default:
    break;
  }
  
  return retval;
}

/**
 * 'constructor' for expression tree nodes
 */
et_node_t *new_etnode(char c, node_t ty) {
  et_node_t *n;

  n = (et_node_t *)malloc(sizeof(et_node_t));
  assert(n != NULL);
  
  n->left = n->right = NULL;
  n->tag = c;
  n->ty = ty;

  return n;
}

/**
 * 'destructor' for expression tree nodes
 */
void destroy_et(et_node_t *root) {
  if (root == NULL) return;

  destroy_et(root->left);
  destroy_et(root->right);
  free(root);
}

/**
 * basic linked list data structure used for unflattening the tree.
 */
struct linklist {
  et_node_t *node; 
  struct linklist *next;
};

/**
 * Create a new linked list node with the given et_node
 */
struct linklist *new_ll(et_node_t *node) {
  struct linklist *ll;
  
  ll = (struct linklist *)malloc(sizeof(struct linklist));
  assert(ll != NULL);
  ll->next = NULL;
  ll->node = node;

  return ll;
}

/**
 * destructor for the linked list.
 */
void destroy_ll(struct linklist *ll) {
  if (ll == NULL) return;
  
  destroy_ll(ll->next);
  free(ll);
}

/**
 * expression function for regression expression trees.
 */
et_node_t *r_express(char *ind, int len) {
  int pos, count;
  struct linklist *ll1_head, *ll1_cur;
  struct linklist *ll2_head, *ll2_cur;
  et_node_t *n1, *root;

  /* null out the lists */
  ll1_head = ll1_cur = NULL;
  ll2_head = ll2_cur = NULL;

  /* don't bother if the chromosome represents a constant. */
  if (ind[0] == 'a') return new_etnode(ind[0],TERMINAL);

  ll1_cur = ll1_head = new_ll(new_etnode(ind[0],FUNCTION));
  root = ll1_cur->node;
  count = 2;
  pos = 1;

  while (count > 0 && pos < len) {

    while (ll1_cur != NULL && count > 0) {
      if (ind[pos] == 'a') {
        n1 = new_etnode(ind[pos],TERMINAL); 
        count--;
      } else {
        n1 = new_etnode(ind[pos],FUNCTION);
        count += 1;
      }

      pos++;

      if (ll1_cur->node->left == NULL) 
        ll1_cur->node->left = n1;
      else {
        ll1_cur->node->right = n1;
        ll1_cur = ll1_cur->next;
      }
      
      if (n1->ty != TERMINAL) {
        if (ll2_cur == NULL) {
          ll2_head = ll2_cur = new_ll(n1);
        } else {
          ll2_cur->next = new_ll(n1);
          ll2_cur = ll2_cur->next;
        }
      } 
    }

    destroy_ll(ll1_head);
    ll1_head = ll1_cur = ll2_head;
    ll2_head = ll2_cur = NULL;
  }

  if (pos == len && count > 0) {
    printf("ERROR!\n");
    return NULL;
  }

  destroy_ll(ll1_head);

  return root;
}

