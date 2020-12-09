/*
 * pwsim.h
 *
 */

#ifndef _PWSIM_H_
#define _PWSIM_H_

extern void
pwdot(
  const int *M_,
  const int *N_,
  const double *X,
  double *C           // [N(N+1)/2] upper triangular
  )
;

extern void
pwscale(
  const int *N_,
  double *C
  )
;

#endif // _PWSIM_H_
