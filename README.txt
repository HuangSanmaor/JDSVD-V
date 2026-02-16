JDSVD_V, JDSVD_V_CP, JDSVD_V_HYBRID Documentation
*************************************************************

JDSVD_V: New Variant of Jacobi--Davidson SVD method.
***********************************************************

"JDSVD_V()" finds k singular values of the matrix A closest to a given target
and/or the corresponding left and right singular vectors using the standard
extraction-based new variant of the Jacobi-Davidson (JD) type SVD algorithm
such that the computed approximate partial SVD (T,U,V) of A satisfies
      ||AV-UT||_F^2+||A^TU-VT||_F^2 <= k||A||_2^2*tol^2,
where the diagonal elements of the diagonal matrix T save the approximate
singular values of A and the columns of U and V save the corresponding left
and right singular vectors, respectively.
 

JDSVD_V_CP: New Variant of Jacobi--Davidson SVD method Based on the
eigenproblem of the cross-product matrix
*****************************************************************************

JDSVD_V_CP() finds k singular values of the matrix A closest to a given target
and/or the corresponding left and right singular vectors using the standard
extraction-based new version of the Jacobi-Davidson (JD) type SVD algorithm
based on the equivalent eigenproblem of the cross product matrix A^TA such
that the computed approximate partial SVD (T,U,V) of A satisfies AV = UT and
      ||(A^TU-VT)T||_F <= k||A||_2*tol,
where the diagonal elements of the diagonal matrix T save the approximate
singular values of A and the columns of U and V save the corresponding
left and right singular vectors, respectively.
 

JDSVD_V_HYBRID: Two-Stage Hybrid Jacobi--Davidson SVD method--New Variant
***************************************************************************************

JDSVD_V_HYBRID() finds k singular values of the matrix A closest to a given target
and/or the corresponding left and right singular vectors using the standard
extraction-based new version of the Jacobi-Davidson (JD) type hybrid SVD algorithm
such that the computed approximate partial SVD (T,U,V) of A satisfies
      ||AV-UT||_F^2+||A^TU-VT||_F^2 <= k||A||_2^2*tol^2,
where the diagonal elements of the diagonal matrix T save the approximate singular
values of A and the columns of U and V save the corresponding left and right singular
vectors, respectively.



License Information
==========================================

JDSVD_V, JDSVD_V_CP, JDSVD_V_HYBRID are licensed under the 
3-clause license BSD.

   Copyright (c) 2018, Tsinghua University and Soochow University.
   All rights reserved.


Contact Information
==========================================

For reporting bugs or questions about functionality contact Zhongxiao
Jia by email, *jiazx* at *tsinghua.edu.cn* or Jinzhi Huang by email
*jzhuang21* at *suda.edu.cn*.


Support
==========================================

* National Science Foundation of China No. 12571404

* Youth Fund of the National Science Foundation of China No. 12301485

* Youth Program of the Natural Science Foundation of Jiangsu Province
  No. BK20220482

* Jiangsu Province Youth Science and Technology Talent Support Program
   No. JSTJ-2025-828


Pre-process
==========================================

Before using JDSVD_V, add the folder (without subfolders) "JDSVD-V-Master"
to the MATLAB search path list on the MATLAB homepage.

 
Implementation information--JDSVD_V
==========================================

T = JDSVD_V(A) returns 6 largest singular values of A.

T = JDSVD_V(A,K) returns K largest singular values of A.

T = JDSVD_V(A,K,SIGMA) returns K singular values of A depending on SIGMA:
       'largest'   - compute K largest singular values. This is the default.
       'smallest' - compute K smallest singular values.
       numeric   - compute K singular values nearest to SIGMA.

 T = JDSVD_V(A,K,SIGMA,NAME,VALUE) configures additional options specified
       by one or more name-value pair arguments. See the documentation below
       for more information.
 
 T = JDSVD_V(A,K,SIGMA,OPTIONS) alternatively configures the additional 
       options using a structure. See the documentation below for more 
       information.

[T,U,V] = JDSVD_V(A,...) computes the singular vectors as well. If A is  M-by-N
       K singular values are computed, then T is K-by-K diagonal, U and V are
       M-by-K and N-By-K orthonormal, respectively.

[T,U,V,RELRES] = JDSVD_V(A,...) also returns a matrix each column of which
       contains all the relative residual norms of the approximate singular triplets
       during computing the relavant desired exact singular triplet.

[T,U,V,RELRES,OUTIT] = JDSVD_V(A,...) also returns a vector of numbers of outer
       iterations uesd to compute each singular triplet.

[T,U,V,RELRES,OUTIT,INNIT] = JDSVD_V(A,...) also returns a vector of numbers of
       total inner iterations uesd to compute each singular triplet.

[T,U,V,RELRES,OUTIT,INNIT,INNITPERSTEP] = JDSVD_V(A,...) also returns a matrix
       each of whose column contains the numbers of inner iterations during each
       outer iteration when compute the relavant singular triplet.

[T,U,V,RELRES,OUTIT,INNIT,INNITPERSTEP,CPUTIME] = JDSVD_V(A,...) also returns
       a vector of CPU time used to compute each singular triplet.

[T,U,V,RELRES,OUTIT,INNIT,INNITPERSTEP,CPUTIME,FLAG] = JDSVD_V(A,...) also
       returns a convergence flag. If the method has converged, then FLAG = 0;
       If the maximun number of outer iterations have been used before convergence,
       then FLAG = 1.

[...] = JDSVD_V(AFUN,MN, ...) accepts function handle AFUN instead of the
       matrix A. AFUN(X,'notransp') must accept a vector input X and return the
       matrix-vector product A*X, while AFUN(U,'transp') must return A'*U. MN
       is a 1-by-2 row vector [M N] where M and N are the numbers of rows and
       columns of A respectively.


Implementation information--JDSVD_V_CP
==========================================

T = JDSVD_V_CP(A) returns 6 largest singular values of A.

T = JDSVD_V_CP(A,K) returns K largest singular values of A.

T = JDSVD_V_CP(A,K,SIGMA) returns K singular values of A depending on SIGMA:
       'largest' - compute K largest singular values. This is the default.
       'smallest' - compute K smallest singular values.
       numeric  - compute K singular values nearest to SIGMA.

T = JDSVD_V_CP(A,K,SIGMA,NAME,VALUE) configures additional options specified by
       one or more name-value pair arguments.  See the documentation below for more
       information.

T = JDSVD_V_CP(A,K,SIGMA,OPTIONS) alternatively configures the additional options 
      using a structure. See the documentation below for more information.

[T,U,V] = JDSVD_V_CP(A,...) computes the singular vectors as well. If A is
      M-by-N K singular values are computed, then T is K-by-K diagonal, U and V
      are M-by-K and N-By-K orthonormal, respectively.

[T,U,V,RELRES] = JDSVD_V_CP(A,...) also returns a vector of relative residual norms 
      of the computed approximate singular triplets.

[T,U,V,RELRES,OUTIT] = JDSVD_V_CP(A,...) also returns a vector of numbers
      of outer iterations uesd to compute each singular triplet.

[T,U,V,RELRES,OUTIT,INNIT] = JDSVD_V_CP(A,...) also returns a vector of
      numbers of total inner iterations uesd to compute each singular triplet.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME] = JDSVD_V_CP(A,...) also returns a
      vector of CPU time used to compute each singular triplet.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG] = JDSVD_V_CP(A,...) also returns
      a convergence flag. If the method has converged, then FLAG = 0; If the
      maximun number of outer iterations have been used before convergence, 
      then FLAG = 1.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV] = JDSVD_V_CP(A,...) also
      returns a vector of relative residual norms at each outer iteration.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC] = JDSVD_V_CP(A,...)
      also returns a vector of inner iterations at each outer iteration.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC,PERCENT] = JDSVD_V_CP(A,...)
      also returns the percentage of outer iterations during which the inner iterations
      achieve the desired convergece.

[...] = JDSVD_V_CP(AFUN,MN, ...) accepts function handle AFUN instead of
      the matrix A. AFUN(X,'notransp') must accept a vector input X and return
      the matrix-vector product A*X, while AFUN(U,'transp') must return A'*U.
      MN is a 1-by-2 row vector [M N] where M and N are the numbers of rows and
      columns of A respectively.


Implementation information--JDSVD_V_HYBRID
==========================================

T = JDSVD_V_HYBRID(A) returns 6 largest singular values of A.

T = JDSVD_V_HYBRID(A,K) returns K largest singular values of A.

T = JDSVD_V_HYBRID(A,K,SIGMA) returns K singular values of A depending on
SIGMA:
       'largest' - compute K largest singular values. This is the default.
      'smallest' - compute K smallest singular values.
        numeric  - compute K singular values nearest to SIGMA.

T = JDSVD_V_HYBRID(A,K,SIGMA,NAME,VALUE) configures additional options
      specified by one or more name-value pair arguments. See the documentation 
      below for more information.

T = JDSVD_V_HYBRID(A,K,SIGMA,OPTIONS) alternatively configures the additional
      options using a structure. See the documentation below for more information.

[T,U,V] = JDSVD_V_HYBRID(A,...) computes the singular vectors as well. If A is
      M-by-N K singular values are computed, then T is K-by-K diagonal, U
      and V are M-by-K and N-By-K orthonormal, respectively.

[T,U,V,RELRES] = JDSVD_V_HYBRID(A,...) also returns a vector of relative
      residual norms of the computed approximate singular triplets.

[T,U,V,RELRES,OUTIT] = JDSVD_V_HYBRID(A,...) also returns the total number
      of outer iterations.

[T,U,V,RELRES,OUTIT,INNIT] = JDSVD_V_HYBRID(A,...) also returns the total
      number of inner iterations.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME] = JDSVD_V_HYBRID(A,...) also returns
      the total CPU time in second.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG] = JDSVD_V_HYBRID(A,...) also returns
      a convergence flag. If the method has converged, then FLAG = 0; If the
      maximun number of outer iterations have been used before convergence, 
      then FLAG = 1.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV] = JDSVD_V_HYBRID(A,...) also
      returns a vector of relative residual norms at each outer iteration.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC] = JDSVD_V_HYBRID(A,...)
      also returns a vector of inner iterations at each outer iteration.

[T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC,PERCENT] = JDSVD_V_HYBRID(A,...)
      also returns the percentage of outer iterations during which the inner iterations
      achieve the desired convergece.

[...] = JDSVD_V_HYBRID(AFUN,MN, ...) accepts function handle AFUN instead
      of the matrix A. AFUN(X,'notransp') must accept a vector input X and return
      the matrix-vector product A*X, while AFUN(U,'transp') must return A'*U.
      MN is a 1-by-2 row vector [M N] where M and N are the numbers of rows and
      columns of A respectively.


Description for the parameters specified by the name value pairs. 
==========================================

                                        NMAE - DESCRIPTION
                                 'Tolerance' - Convergence tolerance
                           'MaxIterations' - Maximum number of iterations
          'MaxSubspaceDimension' - Maximum size of subspaces
          'MinSubspaceDimension' - Minimum size of subspaces
                  'SwitchingTolerance' - Switching tolerance for correction equations
                         'InnerTolerance' - Parameter for inner stopping tolerance
                         'LeftStartVector' - Left starting vector
                       'RightStartVector' - Right starting vector
   'InnerPreconditionTolerance1' - 1st parameter for inner preconditioning
   'InnerPreconditionTolerance1' - 2nd parameter for inner preconditioning
                                      'Display' - Display diagnostic messages


Description for the parameters contained in the structure OPTS.
==========================================

PARAMETER           DESCRIPTION

  OPTS.TOL            Convergence tolerance, same as the name 'Tolerance'.
                              Convergence is determined when
                                    ||[A^Tu-θv;Av-θu|| <= ||A||·TOL,
                              where (θ,u,v) is the current approximate singular
                      	      triplet of A and ||A|| is the 1-norms of A is A is
                              given as a matrix or estimate 2-norm of A if the
                              function handle Afun is given.
                              DEFAULT VALUE    TOL = 1e-8

  OPTS.MAXIT        Maximum number of outer iterations, same as the
                              parameter 'MaxIterations'.
                              DEFAULT VALUE    MAXIT = N

  OPTS.MAXSIZE    Maximum size of searching subspaces, same as
                              'MaxSubspaceDimension'.
                              DEFAULT VALUE    MAXSIZE = 30

  OPTS.MINSIZE     Minimum size of searching subspaces, same as
                              'MinSubspaceDimension'.
                              DEFAULT VALUE    MINSIZE = 3

  OPTS.FIXTOL       Switching tolerance for the inner correction euqations,
                              same as the parameter 'SwitchingTolerance'.
                              DEFAULT VALUE    FIXTOL = 1e-4

  OPTS.INNTOL      Accuracy requirement for the apprixmate solution of
                              the inner correction equations, same as the
                              'SwitchingTolerance'.
                              DEFAULT VALUE    FIXTOL = 1e-4

  OPTS.U0             Left starting vector, same as the parameter
                             'LeftStartVector'.
                             DEFAULT VALUE  U0 = randn(M,1);

  OPTS.V0             Right starting vector, same as the parameter
                             'RightStartVector'.
                             DEFAULT VALUE  V0 = randn(N,1);

  OPTS.INNPRECTOL1    The first inner preconditioning parameter, same as
                                       the parameter 'InnerPreconditionTolerance1'.
                             DEFAULT VALUE  INNPRECTOL1 = 0.05;

  OPTS.INNPRECTOL2    The second inner preconditioning parameter, same as
                                       the parameter 'InnerPreconditionTolerance2'.
                             DEFAULT VALUE  INNPRECTOL1 = 0.01;

  OPTS.DISPS        Indicates if K approximate singular values are to be
                             displayed during the computation. Set DISPS > 1 to
                             display the values at each outer iteration or when
                             the algorithm is stopped. If 0 < DISPS <= 1, then the
                             results are only display after the overall convergence.
                             DEFAULT VALUE   DISPS = 0

 


Test Examples
==========================================

M = 10000;
sigA = [0.8:0.0001:0.9979,1.0001:0.0001:1.002];
N = size(sigA,2);
dense = 0.2;
A = sprand(M,N,dense,sigA);

k = 20;
target = 1;
outtol = 1e-8;
maxit = min(M,N);
maxdim = 30;
mindim = 3;
fixtol = 0;
inntol = 1e-4;
u0 = ones(M,1);
v0 = ones(N,1);

% % Tests the JDSVD_V algorithm with matrices;
[Tp,Up,Vp,RELRESp,OUTITp,INNITp,INNITPERSTEPp,CPUTIMEp,FLAGp,PERCp] = ...
    JDSVD_V(A,k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,...
    'InnerPreconditionTolerance1',innprectol1,...
    'InnerPreconditionTolerance2',innprectol2,'Display',display);
[To,Uo,Vo,RELRESo,OUTITo,INNITo,INNITPERSTEPo,CPUTIMEo,FLAGo,PERCo] =...
    JDSVD_V(A,k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,'Display',display,...
    'InnerPreconditionTolerance1',0,'InnerPreconditionTolerance2',0);
% % Tests the JDSVD_V algorithm with function handles and/or options;
jdsvd_v(A,k,target,options);
jdsvd_v(afun,[M N],k,target,options);


% Tests the JDSVD_V_CP algorithm with matrices;
[Tp,Up,Vp,RELRESp,OUTITp,INNITp,CPUTIMEp,FLAGp,RELRESVEVp,INNERVECp,PERCp] = ...
    JDSVD_V_CP(afun,[M N],k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,...
    'InnerPreconditionTolerance1',innprectol1,...
    'InnerPreconditionTolerance2',innprectol2,'Display',display); 
[To,Uo,Vo,RELRESo,OUTITo,INNITo,CPUTIMEo,FLAGo,RELRESVEVo,INNERVECo,PERCo] =...
    JDSVD_V_CP(A,k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,'Display',display,...
    'InnerPreconditionTolerance1',0,'InnerPreconditionTolerance2',0);
% % Tests the JDSVD_V_CP algorithm with function handles and/or options;
jdsvd_v_cp(A,k,target,options);
jdsvd_v_cp(afun,[M N],k,target,options);


% % Tests the JDSVD_V_HYBRID algorithm with matrices;
[Tp,Up,Vp,RELRESp,OUTITp,INNITp,CPUTIMEp,FLAGp,RELRESVEVp,INNERVECp,PERCp] = ...
    JDSVD_V_HYBRID(A,k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,...
    'InnerPreconditionTolerance1',innprectol1,...
    'InnerPreconditionTolerance2',innprectol2,'Display',display);
[To,Uo,Vo,RELRESo,OUTITo,INNITo,INNITPERSTEPo,CPUTIMEo,FLAGo,PERCo] =...
    JDSVD_V_HYBRID(A,k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,'Display',display,...
    'InnerPreconditionTolerance1',0,'InnerPreconditionTolerance2',0);
% % Tests the IPJDSVD_CP algorithm with function handles and/or options;
jdsvd_v_hybrid(A,k,target,options);
jdsvd_v_hybrid(afun,[M N],k,target,options);


% MatrixFunction turns a matrix to a function handle
% For experimental purpose only
function varargout = MatrixFunction(C)
% Creat the function handle with a given matrix C.

varargout{1} = @matrixfun;
    function y = matrixfun(x,transpornot)
        if strcmp(transpornot,'notransp')
            y = C*x;
        elseif strcmp(transpornot,'transp')
            y = C'*x;
        else
            error('wrong input for transpornot');
        end
    end
end
See also: *testrand.m*.


REFERENCES:
 [1] Jinzhi Huang and Zhongxiao Jia, On inner iterations of Jacobi-Davidson
      type methods for large SVD computations, SIAM J. SCI. COMPUT., 41,3
      (2019), pp. A1574–A1603.

 [2] Jinzhi Huang and Zhongxiao Jia, New Jacobi--Davidson type methods for
      the large SVD computations, (2026), 25 pages.
 