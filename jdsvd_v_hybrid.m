function varargout = jdsvd_v_hybrid(A,varargin)
% JDSVD_V_HYBRID finds k singular values of the matrix A closest to a given 
% target and/or the corresponding left and right singular vectors using the 
% standard extraction-based new version of the Jacobi-Davidson (JD) type 
% hybrid SVD algorithm such that the computed approximate partial SVD (T,U,V) 
% of A satisfies
%       ||AV-UT||_F^2+||A^TU-VT||_F^2 <= k||A||_2^2*tol^2,
% where the diagonal elements of the diagonal matrix T save the approximate
% singular values of A and the columns of U and V save the corresponding
% left and right singular vectors, respectively.
%
% T = JDSVD_V_HYBRID(A) returns 6 largest singular values of A.
% 
% T = JDSVD_V_HYBRID(A,K) returns K largest singular values of A.
% 
% T = JDSVD_V_HYBRID(A,K,SIGMA) returns K singular values of A depending on 
% SIGMA:
% 
%        'largest' - compute K largest singular values. This is the default.
%       'smallest' - compute K smallest singular values.
%         numeric  - compute K singular values nearest to SIGMA.
% 
% T = JDSVD_V_HYBRID(A,K,SIGMA,NAME,VALUE) configures additional options 
% specified by one or more name-value pair arguments:
% 
%                     'Tolerance' - Convergence tolerance
%                 'MaxIterations' - Maximum number of iterations
%          'MaxSubspaceDimension' - Maximum size of subspaces
%          'MinSubspaceDimension' - Minimum size of subspaces
%            'SwitchingTolerance' - Switching tolerance for correction equations
%                'InnerTolerance' - Parameter for inner stopping tolerance
%               'LeftStartVector' - Left starting vector
%              'RightStartVector' - Right starting vector
%   'InnerPreconditionTolerance1' - 1st parameter for inner preconditioning
%   'InnerPreconditionTolerance1' - 2nd parameter for inner preconditioning
%                       'Display' - Display diagnostic messages
% 
% T = JDSVD_V_HYBRID(A,K,SIGMA,OPTIONS) alternatively configures the additional
% options using a structure. See the documentation below for more information.
% 
% [T,U,V] = JDSVD_V_HYBRID(A,...) computes the singular vectors as well. If 
% A is M-by-N K singular values are computed, then T is K-by-K diagonal, U 
% and V are M-by-K and N-By-K orthonormal, respectively.
% 
% [T,U,V,RELRES] = JDSVD_V_HYBRID(A,...) also returns a vector of relative 
% residual norms of the computed approximate singular triplets. 
% 
% [T,U,V,RELRES,OUTIT] = JDSVD_V_HYBRID(A,...) also returns the total number
% of outer iterations.
% 
% [T,U,V,RELRES,OUTIT,INNIT] = JDSVD_V_HYBRID(A,...) also returns the total 
% number of inner iterations.
% 
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME] = JDSVD_V_HYBRID(A,...) also returns 
% the total CPU time in second. 
% 
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG] = JDSVD_V_HYBRID(A,...) also returns 
% a convergence flag. If the method has converged, then FLAG = 0; If the 
% maximun number of outer iterations have been used before convergence, then FLAG = 1. 
% 
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV] = JDSVD_V_HYBRID(A,...) also 
% returns a vector of relative residual norms at each outer iteration.
% 
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC] = JDSVD_V_HYBRID(A,...)
% also returns a vector of inner iterations at each outer iteration.
% 
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC,PERCENT] = JDSVD_V_HYBRID(A,...)
% also returns the percentage of outer iterations during which the inner iterations
% achieve the desired convergece. 
% 
% [...] = JDSVD_V_HYBRID(AFUN,MN, ...) accepts function handle AFUN instead 
% of the matrix A. AFUN(X,'notransp') must accept a vector input X and return
% the matrix-vector product A*X, while AFUN(U,'transp') must return A'*U.
% MN is a 1-by-2 row vector [M N] where M and N are the numbers of rows and
% columns of A respectively.
%
% Description for the parameters contained in the structure OPTS.
%
%  PARAMETER           DESCRIPTION
%
%  OPTS.TOL            Convergence tolerance, same as the name 'Tolerance'.
%                      Convergence is determined when
%                               ||[A^Tu-θv;Av-θu|| <= ||A||·TOL,
%                      where (θ,u,v) is the current approximate singular
%                      triplet of A and ||A|| is the 2-norm of A.
%                      DEFAULT VALUE    TOL = 1e-8.
%
%  OPTS.MAXIT          Maximum number of outer iterations, same as the
%                      parameter 'MaxIterations'.
%                      DEFAULT VALUE    MAXIT = N.
%
%  OPTS.MAXSIZE        Maximum size of searching subspaces, same as
%                      'MaxSubspaceDimension'.
%                      DEFAULT VALUE    MAXSIZE = 30.
%
%  OPTS.MINSIZE        Minimum size of searching subspaces, same as
%                      'MinSubspaceDimension'.
%                      DEFAULT VALUE    MINSIZE = 3.
%
%  OPTS.FIXTOL         Switching tolerance for the inner correction euqations,
%                      same as the parameter 'SwitchingTolerance'.
%                      DEFAULT VALUE    FIXTOL = 1e-4.
%
%  OPTS.INNTOL         Accuracy requirement for the apprixmate solution of
%                      the inner correction equations, same as the
%                      'SwitchingTolerance'.
%                      DEFAULT VALUE    FIXTOL = 1e-4.
%
%  OPTS.U0             Left starting vector, same as the parameter
%                      'LeftStartVector'.
%                      DEFAULT VALUE  U0 = randn(M,1).
%
%  OPTS.V0             Right starting vector, same as the parameter
%                      'RightStartVector'.
%                      DEFAULT VALUE  V0 = randn(N,1).
%
%  OPTS.INNPRECTOL1    The first inner preconditioning parameter, same as
%                      the parameter 'InnerPreconditionTolerance1'.
%                      DEFAULT VALUE  INNPRECTOL1 = 0.05.
%
%  OPTS.INNPRECTOL2    The second inner preconditioning parameter, same as
%                      the parameter 'InnerPreconditionTolerance2'.
%                      DEFAULT VALUE  INNPRECTOL1 = 0.01.
%
%  OPTS.DISPS          Indicates if K approximate singular values are to be
%                      displayed during the computation. Set DISPS > 1 to
%                      display the values at each outer iteration or when
%                      the algorithm is stopped. If 0 < DISPS <= 1, then the
%                      results are only display after the overall convergence.
%                      DEFAULT VALUE   DISPS = 0.
%
% REFERENCES:
% [1] Jinzhi Huang and Zhongxiao Jia, On inner iterations of Jacobi-Davidson
%     type methods for large SVD computations, SIAM J. SCI. COMPUT., 41,3
%     (2019), pp. A1574–A1603.
% [2] Jinzhi Huang and Zhongxiao Jia, New Jacobi--Davidson type methods for 
%     the large SVD computations, (2026), 25 pages.

% Check the number of outputs.
if nargout == 2  || nargout >= 12
    error('Incorrect number of output arguments.');
end
% Initialize a dedicated randstream, to make output reproducible.
randStr = RandStream('dsfmt19937','Seed',0);
% Get inputs and check the items.
VarIn = varargin;
[A,M,N,k,target,u0,v0,Options] = tools.chekInputs(A,VarIn,randStr);
% Reset the stream.
reset(randStr,1);
% Use the function handles if given, otherwise build them from the matrices.
if isa(A,'function_handle')
    Afun = A;  
elseif ismatrix(A)
    Afun = tools.MatrixFunction(A);  
end
if M < N
    AfunT = tools.FunctionT(Afun);
end
if ~isfield(Options,'normAest')
    normAest = svds(Afun,[M,N],1,"largest","Tolerance",1e-2,"MaxIterations",30);
    Options.normAest = normAest;
end
% If target = 0, replace the left and right initial vectors u0 and v0 as
% A^Tu0 and Av0, respectively, in order to discard the information on the
% zero singular value of A and nontrivial singular triplets corresponding
% to nonzero singular values of A are computed. 
if (issparse(target)&&(target == 0)) || strcmp(target,'smallest') 
    if ~isfield(Options,'u0')
        ut = u0;
        u0 = Afun(v0,'notransp');
    end
    if ~isfield(Options,'v0')
        if ~isfield(Options,'u0')
            v0 = Afun(ut,'transp');
        else
            v0 = Afun(u0,'transp');
        end
    end
end 
%--------------------------------%
% BEGIN: CHECK FOR EMPTY MATRIX. %
%--------------------------------%
if M==0 || N==0
    % Nothing need to been done so no singular values to be displayed.
    if nargout == 0
        fprintf('\nAlgorithm: the Standard JDSVD_V_HYBRID algorithm;\n');
        fprintf('Nothing need to be done since the given matrices is empty');
        fprintf('Total number of outer iterations = %d;\n',0);
        fprintf('Total number of inner iterations = %d;\n',0);
        fprintf('Total CPU time in seconds = %s;\n\n',0);
    end
    % Empty vector for computed singular values.
    if nargout == 1
        varargout{1} = [];
    end
    % Empty matrices for computed partial SVD.
    if nargout > 1
        varargout{1} = [];
        varargout{2} = [];
        varargout{3} = [];
    end
    % Empty matrix for the relative residual norms.
    if nargout >= 4
        varargout{4} = [];
    end
    % Empty vector for numbers of outer iterations.
    if nargout >= 5
        varargout{5}= [];
    end
    % Empty vector for numbers of inner iterations.
    if nargout >= 6
        varargout{6} = [];
    end
    % Emtpy vector for CPU time.
    if nargout >= 7
        varargout{7} = 0;
    end
    % The flag is set as zero for empty matrix.
    if nargout >= 8
        varargout{8} = 0;
    end
    % Empty vector for relative residual norms of each step. 
    if nargout >= 9
        varargout{9} = [];
    end
    % Empty vector for numbers of inner iterations during each step. 
    if nargout >= 10
        varargout{10} = [];
    end
    % The percentage of inner convergence is set as 1 for empty matrix.
    if nargout >= 11
        varargout{11} = 1;
    end 
    return;
end
%------------------------------%
% END: CHECK FOR EMPTY MATRIX. %
%------------------------------%
% Converged partial SVD of A that is desired by the users.
T_conve = zeros(k,k);
U_conve = zeros(M,k);
V_conve = zeros(N,k);
% Computational results recorded to illstrate the performence of the algorithm.
RelaResidual = inf(k,1);  
% % Implement the Standard JDSVD_V_CP method based on the eigenproblem of the 
% cross product matrix A^TA at the first stage, with a relatively large 
% outer stopping tolerance 
if M >= N
    Options_cp = Options;
    Options_cp.u0 = u0;
    Options_cp.v0 = v0; 
    [Tcp,Ucp,Vcp,RELREScp,OUTITcp,INNITcp,CPUTimecp,FLAGcp,RELRESVECcp,...
        INNERVECcp,PERCENTcp] = jdsvd_v_cp(Afun,[M N],k,target,Options_cp);
else
    Options_cp = Options;
    Options_cp.u0 = v0;
    Options_cp.v0 = u0;
    [Tcp,Vcp,Ucp,RELREScp,OUTITcp,INNITcp,CPUTimecp,FLAGcp,RELRESVECcp,...
        INNERVECcp,PERCENTcp] = jdsvd_v_cp(AfunT,[N M],k,target,Options_cp);
end
%k_convecp = size(Tcp,1) - FLAGcp;
converges = find(RELREScp <= Options.tol); 
k_convecp = length(converges);
U_conve(:,1:k_convecp) = Ucp(:,converges); 
V_conve(:,1:k_convecp) = Vcp(:,converges); 
T_conve(1:k_convecp,1:k_convecp) = Tcp(converges,converges);
RelaResidual(1:k_convecp) = RELREScp(converges); 
% Update normAest if the largest singular values are computed
if strcmp(target,'largest') && k_convecp >= 1 && k_convecp < k
    Options.normAest = Tcp(1,1);
end 
% % If the total outer stopping tolerance is small, use the results
% computed at the first stage as the initial left and right searching
% subspaces for the second stage and perform the Standard JDSVD_V method
% based on the eigenproblem of the augmented matrix [0 A;A^T 0] yet with
% the left and right searching subspaces dealt separately.
if k_convecp < k
    unconvergs = find(RELREScp > Options.tol); 
    Options_am = Options;
    Options_am.maxit = Options.maxit - sum(OUTITcp);
    Options_am.U_conve = U_conve(:,1:k_convecp);
    Options_am.V_conve = V_conve(:,1:k_convecp);
    Options_am.u0 = Ucp(:,unconvergs);
    Options_am.v0 = Vcp(:,unconvergs); 
    [Tam,Uam,Vam,RELRESam,OUTITam,INNITam,CPUTimeam,FLAGam,RELRESVECam,...
        INNERVECam,PERCENTam] = jdsvd_v(Afun,[M N],k-k_convecp,target,Options_am); 
    k_conveam = size(Tam,1);
    U_conve(:,k_convecp+(1:k_conveam)) = Uam;
    V_conve(:,k_convecp+(1:k_conveam)) = Vam;
    T_conve(k_convecp+(1:k_conveam),k_convecp+(1:k_conveam)) = Tam(1:k_conveam,1:k_conveam);
    RelaResidual(k_convecp+(1:k_conveam)) = RELRESam;  
else
    k_conveam = 0;
    OUTITam = 0;
    INNITam = 0;
    CPUTimeam = 0;
    RELRESVECam = [];
    INNERVECam = [];
    FLAGam = 0;
end
% Handle the outputs. 
OuterIter = [sum(OUTITcp) sum(OUTITam)];
InnerIter = [sum(INNITcp) sum(INNITam)];
CPUTimese = [sum(CPUTimecp) sum(CPUTimeam)];
RELRESVEC = [RELRESVECcp;RELRESVECam];
INNERITERVEC = [INNERVECcp;INNERVECam]; 
FLAGofJDSVDV = FLAGam;
if sum(OUTITam) == 0
    PercentConv = PERCENTcp;
else
    PercentConv = (PERCENTcp*sum(OUTITcp)+...
        PERCENTam*sum(OUTITam))/sum(OuterIter);
end
k_result = k_convecp + k_conveam; 
T_conve = T_conve(1:k_result,1:k_result);
U_conve = U_conve(:,1:k_result);
V_conve = V_conve(:,1:k_result); 
RelaResidual = RelaResidual(1:k_result);
Isort = tools.SortSigmas(diag(T_conve),target);
T_conve = T_conve(Isort,Isort);
U_conve = U_conve(:,Isort);
V_conve = V_conve(:,Isort); 
% Compute the results when needed.
if nargout == 0 || Options.disp > 0 
    Results=zeros(k_result,3);
    Results(:,1) = (1:k_result)';
    Results(:,2) = diag(T_conve);
    Results(:,3) = RelaResidual;  
    k_conve = sum(Results(:,3) <= Options.tol);
end
%------------------------%
% BEGIN: OUTPUT RESULTS. %
%------------------------%
% Output option I: Display the singular values. 
if nargout == 0 || Options.disp > 0 
    fprintf('Number of converged singular triplets: %d;\n',k_conve);
    fprintf('  i:     SingularValues;     RelResNorms;; \n')
    fprintf('%3d:     %14.5e;     %11.5e;\n',Results');
    fprintf('Total number of outer iterations = %8d;%8d;\n',OuterIter);
    fprintf('Total number of inner iterations = %8d;%8d;\n',InnerIter);
    fprintf('Total CPU time in seconds = %8.2e;  %8.2e;\n',CPUTimese);
    fprintf('Percentage of converged outer iterations = %6.2f%%;\n',PercentConv*100);
end
% Output option II: Output the singualr values.
if nargout == 1
    varargout{1} = diag(T_conve);
end
% Output option III: Output the partial SVD.
if nargout > 1
    varargout{1} = T_conve;
    varargout{2} = U_conve;
    varargout{3} = V_conve;
end
% Output option IV: Output the relative residual norm of the approximate 
% singular triplets.
if nargout >= 4
    varargout{4} = RelaResidual;
end
% Output option V: Output the numbers of outer iterations the algorithm 
% takes to computed each approximate singular triplet.
if nargout >= 5
    varargout{5} = OuterIter;
end
% Output option VI: Output the numbers of total inner iterations the 
% algorithm takes to computed each approximate singular triplet.
if nargout >= 6
    varargout{6} = InnerIter;
end
% Output option VII: Output the total CPU time the algorithm takes to 
% compute each approximate singular triplet of A.
if nargout >= 7
    varargout{7} = CPUTimese;
end
% Output option VIII: Oupput the flag: flag=0 if all the desired singular 
% triplets have been successfully computed; flag=1 if the maximum number 
% of outer iterations have been used.
if nargout >= 8
    varargout{8} = FLAGofJDSVDV;
end
% Output Options IX: Output the relative residual norms of the approximate 
% singular triplet during each outer iteration.
if nargout >= 9
    varargout{9} = RELRESVEC; 
end
% Output option X: Output the numbers of inner iterations during each 
% outer iteration.
if nargout >= 10
    varargout{10} = INNERITERVEC;
end
% Output Options XI: Oupput the percentage of total outer iterations during
% which the inner iterations have converged to the desired inner accuracy.
if nargout >= 11
    varargout{11} = PercentConv;
end
%----------------------%
% END: OUTPUT RESULTS. %
%----------------------%
end    


function outitnew = spreading(outit,converges,k)
outitnew = zeros(k,1);
t = setdiff(1:k,converges);
outitnew(t) = [outit;zeros(length(t)-length(outit),1)]; 
end