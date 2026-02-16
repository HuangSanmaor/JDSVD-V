function varargout = jdsvd_v(A,varargin)
% JDSVD_V finds k singular values of the matrix A closest to a given target 
% and/or the corresponding left and right singular vectors using the standard 
% extraction-based new version of the Jacobi-Davidson (JD) type SVD algorithm 
% such that the computed approximate partial SVD (T,U,V) of A satisfies
%       ||AV-UT||_F^2+||A^TU-VT||_F^2 <= k||A||_2^2*tol^2,
% where the diagonal elements of the diagonal matrix T save the approximate
% singular values of A and the columns of U and V save the corresponding
% left and right singular vectors, respectively.
%
% T = JDSVD_V(A) returns 6 largest singular values of A.
%
% T = JDSVD_V(A,K) returns K largest singular values of A.
%
% T = JDSVD_V(A,K,SIGMA) returns K singular values of A depending on SIGMA:
%
%        'largest' - compute K largest singular values. This is the default.
%       'smallest' - compute K smallest singular values.
%         numeric  - compute K singular values nearest to SIGMA.
%
% T = JDSVD_V(A,K,SIGMA,NAME,VALUE) configures additional options specified 
% by one or more name-value pair arguments:
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
% T = JDSVD_V(A,K,SIGMA,OPTIONS) alternatively configures the additional options 
% using a structure. See the documentation below for more information.
%
% [T,U,V] = JDSVD_V(A,...) computes the singular vectors as well. If A is
% M-by-N K singular values are computed, then T is K-by-K diagonal, U and V
% are M-by-K and N-By-K orthonormal, respectively.
%
% [T,U,V,RELRES] = JDSVD_V(A,...) also returns a vector of relative residual 
% norms of the computed approximate singular triplets. 
%
% [T,U,V,RELRES,OUTIT] = JDSVD_V(A,...) also returns a vector of numbers of 
% outer iterations uesd to compute each singular triplet.
%
% [T,U,V,RELRES,OUTIT,INNIT] = JDSVD_V(A,...) also returns a vector of numbers 
% of total inner iterations uesd to compute each singular triplet.
%
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME] = JDSVD_V(A,...) also returns a vector 
% of CPU time used to compute each singular triplet.
% 
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG] = JDSVD_V(A,...) also returns a 
% convergence flag. If the method has converged, then FLAG = 0; If the maximun 
% number of outer iterations have been used before convergence, then FLAG = 1. 
%
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV] = JDSVD_V(A,...) also 
% returns a vector of relative residual norms at each outer iteration.
%
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC] = JDSVD_V(A,...)
% also returns a vector of inner iterations at each outer iteration.
%
% [T,U,V,RELRES,OUTIT,INNIT,CPUTIME,FLAG,RELRESVEV,INNERVEC,PERCENT] = JDSVD_V(A,...)
% also returns the percentage of outer iterations during which the inner iterations
% achieve the desired convergece. 
%
% [...] = JDSVD_V(AFUN,MN, ...) accepts function handle AFUN instead of the 
% matrix A. AFUN(X,'notransp') must accept a vector input X and return the 
% matrix-vector product A*X, while AFUN(U,'transp') must return A'*U. MN is 
% a 1-by-2 row vector [M N] where M and N are the numbers of rows and columns 
% of A, respectively.
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
if isfield(Options,'normAest')
    normAest = Options.normAest;
else
    normAest = svds(Afun,[M,N],1,"largest","Tolerance",1e-2,"MaxIterations",30);
end
% If target = 0, replace the left and right initial vectors u0 and v0 as
% A^Tu0 and Av0, respectively, in order to discard the information on the
% zero singular value of A and nontrivial singular triplets corresponding
% to nonzero singular values of A are computed.
if (issparse(target) && (target == 0)) || strcmp(target,'smallest') 
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
        fprintf('\nAlgorithm: the Standard JDSVD_V algorithm;\n');
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
%-------------------------------------------------%
% BEGIN: DISPLAY THE KEY INFORMATION IF REQUIRED. %
%-------------------------------------------------%
if Options.disp > 1
    fprintf('Size of A: %d×%d; \n',M,N);
    fprintf('Number of desired SVD singular triplets: %d;\n',k);
    if ischar(target)
        fprintf('Target of desired SVD singular triplets: %s;\n',target);
    else
        fprintf('Target of desired SVD singular triplets: %7.5e;\n',target);
    end 
    fprintf('Maximum number of outer iterations: %d;\n',Options.maxit);
    fprintf('Maximum dimension of searching subspaces: %d;\n',Options.maxsize);
    fprintf('Minimum dimension of searching subspaces: %d;\n',Options.minsize);
    fprintf('Switching tolerance of shift: %5.2e;\n',Options.fixtol);
    fprintf('Stopping tolerance of inner iterations: %5.2e;\n',Options.inntol);
    fprintf('1st inner preconditioning parameter: %5.2e;\n',Options.innprectol1);
    fprintf('2nd inner preconditioning parameter: %5.2e;\n',Options.innprectol2); 
    fprintf('Algorithm: the Standard JDSVD_V algorithm;\n');
end
%-----------------------------------------------%
% END: DISPLAY THE KEY INFORMATION IF REQUIRED. %
%-----------------------------------------------%
%-----------------------------------------------------------%
% BEGIN: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%-----------------------------------------------------------%
% Initialization and description of local variables.
% Converged partial SVD of A that is desired by the users.
T_conve = zeros(k,k);
U_conve = zeros(M,k);
V_conve = zeros(N,k);
% Computational results recorded to illstrate the performence of the algorithm.
RelaResidual = inf(k,1); 
OuterIterVec = zeros(k,1);
InnerIterVec = zeros(k,1);
CPUTimeseVec = zeros(k,1);
RelatResiVec = zeros(Options.maxit+k,1); 
InnerIterPer = zeros(Options.maxit,1); 
% Auxiliary items to help recording performance. 
OuterIterCounter = 0;
InnerIterCounter = 0;
k_conve = 0;   % Number of converged singular triplets.
OuterIter = 0; % Count the overall outer iterations.
InnerConv = 0; % Count the overall convergence of inner iterations.  
% Intermediate matrices to extract approximate singular triplets.
U   = zeros(M,Options.maxsize);
V   = zeros(N,Options.maxsize);
AV  = zeros(M,Options.maxsize);
ATU = zeros(N,Options.maxsize);
H   = zeros(Options.maxsize,Options.maxsize); 
% Orthonormalizing the starting left and right vectors.
if ~isempty(Options.Uconve)
    u0 = u0 - Options.Uconve*(Options.Uconve'*u0);
end
if ~isempty(Options.Vconve)
    v0 = v0 - Options.Vconve*(Options.Vconve'*v0);
end 
[U0,~] = qr(u0,'econ'); 
[V0,~] = qr(v0,'econ'); 
% Record the size of the searching subspaces.
m = min(min(size(U0,2),size(V0,2))-1,Options.maxsize);
if m > 0
    U(:,1:m) = U0(:,1:m);
    V(:,1:m) = V0(:,1:m);
    AV(:,1:m) = Afun(V(:,1:m),'notransp');
    ATU(:,1:m) = Afun(U(:,1:m),'transp');
    H(1:m,1:m) = ATU(:,1:m)'*V(:,1:m);
end
uplus = U0(:,m+1);
vplus = V0(:,m+1);
%---------------------------------------------------------%
% END: DESCRIPTION AND INITIALIZATION OF LOCAL VARIABLES. %
%---------------------------------------------------------%
%-----------------------------%
% BEGIN: MAIN ITERATION LOOP. %
%-----------------------------%
tic;
while k_conve < k && OuterIter < Options.maxit
    m = m+1;
    %-------------------------------%
    % BEGIN: Updating the matrices. %
    %-------------------------------%
    U(:,m) = uplus;
    V(:,m) = vplus;
    ATU(:,m) = Afun(uplus,'transp');
    AV(:,m) = Afun(vplus,'notransp');
    H(1:m,m) = ATU(:,1:m)'*V(:,m);
    H(m,1:m-1) = U(:,m)'*AV(:,1:m-1);
    %-----------------------------%
    % END: Updating the matrices. %
    %-----------------------------%
    %--------------------------------------------------------%
    % BEGIN: SEARCH AND EXTRACT APPROXIMATE SINGULAR VALUES. %
    %--------------------------------------------------------%
    search=1;
    while search
        % Compute the SVD of the small projected matrix.
        [Theta,Ctemp,Dtemp] = StandardExtraction(H,target,m);
        % Compute the approximate singular triplet and the corresponding residual.
        theta = Theta(1,1);  
        ufound = U(:,1:m)*Ctemp(:,1);
        vfound = V(:,1:m)*Dtemp(:,1);
        router = [AV(:,1:m)*Dtemp(:,1)-theta*ufound;ATU(:,1:m)*Ctemp(:,1)-theta*vfound];
        % Compute the relative residual norm and add it to the RelatResiVec.
        RelResidual = norm(router)/normAest;   
        RelatResiVec(OuterIter + k_conve + 1) = RelResidual;
        % Check for convergence.
        if RelResidual <= Options.tol
            % Update the number of converged singular triplets
            k_conve=k_conve+1;
            T_conve(k_conve,k_conve) = theta;
            U_conve(:,k_conve) = ufound;
            V_conve(:,k_conve) = vfound;
            % Update the computational results.
            CPUTimeseVec(k_conve) = toc;
            RelaResidual(k_conve) = RelResidual;
            OuterIterVec(k_conve) = OuterIterCounter; 
            InnerIterVec(k_conve) = InnerIterCounter;   
            if k_conve == 1 && strcmp(target,'largest')
                normAest = theta;
                Options.normAest = theta;
            end
            % If all desired eigenpairs have converged, then break the loop;
            if k_conve == k
                RelatResiVec = RelatResiVec(1:(OuterIter+k_conve));
                InnerIterPer = InnerIterPer(1:(OuterIter)); 
                break;
            else
                % Deflat and restart the algorithm efficiently.
                [U,V,ATU,AV,H] = StandardDeflation_AM(U,V,ATU,AV,Ctemp,Dtemp,Theta,m);
                m = m-1;
                % Reset the auxiliary items and the timer. 
                OuterIterCounter = 0;
                InnerIterCounter = 0;
                tic;
            end
        else
            search = 0;
        end
        %-------------------------------------------------%
        % BEGIN: DISPLAY THE KEY INFORMATION IF REQUIRED. %
        %-------------------------------------------------%
        if Options.disp > 1
            fprintf('Number of Converged singular triplets: %d;\n',k_conve);
            if k_conve >= 1
                ConvergedSVD = [(1:1:k_conve)' (diag(T_conve(1:k_conve,1:k_conve)))]';
                fprintf('Converged eigenvalues:\n');
                fprintf('%d %11.4e;\n',ConvergedSVD);  
            end
            fprintf('Relativel residual of the next approximate singular triplet: %s;\n',RelResidual);
        end
        %-----------------------------------------------%
        % END: DISPLAY THE KEY INFORMATION IF REQUIRED. %
        %-----------------------------------------------%
    end
    %------------------------------------------------------%
    % END: SEARCH AND EXTRACT APPROXIMATE SINGULAR VALUES. %
    %------------------------------------------------------%
    if k_conve == k
        break;
    end 
    % Determine the inner target for the correction equation 
    targetin = tools.GetInnerTarget(theta,target,normAest,RelResidual,...
        Options.fixtol,k_conve,T_conve); 
    if Options.innprectol1 == inf && Options.innprectol2 == inf
        iexten = []; 
    else 
    % Determine the indexes of the singular vectors to be exploited in the
    % correction equations for the first time.
    iexten = StandardExtendIndexI_AM(Theta,targetin,Options,m); 
    if isempty(iexten)
        Uexten = [];
        Vexten = [];
    else
        % Compute the singular triplets associated with the extending singular
        % vectors as well as the residual norms.
        Texten = Theta(iexten,iexten);
        Cexten = Ctemp(:,iexten);
        Dexten = Dtemp(:,iexten);
        Uexten = U(:,1:m)*Cexten;
        Vexten = V(:,1:m)*Dexten;
        Rexten = [AV(:,1:m)*Dexten-Uexten*Texten;...
            ATU(:,1:m)*Cexten-Vexten*Texten];
        Resexten = vecnorm(Rexten)/normAest;
        % Determine the indexes of the singular vectors to be finally exploited
        % in the correction equations.
        jexten = StandardExtendIndexII_AM(Resexten,Options);
        % Extract the final singular vectors to be used to precondition the
        % correction equations.
        iexten = iexten(jexten);
        Uexten = Uexten(:,jexten);
        Vexten = Vexten(:,jexten); 
    end 
    end
    % Determine the stopping cretiria epsilonin for inner iterations.
    epsilonin = StandardParamterN_AM(Theta,targetin,iexten,Options,m); 
    % Iteratively solve the correction equations approximately. 
    if Options.innprectol1 == inf && Options.innprectol2 == inf
        [s,t,iter,~,inflag] = tools.InnerSolverPreJDSVD_AM(Afun,[M N],router,...
            [U_conve(:,1:k_conve) Options.Uconve],...
            [V_conve(:,1:k_conve) Options.Vconve],...
            U(:,1),V(:,1),U(:,2:m),V(:,2:m),targetin,epsilonin); 
    else
        [s,t,iter,~,inflag] = tools.InnerSolverPreJDSVD_AM(Afun,[M N],router,...
            [U_conve(:,1:k_conve) Options.Uconve],...
            [V_conve(:,1:k_conve) Options.Vconve],...
            ufound,vfound,Uexten,Vexten,targetin,epsilonin);
    end
    % Update the computational results.
    OuterIter = OuterIter+1;
    InnerConv = InnerConv + (inflag==0);
    OuterIterCounter = OuterIterCounter+1;
    InnerIterCounter = InnerIterCounter+iter;
    InnerIterPer(OuterIter) = iter;
    % If the dimension of the searching subspaces reaches the given maximum
    % number, then efficiently restart the algorithm.
    if m == Options.maxsize 
        if Options.innprectol1 == inf && Options.innprectol2 == inf
            iexten = [];
        end
        [U,V,ATU,AV,H,m] = StandardRestart_AM(U,V,ATU,AV,Ctemp,Dtemp,Theta,...
            m,Options.minsize,iexten);
    end
    % Orthogonalize the approximate solution of the correction equation to
    % the orthonormal basis matrix of the left and right searching subspaces.
    uplus = tools.orthog([U_conve(:,1:k_conve) U(:,1:m)],s);
    vplus = tools.orthog([V_conve(:,1:k_conve) V(:,1:m)],t);
end
%---------------------------%
% END: MAIN ITERATION LOOP. %
%---------------------------%
% Restore the results about the last approximate singular triplet if the
% maximun number of outerations have been used.
if OuterIter >= Options.maxit
    % Set the flag as 1.
    FLAGofJDSVDV = 1;
    % Update the number of converged singular triplets temporarily.
    k_conve = k_conve + 1;
    T_conve(k_conve,k_conve) = theta;
    U_conve(:,k_conve) = ufound;
    V_conve(:,k_conve) = vfound;
    % Update the computational results.
    CPUTimeseVec(k_conve) = toc; 
    OuterIterVec(k_conve) = OuterIterCounter; 
    InnerIterVec(k_conve) = InnerIterCounter; 
    RelatResiVec = RelatResiVec(1:(OuterIter + k_conve));
    InnerIterPer = InnerIterPer(1:(OuterIter)); 
    % Reset the number of converged singular triplets.
    k_conve=k_conve-1;
    % Warn the users.
    warning(['Maximum number of outer iterations have been used ',...
        'before computing all the desired singular triplets']);
else
    FLAGofJDSVDV = 0;
end
% Compute the results when needed.
k_result = k_conve + FLAGofJDSVDV;
if nargout == 0 || Options.disp > 0  
    Results=zeros(k_result,6);
    Results(:,1) = (1:k_result)';
    Results(:,2) = diag(T_conve(1:k_result,1:k_result));
    Results(:,3) = RelaResidual(1:k_result,:);
    Results(:,4) = OuterIterVec(1:k_result,:);
    Results(:,5) = InnerIterVec(1:k_result,:);
    Results(:,6) = CPUTimeseVec(1:k_result,:);
end
%------------------------%
% BEGIN: OUTPUT RESULTS. %
%------------------------%
% Output option I: Display the singular values.
if nargout == 0 || Options.disp > 0 
    fprintf('Number of converged singular triplets: %d;\n',k_conve);
    fprintf('  i:    SingularValues;    RelResNorms;    out-its;    inn-its;    CPUtime; \n')
    fprintf('%3d:    %14.5e;    %11.5e;    %7d;    %7d;    %8.2e;\n',Results');
    fprintf('Total number of outer iterations = %d;\n',sum(OuterIterVec));
    fprintf('Total number of inner iterations = %d;\n',sum(InnerIterVec));
    fprintf('Total CPU time in seconds = %s;\n',sum(CPUTimeseVec));
    if OuterIter > 0
        fprintf('Percentage of converged outer iterations = %4.2f%%;\n',InnerConv/OuterIter*100);
    end
end
% Output option II: Output the singualr values.
if nargout == 1
    varargout{1} = diag(T_conve(1:k_result,1:k_result));
end
% Output option III: Output the partial SVD.
if nargout > 1
    varargout{1} = T_conve(1:k_result,1:k_result);
    varargout{2} = U_conve(:,1:k_result);
    varargout{3} = V_conve(:,1:k_result);
end
% Output option IV: Output the relative residual norm of the approximate 
% singular triplets.
if nargout >= 4
    varargout{4} = RelaResidual;
end
% Output option V: Output the numbers of outer iterations the algorithm 
% takes to computed each approximate singular triplet.
if nargout >= 5
    varargout{5} = OuterIterVec;
end
% Output option VI: Output the numbers of total inner iterations the 
% algorithm takes to computed each approximate singular triplet.
if nargout >= 6
    varargout{6} = InnerIterVec;
end
% Output option VII: Output the total CPU time the algorithm takes to 
% compute each approximate singular triplet of A.
if nargout >= 7
    varargout{7} = CPUTimeseVec;
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
    varargout{9} = RelatResiVec; 
end
% Output option X: Output the numbers of inner iterations during each 
% outer iteration.
if nargout >= 10
    varargout{10} = InnerIterPer;
end
% Output Options XI: Oupput the percentage of total outer iterations during
% which the inner iterations have converged to the desired inner accuracy.
if nargout >= 11
    varargout{11} = InnerConv/OuterIter;
end
%----------------------%
% END: OUTPUT RESULTS. %
%----------------------%
end
 
function [Theta,Ctemp,Dtemp] = StandardExtraction(H,target,m)
% Compute the SVD of the small projected matrix and resort the approximate
%singular values in the desired oder.
% Compute the small scaled SVD.
[Ctemp,Theta,Dtemp] = svd(H(1:m,1:m));
% Resort the approximate singular values and the Ritz vectors.
Isort = tools.SortSigmas(diag(Theta),target);
Theta = Theta(Isort,Isort);
Ctemp = Ctemp(:,Isort);
Dtemp = Dtemp(:,Isort);
end

function [U,V,ATU,AV,H]=StandardDeflation_AM(U,V,ATU,AV,Ctemp,Dtemp,Theta,m)
% Perform one step of restart after the deflation.
U(:,1:m-1) = U(:,1:m)*Ctemp(:,2:m);
V(:,1:m-1) = V(:,1:m)*Dtemp(:,2:m);
ATU(:,1:m-1) = ATU(:,1:m)*Ctemp(:,2:m);
AV(:,1:m-1) = AV(:,1:m)*Dtemp(:,2:m);
H(1:m-1,1:m-1) = Theta(2:m,2:m);
end

function iexten = StandardExtendIndexI_AM(Theta,targetin,Options,m)
% Determine the indexes of the singular vectors to be extended to the
% correction equations for the first time.
if m > 1
    theta = diag(Theta(2:m,2:m));
    iexten = 1+find(abs(theta-targetin)./max(theta,1) <= Options.innprectol1);
else
    iexten = [];
end
iexten=iexten';
end

function iexten = StandardExtendIndexII_AM(Resexten,Options)
% Determine the indexes of the singular vectors to be finally exploited
% in the correction equations.
iexten = find(Resexten <= Options.innprectol2);
end

function epsilonin = StandardParamterN_AM(Theta,targetin,iexten,Options,m)
% Determine the stopping tolerance epsilonin for the inner iterations to
% solve the correction equation.
% Compute the accuracy requirement for solution.
guard = 0.1;
ConstantParameter = 2*sqrt(2);
if m>1
    Theta = diag(Theta);
    theta = Theta(1);   
    Thetaus = Theta(setdiff(2:1:m,iexten));
    if ~isempty(Thetaus) 
        relativesep = 1-(theta-targetin)./(Thetaus-targetin);  
        ConstantParameter = ConstantParameter/min(abs(relativesep));
    end 
end
% Compute inner stopping criteria.
epsilonin = min([Options.inntol*ConstantParameter guard]);
end

function [U,V,ATU,AV,H,m] = StandardRestart_AM(U,V,ATU,AV,Ctemp,Dtemp,Theta,maxsize,minsize,iexten)
% Perform the thick restart techenique for the standard extraction approach.
if length(iexten)+1<=minsize
    inew = 1:1:minsize;
else
    inew = [1 iexten];
end
m = length(inew);
% Update the indermediate matrices
U(:,1:m) = U(:,1:maxsize)*Ctemp(:,inew);
V(:,1:m) = V(:,1:maxsize)*Dtemp(:,inew);
ATU(:,1:m) = ATU(:,1:maxsize)*Ctemp(:,inew);
AV(:,1:m) = AV(:,1:maxsize)*Dtemp(:,inew);
H(1:m,1:m) = Theta(inew,inew);
end