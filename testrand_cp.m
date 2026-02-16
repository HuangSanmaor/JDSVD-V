M = 3000; 
sigA = [0.8:0.0001:0.9979,1.0001:0.0001:1.002];
%sigA = [0.8001:0.0001:0.802,1.0021:0.0001:1.2];

N = size(sigA,2); 
dense = 0.02; 
A = sprand(M,N,dense,sigA);
afun = Matrix2Function(A); 
fprintf('done generating the problem;\n');

k = 20;
target = 'smallest';
options = struct;
outtol = 1e-12;          options.tol = outtol;
maxit = min(M,N);        options.maxit = maxit;
maxsize = 30;            options.maxsize = maxsize;
minsize = 3;             options.minsize = minsize;
fixtol = 1e-4;           options.fixtol = fixtol; 
inntol = 1e-3;           options.inntol = inntol;
if (isscalar(target) && target == 0) || strcmp(target,'smallest')
    u0 = A*ones(N,1);    options.u0 = u0;
    v0 = A'*ones(M,1);   options.v0 = v0;
else
    u0 = ones(M,1);      options.u0 = u0;
    v0 = ones(N,1);      options.v0 = v0;
end
display = 1;             options.disp = display;
innprectol1 = 0.05;      options.innprectol1 = innprectol1;
innprectol2 = 0.01;      options.innprectol2 = innprectol2;

innprectol1 = 0.05;
innprectol2 = 0.01;
% u0=randn(M,3);
% v0=randn(N,5);
 

% Tests the JDSVD_V_CP algorithm with matrices;
[Tp,Up,Vp,RELRESp,OUTITp,INNITp,CPUTIMEp,FLAGp,RELRESVEVp,INNERVECp,PERCp] = ...
    jdsvd_v_cp(afun,[M N],k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,...
    'InnerPreconditionTolerance1',innprectol1,...
    'InnerPreconditionTolerance2',innprectol2,'Display',display); 
[To,Uo,Vo,RELRESo,OUTITo,INNITo,CPUTIMEo,FLAGo,RELRESVEVo,INNERVECo,PERCo] =...
    jdsvd_v_cp(A,k,target,'Tolerance',outtol,'MaxIterations',maxit,...
    'MaxSubspaceDimension',maxsize,'MinSubspaceDimension',minsize,...
    'SwitchingTolerance',fixtol,'InnerTolerance',inntol,...
    'LeftStartVector',u0,'RightStartVector',v0,'Display',display,...
    'InnerPreconditionTolerance1',0,'InnerPreconditionTolerance2',0);
% % Tests the JDSVD_V_CP algorithm with function handles and/or options;
jdsvd_v_cp(A,k,target,options);
jdsvd_v_cp(afun,[M N],k,target,options);

% MatrixFunction turns a matrix to a function handle, for experimental purpose only
function varargout = Matrix2Function(C)
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