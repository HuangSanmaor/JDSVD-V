function [A,M,N,k,target,u0,v0,Options] = chekInputs(A,VarIn,randStr)
% Get the inputs and do check errors.

% Get A and M, N.
[M,N,VarIn] = getSizeA(A,VarIn);
% Get the number k of desired singular triplets.
if numel(VarIn)<1
    k = 6;
else
    k = VarIn{1};
    if ~isPosInt(k) || ~isscalar(k)
        error('Wrong input of k');
    end
    k = double(full(k));
end
k = min([M,N,k]);
% Get the target of the desired singular values.
target = getTarget(VarIn);
% Get the starting vector x0 and optional parameters.
[u0,v0,Options] = getOptions(VarIn,M,N,k,randStr);
end

function [M,N,VarIn] = getSizeA(A,VarIn)
% Error check and get the size of A in checkInputs.
if isa(A,'function_handle')
    if numel(VarIn) < 1
        error('For a function handle A, its sizes should be given');
    end
    % MN gives the sizes of A.
    MN=VarIn{1};
    % Error check M and N.
    if  ~isPosInt(MN) || ~isrow(MN) || length(MN) ~= 2 || ~all(isfinite(MN))
        error('Wrong inputs for the size numbers M, N');
    else
        M = double(full(MN(1)));
        N = double(full(MN(2)));
    end
    % Remove MN from VarIn. The remaining entries are k, sigma, and Options
    % which matches VarIn when A is given as a matrix.
    VarIn=VarIn(2:numel(VarIn));
elseif ismatrix(A)
    % Save size of A in [M,N].
    [M,N] = size(A);
else
    error('A should be either a function handle or a matrix');
end
end

function target = getTarget(VarIn)
% Error check and get the target in checkInputs.
if length(VarIn) < 2
    target = 'largest';
else
    target = VarIn{2};
    % Error Check sigma.
    if (ischar(target) && isrow(target)) || (isstring(target) && isscalar(target))
        ValidTarget = {'largest','smallest'};
        match = startsWith(ValidTarget, target, 'IgnoreCase', true);
        j = find(match,1);
        if isempty(j) || (strlength(target) == 0)
            error('The target should be one of the strings: "largest" and "smallest"');
        else
            % Reset sigma to the correct valid sigma for cheaper checking.
            target = ValidTarget{j};
        end
    elseif isfloat(target)
        if ~isreal(target) || ~isscalar(target) || ~isfinite(target)
            error('The target should be a real finite scalar');
        end
    end
end
end

function [u0,v0,Options] = getOptions(VarIn,M,N,k,randStr)
% Get x0 and options, and set defaults if they are not provided.
Options = struct;
nVarIn = numel(VarIn);
NameValueFlag = false;
% Get Options, if provided.
if nVarIn >=3
    if isstruct(VarIn{3})
        Options = VarIn{3};
        % VarIn should be {k, sigma, Options}.
        if nVarIn > 3
            error('Only one of the structure "opts" and name-value pairs pattern can be used to specify the parameters' );
        end
    else
        % Convert the Name-Value pairs to a struct for ease of error checking.
        NameValueFlag = true;
        for j=3:2:nVarIn
            name = VarIn{j};
            if (~(ischar(name) && isrow(name)) && ~(isstring(name) && isscalar(name))) ...
                    || (isstring(name) && strlength(name) == 0)
                error('Wtong input for the names');
            end
            nvNames=["Tolerance","MaxIterations","MaxSubspaceDimension",...
                "MinSubspaceDimension","SwitchingTolerance","InnerTolerance",...
                "InnerPreconditionTolerance1","InnerPreconditionTolerance2",...
                "LeftConvergedVectors","RightConvergedVectors","LeftStartVector",...
                "RightStartVector","Display","NormAEstimation"];
            ind = matches(nvNames, name, 'IgnoreCase', true);
            if nnz(ind) ~= 1
                error('There is a repeated Name-Value parameter pair;\n')
            end
            if j+1 > nVarIn
                error('The Name-Value parameters should appear in pair;\n');
            end
            structNames = {'tol','maxit','maxsize','minsize','fixtol','inntol',...
                'innprectol1','innprectol2','Uconve','Vconve','u0','v0',...
                'disp','normAest'};
            Options.(structNames{ind}) = VarIn{j+1};
        end
    end
end
% Error check Options.tol if provided or set the default value for it if not.
if isfield(Options,'tol')
    if ~isnumeric(Options.tol) || ~isscalar(Options.tol) ...
            || ~isreal(Options.tol) || ~(Options.tol >= 0)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.tol = 1e-8;
end
% Error check Options.maxit if provided or set the default value for it if not.
if isfield(Options,'maxit')
    if ~isPosInt(Options.maxit) || ~isscalar(Options.maxit)...
            || Options.maxit == 0
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.maxit = min([M,N]);
end
% Error check Options.maxsize if provided or set the default value for it if not.
if isfield(Options,'maxsize')
    if ~isPosInt(Options.maxsize) || ~isscalar(Options.maxsize)...
            || Options.maxsize == 0
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
    if Options.maxsize < k+2
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.maxsize = max([30 k+3]);
end
% Error check Options.minsize if provided or set the default value for it if not.
if isfield(Options,'minsize')
    if ~isPosInt(Options.minsize) || ~isscalar(Options.minsize)...
            || Options.minsize == 0
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.minsize = 3;
end
% Error check Options.fixtol if provided or set the default value for it if not.
if isfield(Options,'fixtol')
    if ~isnumeric(Options.fixtol) || ~isscalar(Options.fixtol) ...
            || ~isreal(Options.fixtol) || ~(Options.fixtol >= 0)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.fixtol = 1e-4;
end
% Error check Options.inntol if provided or set the default value for it if not.
if isfield(Options,'inntol')
    if ~isnumeric(Options.inntol) || ~isscalar(Options.inntol) ...
            || ~isreal(Options.inntol) || ~(Options.inntol >= 0)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.inntol = 1e-3;
end
% Error check Options.innprec if provided or set the default value for it if not.
if isfield(Options,'innprectol1')
    if ~isnumeric(Options.innprectol1) || ~isscalar(Options.innprectol1) || ~isreal(Options.innprectol1)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.innprectol1 = 0.05;
end
% Error check Options.inntol if provided or set the default value for it if not.
if isfield(Options,'innprectol2')
    if ~isnumeric(Options.innprectol2) || ~isscalar(Options.innprectol2) ...
            || ~isreal(Options.innprectol2) || ~(Options.innprectol2 >= 0)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
elseif Options.innprectol1 > 0
    Options.innprectol2 = 0.01;
end
% Error check Options.Uconve if provided or set the default value for it if not.
if isfield(Options,'Uconve')
    if ~isempty(Options.Uconve) && (~ismatrix(Options.Uconve) || ...
            size(Options.Uconve,1) ~= M || ~isfloat(Options.Uconve)) 
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.Uconve = []; 
end
% Error check Options.Vconve if provided or set the default value for it if not.
if isfield(Options,'Vconve') 
    if ~isempty(Options.Vconve) && (~ismatrix(Options.Vconve) || ...
            size(Options.Vconve,1) ~= N || ~isfloat(Options.Vconve))
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end

    end
else
    Options.Vconve = []; 
end
% Error check Options.u0 if provided or set the default value for it if not.
if isfield(Options,'u0')
    u0 = Options.u0;
    if ~ismatrix(u0) || size(u0,1) ~= M || ~isfloat(u0)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    u0 = randn(randStr,M,1); 
end 
% Error check Options.v0 if provided or set the default value for it if not.
if isfield(Options,'v0')
    v0 = Options.v0;
    if ~ismatrix(v0) || size(v0,1) ~= N || ~isfloat(v0)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    v0 = randn(randStr,N,1);
end 
% Error check Options.disp if provided or set the default value for it if not.
if isfield(Options,'disp')
    if ~isnumeric(Options.disp) || ~isscalar(Options.disp) || ~isreal(Options.disp)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
else
    Options.disp = 0;
end 
% Error check Options.normAest if provided
if isfield(Options,'normAest')
    if ~isnumeric(Options.normAest) || ~isscalar(Options.normAest) || ~isreal(Options.normAest)
        if NameValueFlag
            error('Wrong input in the name-value pairs');
        else
            error('Wrong input in the structure "opts"');
        end
    end
end
end

function  [tf] = isPosInt(X)
% Check if X is a non-negative integer vector.
tf = isnumeric(X) && isreal(X) && all(X(:) >= 0) && all(fix(X(:)) == X(:));
end