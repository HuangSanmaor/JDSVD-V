function [s,t,iter,relres,flag] = InnerSolverPreJDSVD_AM(Afun,MN,router,...
    Uconve,Vconve,ufound,vfound,Uexten,Vexten,target,inepsilon)
% Solve the correction equations involved in the JD type GSVD algorithms. 
 
M = MN(1);  Up = [Uconve ufound Uexten]; 
N = MN(2);  Vp = [Vconve vfound Vexten];  
% Calculate the projected right hand side to make the equation consistant.
router_projection(1:M,:) = orthog(Uconve,router(1:M,:)); 
router_projection((M+1):(M+N),:) = orthog(Vconve,router((M+1):(M+N),:));
% Creat the function handle for applying matrix-vector multiplications with
% the coefficient matrix of the correction equation.
mvpwithB = MatricesFunction(Afun,M,N,Up,Vp,target); 
% Turn off the warning.
warning off; 
% Solve the correction equation using the MINRES algorithm.
[st,flag,relres,iter] = minres(mvpwithB,router_projection,inepsilon,M+N);  
% Turn on the warning.
warning on; 
s = st(1:M,:); 
t = st((1+M):(M+N),:); 
end
 
function varargout = MatricesFunction(Afun,M,N,Up,Vp,target)
% Creat the function handle involved in innersolver.
varargout{1}=@matricesfunction;
function y = matricesfunction(x)  
   yupper = Afun(x((M+1):(M+N),:),'notransp')-target*x(1:M,:);
   ylower = Afun(x(1:M,:),'transp')-target*x((M+1):(M+N),:);
   y(1:M,:) = orthog(Up,yupper);
   y((M+1):(M+N),:) = orthog(Vp,ylower);
end 
% function y = matricesfunction(x) 
%     xupper = orthog(Up,x(1:M,:));
%     xlower = orthog(Vp,x((M+1):(M+N),:));
%    yupper = Afun(xlower,'notransp')-target*xupper;
%    ylower = Afun(xupper,'transp')-target*xlower;
%    y(1:M,:) = orthog(Up,yupper);
%    y((M+1):(M+N),:) = orthog(Vp,ylower);
% end 
end
 
function y = orthog(X,y)
% Orthogonalize a vector to a orthonormal matrix. 
if isempty(X) 
    return;
end
alpha = 0.5;
normy_ori = norm(y);   
normy_old = normy_ori;
y = y - X*(X'*y);   
normy = norm(y);
while normy/normy_old < alpha
    if normy > normy_ori*0.1*eps
        normy_old = normy;
    else
        normy_ori = 0.1*normy_old*eps;
        normy_old = normy_ori;
        vecnormXi = vecnorm(X',1,1); 
        e = zeros(size(y));
        e(find(vecnormXi == min(vecnormXi),1)) = 1;
        y = normy_old*e;
    end
    y = y - X*(X'*y);
    normy = norm(y);
end  
end