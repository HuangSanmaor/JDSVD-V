function [t,iter,relres,flag] = InnerSolverPreJDSVD_CP(Afun,N,router,...
    Vconve,vfound,Vexten,target,inepsilon)
% Solve the correction equations involved in the JD type GSVD algorithms.
%Afun,MN,router,U_conve,V_conve,ufound,vfound,targetin,epsilonin
 
Vp = [Vconve vfound Vexten]; 
% Calculate the projected right hand side to make the equation consistant.
router_projection = orthog([Vconve Vexten],router); 
% Creat the function handle of applying matrix-vector multiplication with
% the coefficient matrix.
mvpwithB = MatricesFunction(Afun,Vp,target); 
% Turn off the warning.
warning off; 
% Solve the correction equation using the MINRES algorithm.
[t,flag,relres,iter] = minres(mvpwithB,router_projection,inepsilon,N); 
% Turn on the warning.
warning on;   
end
 
function varargout = MatricesFunction(Afun,Vp,target)
% Creat the function handle involved in innersolver.
varargout{1}=@matricesfunction;
function y = matricesfunction(x)
   %x = orthog(Vp,x);
   y = Afun(Afun(x,'notransp'),'transp')-target^2*x;
   y = orthog(Vp,y);
end
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