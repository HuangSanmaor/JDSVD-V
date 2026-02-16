function y = orthog(X,y)
% Orthogonalize a vector to a orthonormal matrix. 

if isempty(X)
    normy = norm(y);
    y = (1/normy)*y;
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
y = (1/normy)*y;  
end