function targetin = GetInnerTarget(theta,target,normAest,RelResidual,...
    fixtol,k_conve,T_conve)
% Determine the parameter targetin used in the followed correction equation
if RelResidual <= fixtol
    targetin = theta;
    return;
end
if isscalar(target)
    targetin = target;
    return;
end 
if strcmp(target,'largest') 
    if k_conve <= 0
        targetin = 1.1*normAest;
    else
        targetin = T_conve(k_conve,k_conve);
    end
    return;
end
if strcmp(target,'smallest') || (isscalar(target) && target ==0)  
    if k_conve <= 0
        targetin = 0;
    else
        targetin = T_conve(k_conve,k_conve);
    end 
    return;
end
end