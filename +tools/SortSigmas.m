function Isort = SortSigmas(theta,target)
% Reorder the approximate singular values.
if strcmp(target,'smallest')
    [~,Isort] = sort(theta);
elseif strcmp(target,'largest')
    [~,Isort] = sort(theta,'descend');
elseif isscalar(target)  
    [~,Isort] = sort(abs(theta-target));
else
    error('wrong input for target');
end
end