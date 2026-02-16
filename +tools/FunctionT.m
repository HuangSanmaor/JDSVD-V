function varargout = FunctionT(C)
% Creat the function handle for the transepose C^T of C with a given 
% function handle for C.
varargout{1} = @matrixfun;
    function y = matrixfun(x,transpornot)
        if strcmp(transpornot,'notransp')
            y = C(x,'transp');
        elseif strcmp(transpornot,'transp')
            y = C(x,'notransp');
        else
            error('wrong input for transpornot');
        end
    end
end