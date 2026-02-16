function varargout = MatrixFunction(C)
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