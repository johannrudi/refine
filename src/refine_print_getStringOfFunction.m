%REFINE_PRINT_GETSTRINGOFFUNCTION Generate human readable string for function.
% Generates a formatted string that describes a function with derivatives.
%
% Syntax:
% string = REFINE_PRINT_GETSTRINGOFFUNCTION (fnc_name, derivatives)
%
% Input:
% fnc_name       name of the function
% derivatives[]  function derivatives
%
% Output:
% string         formatted string that describes function with derivatives
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function string = refine_print_getStringOfFunction (fnc_name, derivatives)

% set dimension
dim = size(derivatives, 2);

% write derivatives to string
str_deriv = num2str(derivatives, '%d,');

% set name of function
if isequal(derivatives, zeros(1, dim)) % if no derivatives
    string = fnc_name;
else
    if dim == 1 % if 1-dimensional
        string = ['(D^' str_deriv(1:end-1) ' ' fnc_name ')'];
    else
        string = ['[D^(' str_deriv(1:end-1) ') ' fnc_name ']'];
    end
end

% end function
end
