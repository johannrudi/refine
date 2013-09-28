%REFINE_PRINTL Print a line of refine function on screen.
%
% Syntax:
% REFINE_PRINTL(string)
%
% Input:
% string  string that will be printed
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function refine_printl (string)

PREFIX = '[refine]';
fprintf([PREFIX ' ' string '\n']);

% end function
end
