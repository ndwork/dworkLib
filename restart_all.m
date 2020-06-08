
% restart_all
%
% closes all windows, clears variables, clears the command window, ends debug
% mode, and deletes any parforProgess files that were created by the current
% instance of Matlab
%
% Written by Nicholas Dwork - Copyright 2020
%
% https://github.com/ndwork/dworkLib.git
%
% This software is offered under the GNU General Public License 3.0.  It
% is offered without any warranty expressed or implied, including the
% implied warranties of merchantability or fitness for a particular
% purpose.


clear all;  %#ok<CLALL>
evalin( 'base', 'clear all' );  % delete all variables in global scope

parforFiles = dir( 'parforProgress*.txt' );
for fileIndx = 1 : numel( parforFiles )
  delete( parforFiles(fileIndx).name );
end

restart

