
function restart( varargin )
  % restart
  %
  % closes all windows, clears variables, clears the command window, and
  % removes the parforProgress.txt file (if it exists)
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  type = [];
  if nargin > 0, type = varargin{1}; end;

  if strcmp( type, 'all' ), clear all; end;                                                %#ok<CLALL>

  parforFiles = dir( 'parforProgress*.txt' );
  for fileIndx = 1 : numel( parforFiles )
    delete( parforFiles(fileIndx).name );
  end

  close all;  clear;  clc;
end
