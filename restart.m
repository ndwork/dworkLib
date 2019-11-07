
function restart( varargin )
  % restart( [ 'all' ] )
  %
  % closes all windows, clears variables, clears the command window, ends debug
  % mode, and deletes any parforProgess files that were created by the current
  % instance of Matlab
  %
  % Optional Inputs:
  % if 'all' is supplied, restart removes any parforProgress.txt files that exist
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  type = [];
  if nargin > 0, type = varargin{1}; end
  
  if strcmp( type, 'all' )
    clear all;  %#ok<CLALL>
    evalin( 'base', 'clear all' );  % delete all variables in global scope

    parforFiles = dir( 'parforProgress*.txt' );
    for fileIndx = 1 : numel( parforFiles )
      delete( parforFiles(fileIndx).name );
    end

  else

    pid = feature('getpid');
    parforProgressFile = ['parforProgress_', num2str(pid), '.txt'];
    if exist( parforProgressFile, 'file' )
      delete( parforProgressFile );
    end

  end

  close all;  clc;  clear;  clear global;

  % I would like to be able to call 'dbquit' here, but Matlab
  % prevents it.  So instead, I must simulate keyboard commands,
  % which is what is done here:
  for i=1:5, pressShiftF5;  pause( 0.05 ); end
end


function pressShiftF5
  import java.awt.Robot;
  import java.awt.event.KeyEvent;

  rob = Robot;  %Create a Robot-object to do the key-pressing
  rob.keyPress( KeyEvent.VK_SHIFT );
  rob.keyPress( KeyEvent.VK_F5 );
  rob.keyRelease( KeyEvent.VK_F5 );
  rob.keyRelease( KeyEvent.VK_SHIFT );
end


