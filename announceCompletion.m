
function announceCompletion( varargin )
  % announceCompletion( [ 'earlyTime', earlyTime, 'lateTime', lateTime ] )
  %
  % This function makes an audible sound when the function is called.
  % Intended to be run at the end of the program to let you know it's done.
  %
  % Inputs:
  % earlyTime - specified as hh:mm 24 hour format.  If specified, won't make a sound
  %   before this time in the morning.
  % lateTime - specified as hh:mm 24 hour format.  If sepcified, won't make a sound 
  %   after this time in the evening.
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'earlyTime', [], @(x) true );
  p.addParameter( 'lateTime', [], @(x) true );
  p.parse( varargin{:} );
  earlyTime = p.Results.earlyTime;
  lateTime = p.Results.lateTime;

  t = clock();

  if numel( earlyTime ) > 0
    earlyTimeParts = strsplit( earlyTime, ':' );
    earlyHour = str2double( earlyTimeParts(1) );
    earlyMin = str2double( earlyTimeParts(2) );
    if( t(4) <= earlyHour )
      if( t(4) < earlyHour ), return; end;
      if( t(4) == earlyHour && t(5) <= earlyMin ), return; end;
    end
  end

  if numel( lateTime ) > 0
    lateTimeParts = strsplit( lateTime, ':' );
    lateHour = str2double( lateTimeParts(1) );
    lateMin = str2double( lateTimeParts(2) );
    if( t(4) >= lateHour )
      if( t(4) == lateHour && t(5) >= lateMin ), return; end;
      if( t(4) > lateHour ), return; end;
    end
  end

  if exist( 'programComplete.wav', 'file' )
    [y,Fs] = audioread( 'programComplete.wav' );
  else
    load( 'gong' );
  end
  soundsc( y, Fs );

end

