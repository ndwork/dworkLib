
classdef parforProgress
  % Function measures progress of an iteration implemented with parfor
  %
  % Usage:
  %
  % N = 100;
  % p = parforProgress( N, tmpFile );  % Note: if tmpFile is not supplied
  %                                    % then parforProgress.txt is used.
  % parfor n=1:N
  %   p.progress( n );
  %   pause(rand*10); % Replace with real code
  % end
  % p.clean;
  %
  % Notes:
  % progress member function accepts optional downsampling input.
  %   Ex: p.progress( n, 10 );  % displays only when mod(n,10) == 0
  %
  % Written by Nicholas Dwork - Copyright 2017
  % Based on parfor_progress written by Jeremy Scheff
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  properties
    nTotal
    tmpFile
  end

  methods
    % Constructor
    function obj = parforProgress( nTotal, tmpFile )
      % nTotal is the total number of iterations for completion
      if nargin < 1, error('Must supply total number of iterations.'); end;

      pid = feature('getpid');
      if nargin < 2
        tmpFile = ['parforProgress_', num2str(pid), '.txt'];
      end

      obj.nTotal = nTotal;
      obj.tmpFile = tmpFile;
      fid = fopen( obj.tmpFile, 'w' );
      fclose(fid);
    end

    % Destructor
    function delete( obj )
      delete( obj.tmpFile );
    end

    % Member functions
    function clean( obj )
      delete( obj.tmpFile );
    end

    function progress( obj, n, varargin )
      p = inputParser;
      p.addOptional( 'downsample', 1, @isnumeric );
      p.parse( varargin{:} );
      downsample = p.Results.downsample;

      fid = fopen( obj.tmpFile, 'a' );
      if fid<0, error( 'Unable to open parforProgress.txt temporary file' ); end;
      nowTime = datestr( now, 'yyyy-mm-dd HH:MM:SS.FFF');
      fprintf( fid, [ '%d: ', nowTime, '\n' ], n );
        % n is the index of the current iteration
      fclose( fid );
      nLines = findNumLinesInFile( obj.tmpFile );
      if mod( n, downsample )==0
        disp(['Working on ', num2str(n), ' of ', num2str(obj.nTotal), ...
          ': ', num2str( nLines / obj.nTotal * 100 ), '%' ]);
        drawnow( 'update' );
      end
    end

  end

end
