
classdef parforProgress
  % Function measures progress of an iteration implemented with parfor
  %
  % Usage:
  %
  % N = 100;
  % p = parforProgress( N );
  % parfor n=1:N
  %   p.progress( n );
  %   pause(rand*10); % Replace with real code
  % end
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
    function obj = parforProgress( nTotal )
      % nTotal is the total number of iterations for completion
      if nargin < 1, error('Must supply total number of iterations.'); end;
      obj.nTotal = nTotal;
      obj.tmpFile = 'parforProgress.txt';
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
      obj.nTotal = 0;
    end
    
    function progress( obj, n, varargin )
      defaultDownsample = 1;
      p = inputParser;
      p.addOptional( 'downsample', defaultDownsample, @isnumeric );
      p.parse( varargin{:} );
      downsample = p.Results.downsample;
      
      % n is the index of the current iteration
      fid = fopen( obj.tmpFile, 'a' );
      if fid<0, error( 'Unable to open parforProgress.txt temporary file' ); end;
      fprintf( fid, '%d\n', n ); % Save n at the top of progress.txt
      fclose( fid );
      nLines = findNumLinesInFile( obj.tmpFile );
      if mod( n, downsample )==0
        disp(['Working on ', num2str(n), ' of ', num2str(obj.nTotal), ...
          ': ', num2str( nLines / obj.nTotal * 100 ), '%' ]);
        drawnow('update');
      end
    end

  end

end
