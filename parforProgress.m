
classdef parforProgress
  % p = parforProgress( [ N, tmpFile, 'msgHdr', msgHdr, 'msgFtr', msgFtr ] );
  % p.progress( n, D, 'msgHdr', msgHdr, 'msgFtr', 'msgFtr' );
  % p.clean;
  % 
  % Function measures progress of an iteration implemented with parfor
  %
  % Usage:
  %
  % N = 100;
  % p = parforProgress( N, tmpFile );  % Note: if tmpFile is not supplied
  %                                    % then parforProgress.txt is used.
  % parfor n=1:N
  %   p.progress( n, D );
  %   pause(rand*10); % Replace with real code
  % end
  % p.clean;
  %
  % Optional Inputs:
  % D - downsampling input.
  %   Ex: p.progress( n, 10 );  % displays only when mod(n,10) == 0
  % tmpFile - the file to use to store progress information
  % msgHdr - a string to append to the beginning of progress statements.
  %   By default: this string is the calling function's name.
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
    msgHdr
    msgFtr
  end

  methods
    % Constructor
    function obj = parforProgress( nTotal, varargin )
      % nTotal is the total number of iterations for completion
      if nargin < 1, error('Must supply total number of iterations.'); end

      pid = feature('getpid');
      defaultTmpFile = ['parforProgress_', num2str(pid), '.txt'];

      st = dbstack;
      defaultMsgHdr = [ st(2).name, ':  ' ]; % The function caller's name (parent)
      defaultMsgFtr = [];

      p = inputParser;
      p.addOptional( 'tmpFile', defaultTmpFile, ...
        @(x) validateattributes(x,{'char'}) );
      p.addParameter( 'msgHdr', defaultMsgHdr, @(x) true );
      p.addParameter( 'msgFtr', defaultMsgFtr, @(x) true );
      p.parse( varargin{:} );
      obj.tmpFile = p.Results.tmpFile;
      obj.msgHdr = p.Results.msgHdr;
      obj.msgFtr = p.Results.msgFtr;

      obj.nTotal = nTotal;
      fid = fopen( obj.tmpFile, 'w' );
      fclose(fid);
    end

    % Destructor
    function delete( obj )
      delete( obj.tmpFile );
    end

    % Member functions
    function clean( obj, varargin )
      p = inputParser;
      p.addOptional( 'type', [], @(x) true );
      p.parse( varargin{:} );
      type = p.Results.type;

      if strcmp( type, 'all' )
        parforProgressFiles = dir( './parforProgress*.txt' );
        for fileIndx = 1 : numel( parforProgressFiles )
          delete( parforProgressFiles(fileIndx).name );
        end
      end

      if exist( obj.tmpFile, 'file' )
        delete( obj.tmpFile );
      end
    end

    function progress( obj, n, varargin )
      p = inputParser;
      p.addOptional( 'downsample', 1, @isnumeric );
      p.addParameter( 'msgHdr', [], @(x) true );
      p.addParameter( 'msgFtr', [], @(x) true );
      p.parse( varargin{:} );
      downsample = p.Results.downsample;
      moreMsgHdr = p.Results.msgHdr;
      moreMsgFtr = p.Results.msgFtr;

      fid = fopen( obj.tmpFile, 'a' );
      if fid<0, error( 'Unable to open parforProgress.txt temporary file' ); end
      nowTime = datestr( now, 'yyyy-mm-dd HH:MM:SS.FFF');
      fprintf( fid, [ '%d: ', nowTime, '\n' ], n );
        % n is the index of the current iteration
      fclose( fid );
      nLines = findNumLinesInFile( obj.tmpFile );
      if mod( nLines+1, downsample ) == 0
        disp([ obj.msgHdr, moreMsgHdr, 'Working on ', indx2str( nLines, obj.nTotal ), ' of ', ...
          num2str(obj.nTotal), ': ', num2str( nLines / obj.nTotal * 100 ), '%', ...
          '   index: ', indx2str( n, obj.nTotal ), '  ', obj.msgFtr, moreMsgFtr ]);
        drawnow( 'update' );
      end
    end

  end

end
