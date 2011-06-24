%FFTNPAD Pad an fft ordered array.
%   B = FFTNPAD(A) pads an array A such that size(B) == 2*size(A).
%
%   B = FFTNPAD(A,PADSIZE) pads an array A with PADSIZE(k) zeros along the
%   k-th dimension of A.  
%
%   B = FFTNPAD(A,PADSIZE,PADVAL) pads array A with PADVAL (a scalar)
%   instead of with zeros. 
%
%   If the dimensions specified in PADSIZE are smaller than the initial
%   dimensions of the signal it us un-padded.  
%
%See also: FFTSPACE

% $Author: graceej $ $Date: 2009/05/06 20:09:45 $
% $Revision: 1.3 $

% (C) Copyright 2009 Edward J. Grace
% 
% This file is part of Gaffe.
% 
% Gaffe is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Gaffe is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Gaffe.  If not, see <http://www.gnu.org/licenses/>.

function B=fftnpad(A,varargin)
%% Manage arguments.
if nargin < 2
    PADSIZE=size(A)*2;
else
    PADSIZE=varargin{1};
end
if nargin < 3
    VALUE=0;
else
    VALUE=varargin{2};
end

SZ1=size(A);
ND=numel(SZ1);
%% Determine if padding or truncation is required for each dimension.
pad = zeros(ND,1);
for d=1:ND
    if PADSIZE(d) > SZ1(d)
        pad(d) = 1;
    else
        pad(d) = 0;
    end
end


% Setup output array.
B = ones(PADSIZE)*VALUE;

% Swap source and target depending if we are padding or unpadding.
for d=1:ND
    if ~pad(d)
        tmp=SZ1(d); SZ1(d) = PADSIZE(d); PADSIZE(d) = tmp;
    end
end
%% Determine source and target ranges for each dimension.
for d=1:ND
    InHalf=floor(SZ1(d)/2);
    Offset = mod(SZ1(d),2);
    
    % Determine source and target ranges for M
    src{d} = [1:(InHalf+Offset) (InHalf+Offset+1):SZ1(d)];
    tgt{d} = [1:(InHalf+Offset) (PADSIZE(d)-InHalf+1):PADSIZE(d)];
    
    % If we are unpadding, swap the source and target ranges.
    if ~pad(d)
        tmp=src{d}; src{d}=tgt{d}; tgt{d}=tmp;
    end
   
end
%% Assign source to target using cell {:} notation of indices (see fftshift)
B(tgt{:}) = A(src{:});
%% Log of changes
%
% $Log: fftnpad.m,v $
% Revision 1.3  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/05/06 17:57:10  graceej
% * Added CVS Revision information.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:25  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.1  2008/09/16 11:28:29  graceej
% * Name change from fftpadn.
%
% Revision 1.1.2.1  2008/09/05 20:38:57  graceej
% * Initial checkin of fftnpad.  It appears to work as advertised.
%
%