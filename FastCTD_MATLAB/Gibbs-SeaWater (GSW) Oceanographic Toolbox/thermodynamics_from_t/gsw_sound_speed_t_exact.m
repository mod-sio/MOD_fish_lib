function sound_speed_t_exact = gsw_sound_speed_t_exact(SA,t,p)

% gsw_sound_speed_t_exact                                       sound speed
%==========================================================================
%
% USAGE:  
%  sound_speed_t_exact = gsw_sound_speed_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the speed of sound in seawater. 
%  
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  sound_speed_t_exact  =  speed of sound in seawater               [ m/s ]
%
% AUTHOR: 
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]   
%
% VERSION NUMBER: 3.04 (10th December, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.17.1) of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_sound_speed_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_sound_speed_t_exact: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_sound_speed_t_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

n0 = 0; 
n1 = 1; 
n2 = 2;

g_tt = gsw_gibbs(n0,n2,n0,SA,t,p); 
g_tp = gsw_gibbs(n0,n1,n1,SA,t,p);

sound_speed_t_exact = gsw_gibbs(n0,n0,n1,SA,t,p) .* ...
               sqrt(g_tt./(g_tp.*g_tp - g_tt.*gsw_gibbs(n0,n0,n2,SA,t,p)));

if transposed
    sound_speed_t_exact = sound_speed_t_exact.';
end

end