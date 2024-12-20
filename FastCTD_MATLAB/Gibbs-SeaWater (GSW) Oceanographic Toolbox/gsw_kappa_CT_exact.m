function kappa_CT_exact = gsw_kappa_CT_exact(SA,CT,p)

% gsw_kappa_CT_exact                              isentropic compressibility
%==========================================================================
%
% USAGE:  
%  kappa_CT_exact = gsw_kappa_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the isentropic compressibility of seawater. 
%  
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  kappa_CT_exact  =  isentropic compressibility                   [ 1/Pa ]
%   Note. The output units are 1/Pa not 1/dbar.
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.04 (10th December, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (2.16.1) and the row for kappa in Table P.1 of appendix P  
%    of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_kappa_CT_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_kappa_CT_exact: SA and CT must have same dimensions')
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
    error('gsw_kappa_CT_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
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

t = gsw_t_from_CT(SA,CT,p);

g_tt = gsw_gibbs(n0,n2,n0,SA,t,p); 
g_tp = gsw_gibbs(n0,n1,n1,SA,t,p);

kappa_CT_exact = (g_tp.*g_tp - g_tt.*gsw_gibbs(n0,n0,n2,SA,t,p))./ ...
                  (gsw_gibbs(n0,n0,n1,SA,t,p).*g_tt);

if transposed
    kappa_CT_exact = kappa_CT_exact.';
end

end