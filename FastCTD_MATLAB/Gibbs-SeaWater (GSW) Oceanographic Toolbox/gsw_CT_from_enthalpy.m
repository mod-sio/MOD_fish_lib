function CT = gsw_CT_from_enthalpy(SA,h,p)

% gsw_CT_from_enthalpy               Conservative Temperature from specific
%                                   enthalpy of seawater (48-term equation)
%==========================================================================
%
% USAGE:
%  CT = gsw_CT_from_enthalpy(SA,h,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of seawater, given the Absolute 
%  Salinity, specific enthalpy, h, and pressure p.  The specific enthalpy 
%  input is the one calculated from the computationally-efficient 48-term
%  expression for specific volume in terms of SA, CT and p (IOC et al., 
%  2010).
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  h   =  specific enthalpy                                        [ J/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & h need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & h are MxN.
%
% OUTPUT:
%  CT  =  Conservative Temperature ( ITS-90)                      [ deg C ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.04 (10th December, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T. J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_CT_from_enthalpy: requires three inputs')
end

[ms,ns] = size(SA);
[mh,nh] = size(h); 
[mp,np] = size(p);

if (mh ~= ms | nh ~= ns)
    error('gsw_CT_from_enthalpy: SA and h must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(ms,ns);
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
    error('gsw_CT_from_enthalpy: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    h = h.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensures that SA is non-negative

CT_freezing = gsw_CT_freezing(SA,p,0); % This is the CT freezing temperature
h_freezing = gsw_enthalpy(SA,CT_freezing,p);
if (h < h_freezing)
    error('gsw_CT_from_enthalpy: The seawater enthalpy h that was inputed is less than the enthalpy at the freezing temperature, i.e. the water is frozen')
end 

h_40 = gsw_enthalpy(SA,40*ones(size(SA)),p);
if (h > h_40)
    error('gsw_CT_from_enthalpy: The seawater enthalpy h that was inputed is greater than the enthalpy when CT is 40 C')
end

CT = CT_freezing + (40 - CT_freezing).*(h - h_freezing)./(h_40 - h_freezing); % First guess of CT
[dummy, h_CT] = gsw_enthalpy_first_derivatives(SA,CT,p);

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------

CT_old = CT;
f = gsw_enthalpy(SA,CT_old,p) - h;
CT = CT_old - f./h_CT ; % this is half way through the modified Newton's method (McDougall and Wotherspoon, 2013)
CT_mean = 0.5*(CT + CT_old);
[dummy, h_CT] = gsw_enthalpy_first_derivatives(SA,CT_mean,p);
CT = CT_old - f./h_CT ; % this is the end of one full iteration of the modified Newton's method

CT_old = CT;
f = gsw_enthalpy(SA,CT_old,p) - h;
CT = CT_old - f./h_CT ; % this is half way through the modified Newton's method (McDougall and Wotherspoon, 2013)

% After 1.5 iterations of this modified Newton-Raphson iteration,
% the error in CT is no larger than 4x10^-13 degrees C, which 
% is machine precision for this calculation. 

if transposed
    CT = CT.';
end

end
