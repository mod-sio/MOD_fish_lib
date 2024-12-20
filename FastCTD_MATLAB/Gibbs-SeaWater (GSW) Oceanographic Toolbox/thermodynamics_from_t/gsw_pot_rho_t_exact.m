function pot_rho_t_exact = gsw_pot_rho_t_exact(SA,t,p,p_ref)

% gsw_pot_rho_t_exact                                     potential density
%==========================================================================
%
% USAGE:
%  pot_rho_t_exact = gsw_pot_rho_t_exact(SA,t,p,p_ref)
%
% DESCRIPTION:
%  Calculates potential density of seawater.  Note. This function outputs
%  potential density, not potential density anomaly; that is, 1000 kg/m^3
%  is not subtracted.  
%
% INPUT:
%  SA     =  Absolute Salinity                                     [ g/kg ]
%  t      =  in-situ temperature (ITS-90)                         [ deg C ]
%  p      =  sea pressure                                          [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref  =  reference pressure                                    [ dbar ]
%            ( i.e. reference absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p & p_ref may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t 
%  are MxN
%
% OUTPUT:
%  pot_rho_t_exact  =  potential density (not potential density anomaly)
%                                                                [ kg/m^3 ]
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.04 (10th December, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.4 of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4 )
   error('gsw_pot_rho_t_exact:  Requires four inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_pot_rho_t_exact: SA and t must have same dimensions')
end

if ~isscalar(unique(p_ref))
    error('gsw_pot_rho_t_exact: The reference pressures differ, they should be unique')
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
    error('gsw_pot_rho_t_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

upr = unique(p_ref);
p_ref = upr*ones(size(SA));

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pt = gsw_pt_from_t(SA,t,p,p_ref);

pot_rho_t_exact = gsw_rho_t_exact(SA,pt,p_ref);

if transposed
    pot_rho_t_exact = pot_rho_t_exact.';
end

end
