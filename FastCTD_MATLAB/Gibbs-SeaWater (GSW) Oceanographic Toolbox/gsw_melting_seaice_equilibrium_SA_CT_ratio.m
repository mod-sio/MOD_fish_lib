function melting_seaice_equilibrium_SA_CT_ratio = gsw_melting_seaice_equilibrium_SA_CT_ratio(SA,p,saturation_fraction)

% gsw_melting_seaice_equilibrium_SA_CT_ratio      ratio of SA to CT changes 
%                                   when sea ice melts into a large mass of 
%                                      seawater, with both the seawater and 
%                                   sea ice temperatures being almost equal
%                                   to the equilibrium freezing temperature                                   
%==========================================================================
%
% USAGE:
%  melting_seaice_equilibrium_SA_CT_ratio = ...
%      gsw_melting_seaice_equilibrium_SA_CT_ratio(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when sea ice melts into 
%  seawater with both the seawater and the sea ice temperatures being  
%  almost equal to the equilibrium freezing temperature.  It is assumed  
%  that a small mass of seaice melts into an infinite mass of seawater.  If 
%  indeed the temperature of the seawater and the sea ice were both equal  
%  to the freezing temperature, then no melting or freezing would occur; an  
%  imbalance between these three temperatures is needed for freezing or 
%  melting to occur (the three temperatures being (1) the seawater 
%  temperature, (2) the sea ice temperature, and (3) the freezing 
%  temperature.  
%
%  Note that the output of this function, dSA/dCT is independent of the 
%  sea ice salinity, SA_seaice.  That is, the output applies equally to  
%  pure ice Ih and to sea ice with seaice salinity, SA_seaice.  This result 
%  is proven in the manuscript, McDougall et al. (2013).  
%
%  The output, melting_seaice_equilibrium_SA_CT_ratio, is dSA/dCT rather  
%  than dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is 
%  zero whereas dCT/dSA would be infinite. 
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%  saturation_fraction = the saturation fraction of dissolved air in 
%               seawater.  The saturation_fraction must be between 0 and 1.
%
% p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  melting_seaice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT  
%                            changes when sea ice melts into seawater, with   
%                            the seawater and sea ice being close to the  
%                            freezing temperature.             [ g/(kg K) ]               
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.04 (6th December, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., P.M. Barker and R. Feistel, 2013: Melting of ice and 
%   sea ice into seawater and frazil ice formation. Journal of Physical 
%   Oceanography, (Submitted).
%    See Eqn. (29) of this manuscript.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3) 
   error('gsw_melting_seaice_equilibrium_SA_CT_ratio: Requires three inputs')
end 

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_melting_seaice_equilibrium_SA_CT_ratio: saturation fraction MUST be between zero and one.')
end

[ms,ns] = size(SA);
[mp,np] = size(p);
[msf,nsf] = size(saturation_fraction);

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                              % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_melting_seaice_equilibrium_SA_CT_ratio: Inputs array dimensions arguments do not agree; check p')
end 

if (msf == 1) & (nsf == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(SA));         % fill to size of SA
elseif (ns == nsf) & (msf == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (nsf == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,ns));        % copy across each row.
elseif (ns == msf) & (nsf == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (ns == nsf)
    % ok
else
    error('gsw_melting_seaice_equilibrium_SA_CT_ratio: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

CTf = gsw_CT_freezing(SA,p,saturation_fraction);
t_seaice = gsw_t_freezing(SA,p,saturation_fraction);

h = gsw_enthalpy_CT_exact(SA,CTf,p);
h_Ih = gsw_enthalpy_ice(t_seaice,p);
[h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CTf,p);
          % Note that h_hat_CT is equal to cp0*(273.15 + t)./(273.15 + pt0)

denominator = h - h_Ih - SA.*h_hat_SA;
melting_seaice_equilibrium_SA_CT_ratio = SA.*h_hat_CT./denominator;
         
if transposed
    melting_seaice_equilibrium_SA_CT_ratio = melting_seaice_equilibrium_SA_CT_ratio.';
end

end