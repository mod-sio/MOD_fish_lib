function [SA, in_ocean] = gsw_SA_from_Sstar(Sstar,p,long,lat)

% gsw_SA_from_Sstar               Absolute Salinity from Preformed Salinity 
%==========================================================================
%
% USAGE:
%  [SA, in_ocean] = gsw_SA_from_Sstar(Sstar,p,long,lat)
%
% DESCRIPTION:
%  Calculates Absolute Salinity from Preformed Salinity. 
%
% INPUT:
%  Sstar  =  Preformed Salinity                                    [ g/kg ]
%  p      =  sea pressure                                          [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  long   =  longitude in decimal degrees                    [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat    =  latitude in decimal degrees north              [ -90 ... +90 ]
%
%  p, lat and long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where Sstar is MxN.
%
% OUTPUT:
%  SA        =  Absolute Salinity                                  [ g/kg ]
%  in_ocean  =  0, if long and lat are a long way from the ocean 
%            =  1, if long and lat are in the ocean
%  Note. This flag is only set when the observation is well and truly on
%    dry land; often the warning flag is not set until one is several 
%    hundred kilometres inland from the coast. 
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
%
%  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and 
%   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
%   Ocean Science, 8, 1123-1134.  
%   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_SA_from_Sstar:  Requires four inputs')
end %if

[ms,ns] = size(Sstar);
[mp,np] = size(p);

if (mp == 1) & (np == 1)              % p scalar - fill to size of Sstar
    p = p*ones(size(Sstar));
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
    error('gsw_SA_from_Sstar: Inputs array dimensions arguments do not agree')
end %if

[mla,nla] = size(lat);

if (mla == 1) & (nla == 1)             % lat is a scalar - fill to size of Sstar
    lat = lat*ones(size(Sstar));
elseif (ns == nla) & (mla == 1)        % lat is a row vector,
    lat = lat(ones(1,ms), :);          % copy down each column.
elseif (ms == mla) & (nla == 1)        % lat is a column vector,
    lat = lat(:,ones(1,ns));           % copy across each row.
elseif (ns == mla) & (nla == 1)        % lat is a transposed row vector,
    lat = lat.';                         % transposed then
    lat = lat(ones(1,ms), :);           % copy down each column.
elseif (ms == mla) & (ns == nla)
    % ok
else
    error('gsw_SA_from_Sstar: Inputs array dimensions arguments do not agree')
end %if

[mlo,nlo] = size(long);
long(long < 0) = long(long < 0) + 360; 

if (mlo == 1) & (nlo == 1)            % long is a scalar - fill to size of Sstar
    long = long*ones(size(Sstar));
elseif (ns == nlo) & (mlo == 1)       % long is a row vector,
    long = long(ones(1,ms), :);       % copy down each column.
elseif (ms == mlo) & (nlo == 1)       % long is a column vector,
    long = long(:,ones(1,ns));        % copy across each row.
elseif (ns == mlo) & (nlo == 1)       % long is a transposed row vector,
    long = long.';                      % transposed then
    long = long(ones(1,ms), :);        % copy down each column.
elseif (ms == mlo) & (ns == nlo)
    % ok
else
    error('gsw_SA_from_Sstar: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    Sstar = Sstar.';
    p = p.';
    lat = lat.';
    long = long.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

r_1 = 0.35;

[Iocean] = find(~isnan(Sstar + p + lat + long));

SA = nan(size(Sstar));
SAAR = SA;
in_ocean = SA;

[SAAR(Iocean), in_ocean(Iocean)] = gsw_SAAR(p(Iocean),long(Iocean),lat(Iocean));

SA(Iocean) = Sstar(Iocean).*(1 + SAAR(Iocean))./(1 - r_1*SAAR(Iocean));

% In the Baltic Sea, SA = Sstar, and note that gsw_delta_SA returns zero
% for dSA in the Baltic.

if transposed
    SA = SA.';
    in_ocean = in_ocean.';
end

end
