function DataGrid = FastCTD_GridData_20240829(FCTD,varargin)
% Version of FastCTD_GridData used on TFO Seamounts. Adds uConductivity and
% fluorometer data (lines 255- )
%
%  DataGrid = FCTDMakeGrid(FCTD);
%
%  INPUTS:
%       FCTD - structure with timeseries of time, pressure, temperature,
%       etc
%
%   Temperature, Conductivity, Pressure, Salinity and Potential Density
%   will be placed on a grid of temperature versus depth
%
%  DataGrid = FCTDMakeGrid(FCTD,downcast,todo);
%   in addition to the variables described above, one can specify other
%   variable by defining TODO using a cell structure of strings
%
%  DataGrid is a structure
%
% Written by Jody Klymak
% Updated 2011 07 14 by San Nguyen
% used 20 point median filter to smooth the pressure

if ~isfield(FCTD,'pressure')
    disp(FCTD);
    DataGrid=[];
    return;
end

vars2Grid = {'pressure','temperature','conductivity','altDist','fluor','drop','longitude','latitude','chi','chi2','w','bb','chla'};
vars2Grid_default = {'pressure','temperature','conductivity','density','salinity'};
pickDownCast = true;
zInterval = 0.5;
zMin = 0;
zMax = 2000;

persistent argsNameToCheck;
if isempty(argsNameToCheck)
    argsNameToCheck = {'VarsToGrid','upcast','downcast','zMin','zMax','zInterval'};
end

index = 1;
n_items = nargin-1;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:FastCTD_GridData:wrongOption','Incorrect option specified: %s',varargin{index});
    end
    
    switch i
        case 1 % varsToGrid
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            vars2Grid = varargin{index+1};
            if iscellstr(vars2Grid)
                for i = 1:length(vars2Grid)
                    if ~isfield(FCTD,vars2Grid{i})
                        error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid{i});
                    else
                        switch lower(vars2Grid{i})
                            case vars2Grid_default
                                vars2Grid{i} = lower(vars2Grid{i});
                                continue;
                            otherwise
                                error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid{i});
                        end
                    end
                end
            elseif ischar(vars2Grid)
                switch lower(vars2Grid)
                    case vars2Grid_default
                        vars2Grid = {lower(vars2Grid)};
                        continue;
                    otherwise
                        error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid);
                end
            else
                error('MATLAB:FastCTD_GridData:wrongVar2Grid','Variable to grid must be specified as a cell strings of variables: pressure, temperature, conductivity');
            end
            
            index = index +2;
            n_items = n_items-2;
        case 2 % upcast
            pickDownCast = false;
            
            index = index + 1;
            n_items = n_items - 1;
        case 3 % downcast
            pickDownCast = true;
            
            index = index + 1;
            n_items = n_items - 1;
        case 4 % zMin
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zMin = varargin{index+1};
            if ~isnumeric(zMin)
                error('MATLAB:FastCTD_GridData:zMinNumeric','zMin must be numeric');
            elseif length(zMin) > 1 || isempty(zMin)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMin must be a scalar');
            end
            index = index +2;
            n_items = n_items-2;
        case 5 % zMax
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zMax = varargin{index+1};
            if ~isnumeric(zMax)
                error('MATLAB:FastCTD_GridData:zMaxNumeric','zMin must be numeric');
            elseif length(zMax) > 1 || isempty(zMax)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMax must be a scalar');
            end
            index = index +2;
            n_items = n_items-2;
        case 6 % zInterval
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zInterval = varargin{index+1};
            if ~isnumeric(zInterval)
                error('MATLAB:FastCTD_GridData:zMaxNumeric','zMin must be numeric');
            elseif length(zInterval) > 1 || isempty(zInterval)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMax must be a scalar');
            elseif zInterval <= 0
                error('MATLAB:FastCTD_GridData:zMinZero','zInterval must be greater than zero');
            end
            index = index +2;
            n_items = n_items-2;
    end
end

if zMin > zMax
    zTemp = zMin;
    zMin = zMax;
    zMax = zTemp;
end
clear zTemp;

if zInterval > (zMax-zMin)
    error('MATLAB:FastCTD_GridData:zIntervalTooBig','zInterval must less than the range of z');
end

% ALB DownLim for downcast 
% ALB DownLim for upcast could /should be different

downLim=.01;
upLim=.001;
if pickDownCast
    FCTD = FastCTD_FindCasts(FCTD,'downLim',downLim);
else
    FCTD = FastCTD_FindCasts(FCTD,'upcast','downLim',upLim);
end

if ~isfield(FCTD,'drop')
    DataGrid = [];
    return;
end

zMin = zMin - zInterval/2;
zMax = zMax + zInterval/2;
DataGrid.depth=(zMin:zInterval:zMax)';



% % time correction for the Conductivity cell
% df = 2.2;
% dpha = 1.05;
%
% dt = dpha/df/2/pi;
% t = FCTD.time-dt/24/3600;
% conductivity = FCTD.conductivity;
%
% good = find(FCTD.time>0);
% if isempty(good);
%     DataGrid = [];
%     return;
% end
%
% [it,ind] = unique(FCTD.time(good));
% good = good(ind);
% if length(good) > 3
%     conductivity(good) = interp1(FCTD.time(good),FCTD.conductivity(good),t(good));
% end
% FCTD.conductivity = conductivity;

drops = unique(FCTD.drop);
drops = drops(drops>0);

num = 0;
for i=1:length(drops)
    %select idx for individual cast
    ind = find(FCTD.drop==drops(i));
    % select casts with at least 10 dbar
    if max(FCTD.pressure(ind))-min(FCTD.pressure(ind))>10
        num = num+1;
    end
end

% allocate space to grid data
DataGrid.time = NaN(1,num);
for i = 1:length(vars2Grid)
    DataGrid.(vars2Grid{i}) = NaN(length(DataGrid.depth)-1,num);
end

% Add longitude and latitude
if isfield(FCTD,'GPS')
    FCTD.longitude = FCTD.GPS.longitude;
    FCTD.latitude = FCTD.GPS.latitude;
end

% load correction factors for FCTD
% ALB how are these Factors computed
FCTD_SalCorr = load('FCTD_SalinityCorrectionFactors_toCond.mat');
% FCTD.depth = sw_dpth(FCTD.pressure,median(FCTD.GPS.latitude));
% FCTD.depth = sw_dpth(FCTD.pressure,20);
% for i = 1:length(vars2Grid)
%     FCTD.(vars2Grid{i}) = nanmedfilt1(FCTD.(vars2Grid{i}),10);
% end
FCTD_SalCorr.GainPFit = FCTD_SalCorr.GainPFit_Dn;
FCTD_SalCorr.PhsPFit = FCTD_SalCorr.PhsPFit_Dn;

% Loop through the drops and concatenate them
num = 0;
for i=1:length(drops)
    ind = find(FCTD.drop==drops(i));
    if max(FCTD.pressure(ind))-min(FCTD.pressure(ind))>10
        for j=1:length(vars2Grid)
            %ALB add try/catch because it break in the fluor during TFO seamount 2023 
            try
                myFCTD.(vars2Grid{j}) = FCTD.(vars2Grid{j})(ind);
            catch
                myFCTD.(vars2Grid{j})=nan.*ind;
            end
            % ALB end of hack
        end
        %         myFCTD.temperature = FCTD.temperature(ind);
        %         myFCTD.pressure = FCTD.pressure(ind);
        %         myFCTD.conductivity = FCTD.conductivity(ind);
        
        % do salinity despiking corrections
        good_ind = ~isnan(myFCTD.pressure);
        % try
        % test_salinity_spiking;
        % catch
        % end
        npts = sum(good_ind);
        df = FCTD_SalCorr.f_Ny/floor(npts/2);
        % creating the frequency axis
        myFCTD.f = (0:npts-1)'*df;
        myFCTD.f(myFCTD.f>FCTD_SalCorr.f_Ny) = myFCTD.f(myFCTD.f>FCTD_SalCorr.f_Ny)-2*FCTD_SalCorr.f_Ny;
        FCTD_SalCorr.GainFit = polyval(FCTD_SalCorr.GainPFit,myFCTD.f);
        FCTD_SalCorr.GainFit = FCTD_SalCorr.GainFit/FCTD_SalCorr.GainFit(1); % need to normalize Gain
        FCTD_SalCorr.PhsFit = polyval(FCTD_SalCorr.PhsPFit,myFCTD.f);
        
        % FFT
        myFCTD.T = fft(myFCTD.temperature(good_ind),npts);
        myFCTD.C = fft(myFCTD.conductivity(good_ind),npts);
        myFCTD.P = fft(myFCTD.pressure(good_ind),npts);
        
        % correct the conductivity
        myFCTD.CCorr = myFCTD.C.*FCTD_SalCorr.GainFit.*exp(-1i*FCTD_SalCorr.PhsFit);
        % Low Pass filter
        myFCTD.TCorr = myFCTD.T.*FCTD_SalCorr.LPfilter(myFCTD.f);
        myFCTD.CCorr = myFCTD.CCorr.*FCTD_SalCorr.LPfilter(myFCTD.f);
        myFCTD.PCorr = myFCTD.P.*FCTD_SalCorr.LPfilter(myFCTD.f);
        
        % get back to physical units
        myFCTD.tCorr = real(ifft(myFCTD.TCorr));
        myFCTD.cCorr = real(ifft(myFCTD.CCorr));
        myFCTD.pCorr = real(ifft(myFCTD.PCorr));
        
        % patching up the pressure because the pressure doesn't have a
        % phaseshift
        myFCTD.pCorr(1:13) = NaN;%myFCTD.pressure(1:13);
        myFCTD.pCorr(end-12:end) = NaN;%myFCTD.pressure(end-12:end);
        
        myFCTD.pressure = myFCTD.pCorr;
        myFCTD.temperature = myFCTD.tCorr;
        myFCTD.conductivity = myFCTD.cCorr;
        myFCTD.depth = sw_dpth(myFCTD.pressure,mean(myFCTD.latitude,'omitmissing'));
        
        num = num+1;
        DataGrid.time(num) = mean(FCTD.time(ind),'omitmissing');
        for j=1:length(vars2Grid)
            DataGrid.(vars2Grid{j})(:,num) = bindata1d(DataGrid.depth,...
                myFCTD.depth, ...
                myFCTD.(vars2Grid{j}));
        end
    end
end

DataGrid.depth = midpoints(DataGrid.depth);

if ~isfield(DataGrid,'temperature') || ~isfield(DataGrid,'pressure') || ~isfield(DataGrid,'conductivity')
    return;
end

DataGrid.salinity = sw_salt(DataGrid.conductivity*10/sw_c3515,DataGrid.temperature,DataGrid.pressure);
DataGrid.density = sw_pden(DataGrid.salinity,DataGrid.temperature,DataGrid.pressure,0);

% gridding in time
mintime = min(DataGrid.time);
maxtime = max(DataGrid.time);

DataGrid.tGrid.time = mintime:2*median(diff(DataGrid.time),'omitmissing'):maxtime; % every minute
DataGrid.tGrid.depth = DataGrid.depth;


% allocate space to grid data
for i = 1:length(vars2Grid)
    DataGrid.tGrid.(vars2Grid{i}) = NaN(length(DataGrid.tGrid.depth),length(DataGrid.tGrid.time)-1);
end

for i = 1:length(DataGrid.tGrid.depth)
    for j = 1:length(vars2Grid)
        DataGrid.tGrid.(vars2Grid{j})(i,:) = bindata1d(DataGrid.tGrid.time,...
            DataGrid.time, DataGrid.(vars2Grid{j})(i,:));
    end
end

DataGrid.tGrid.salinity = sw_salt(DataGrid.tGrid.conductivity*10/sw_c3515,DataGrid.tGrid.temperature,DataGrid.tGrid.pressure);
DataGrid.tGrid.density = sw_pden(DataGrid.tGrid.salinity,DataGrid.tGrid.temperature,DataGrid.tGrid.pressure,0);

DataGrid.tGrid.time = midpoints(DataGrid.tGrid.time);

return;
end

function FCTD = FastCTD_FindCasts(FCTD,varargin)
% FCTD = FastCTD_FindCasts(FCTD);
%   finds when the FastCTD is going up and when it is going down
%   if 'UPCAST' is specified then the profile will search for upcast but by
%   defalt it would search for downcasts
%
%   In this data set we assume the data is in descending order of time
%
% Written by Jody Klymak
% Updated 2011 07 14 by San Nguyen
% Updated 2012 09 29 by San Nguyen



if ~isfield(FCTD,'pressure')
    return;
end

% use 20 point median filter to smooth out the pressure field
p = medfilt1(FCTD.pressure,256);
% try to smooth out the data a bit
dp = conv2(diff(conv2(p,ones(256,1)/256,'same'),1,1)',ones(1,256)/256,'same');


% ALB super important coef that should be defined in above levels (e.g., )
%downLim = 0.1;
downLim = 0.025;
downCast = true;

persistent argsNameToCheck;
if isempty(argsNameToCheck)
    argsNameToCheck = {'downLim','threshold','upcast','downcast'};
end

index = 1;
n_items = nargin-1;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:FastCTD_FindCasts:wrongOption','Incorrect option specified: %s',varargin{index});
    end
    
    switch i
        case 1 % downLim
            if n_items == 1
                error('MATLAB:FastCTD_FindCasts:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0
                error('MATLAB:FastCTD_FindCasts:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 2 % downLim
            if n_items == 1
                error('MATLAB:FastCTD_FindCasts:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0
                error('MATLAB:FastCTD_FindCasts:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 3 % upcast
            downCast = false;
            
            index = index + 1;
            n_items = n_items - 1;
        case 4 % downcast
            downCast = true;
            
            index = index + 1;
            n_items = n_items - 1;
    end
end

%defining the threshold for going up and down
if downCast && downLim < 0 % case where we want downcast but downlim is set to a negative speed through water.
    downLim = -downLim;
elseif (~downCast) && downLim > 0% going up
    downLim = -downLim;
end

if downCast
    dn = find(dp>downLim);
else
    dn = find(dp<downLim);
end

% find all indices of going down

if isempty(dn)
    return;
end

dn = [0, dn];


% find jumps in indices to indicate a start of a profile
startdown = dn(find(diff(dn)>1)+1);

if isempty(startdown)
    return;
end

dn = dn(2:end);
FCTD.drop = 0*FCTD.time;

if dn(1)<startdown(1)
    startdown=[dn(1) startdown];
end

if startdown(end)<dn(end)
    startdown = [startdown dn(end)];
end


for i=1:(length(startdown)-1)
    in = intersect(startdown(i):startdown(i+1)-1,dn);
    FCTD.drop(in) = i;
end
end