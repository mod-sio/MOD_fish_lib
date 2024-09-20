% merge N grid the FCTD for miso-bob NOW FOR 2019 (based on 2018) !!!!!
% this the main file that I used to go through and process the FCTD for the
% grids and figures used in the script drews_mbob_log.m
% the files I produced: FCTD_grid.mat (all downcasts), FCTD_gridup.mat (upcasts
% from only the FCTDF period) and FCTD_all.mat (full timeseries, ungridded)
% are saved locally at ~/miso-bob/cruise/FCTD/MAT/ also see mbob_fctd
% curtain.


%% fluorometry
% after the addition of ECOPUCK (now with different ecopuck call, paths for
% the Sally Ride)

matDir = ('/Users/Shared/FCTD/MAT/');

% clear FastCTD_MATfile_TimeIndex_all

% filenames=FastCTD_MATfile_TimeIndex.filenames;
% timeStart=FastCTD_MATfile_TimeIndex.timeStart;
% timeEnd=FastCTD_MATfile_TimeIndex.timeEnd;
%
% FastCTD_MATfile_TimeIndex_all.filenames=filenames;
% FastCTD_MATfile_TimeIndex_all.timeStart=timeStart;
% FastCTD_MATfile_TimeIndex_all.timeEnd=timeEnd;


% FastCTD_MATfile_TimeIndex=FastCTD_MATfile_TimeIndex_all;
clear FCTD

timeStart = datenum([2019 7 08 00 00 00]);
timeEnd = now;

% we spaced by 59 mins to calculate hours is to ensure that there is no
% rounding errors to find all hours between timestart and timeend
time_line_59_mins_space = timeStart:59/60/24:timeEnd;

time_line_59_mins_space_dv = datevec(time_line_59_mins_space);

time_line_59_mins_space_dv(:,5:6) = 0;

time_line_1_hr_space = datenum(time_line_59_mins_space_dv);

time_line_1_hr_space = unique(time_line_1_hr_space);

%%

for i_hour = 1:numel(time_line_1_hr_space)
    time_now = datevec(now);
    time_now(5:6) = 0;
    time_now = datenum(time_now);
    
    new_fname = datestr(time_line_1_hr_space(i_hour),'yyyymmdd_HHMMSS');
    
    if(exist(fullfile(matDir,'_hourly','combined',['FCTD_' new_fname '.mat']),'file') && ...
            exist(fullfile(matDir,'_hourly','gridded',['FCTD_' new_fname '.mat']),'file') &&...
            exist(fullfile(matDir,'_hourly','gridded','upcasts',['FCTD_' new_fname '.mat']),'file'))
        % make sure to reprocess the last two hours
        if(time_now - time_line_1_hr_space(i_hour))> 2/24
            continue;
        end
    end
    disp(new_fname);
    time_start = time_line_1_hr_space(i_hour);
    time_end = time_line_1_hr_space(i_hour) + 1/24;
    
    load(fullfile(matDir, 'FastCTD_MATfile_TimeIndex.mat'));
    
    myFCTD = [];
    
    % looking at 15 mins back and 15 mins forward
    % might want to change this for different sampling scheme
    ind = find(FastCTD_MATfile_TimeIndex.timeStart < (time_end+15/60/24) &...
        FastCTD_MATfile_TimeIndex.timeEnd >= (time_start - 15/60/24));
    if(isempty(ind))
        continue;
    end
    
    for i = 1:length(ind)
        %load data
        load([matDir '/' FastCTD_MATfile_TimeIndex.filenames{ind(i)} '.mat']);
        if exist('FCTD','var')
            if isstruct(FCTD) && ~isempty(FCTD.time) && (FCTD.time(end)>=timeStart && FCTD.time(1) <= timeEnd)
                myFCTD = FastCTD_MergeFCTD(myFCTD,FCTD);
            end
            clear FCTD;
        end
    end
    
    if(isempty(myFCTD))
        return;
    end
    myFCTD.bs=myFCTD.eco_triplet(:,1);
    myFCTD.chla=myFCTD.eco_triplet(:,2);
    myFCTD.cdom=myFCTD.eco_triplet(:,3);
    
    %smooth the chla and bs a bit to beat down some noise
    
    %keyboard
    chla=medfilt1(myFCTD.chla,16,'omitnan');
    chlaN=naninterp(chla);
    % chlaN=naninterp(myFCTD.chla);
    % chlaN(1)=chlaN(2);
    % chlaN(end)=chlaN(end-1);
    % [b,a]=butter(6,0.5/8); %2 second low pass
    % chlaS=filtfilt(b,a,chlaN);
    myFCTD.chla=chlaN;
    
    bs=medfilt1(myFCTD.bs,16,'omitnan');
    bsN=naninterp(bs);
    % bsN=naninterp(myFCTD.bs);
    % bsN(1)=bsN(2);
    % bsN(end)=bsN(end-1);
    % [b,a]=butter(6,0.5/8); %2 second low pass
    % bsS=filtfilt(b,a,bsN);
    myFCTD.bs=bsN(:);
    
    cdom=medfilt1(myFCTD.cdom,16,'omitnan');
    cdomN=naninterp(cdom);
    % bsN=naninterp(myFCTD.bs);
    % bsN(1)=bsN(2);
    % bsN(end)=bsN(end-1);
    % [b,a]=butter(6,0.5/8); %2 second low pass
    % bsS=filtfilt(b,a,bsN);
    myFCTD.cdom=cdom;
    
    FCTD_GridData = FastCTD_GridDataF(myFCTD,'zMin',0,'zMax',210,'zInterval',1,'downcast');
    FCTD_GridData_up = FastCTD_GridDataF(myFCTD,'zMin',0,'zMax',210,'zInterval',1,'upcast'); %check upcast for fluorometry
    
    %calc N2
    [m,n]=size(FCTD_GridData.temperature);
    FCTD_GridData.N2=nan(m,n);
    
    for kk=1:n
        id=find(~isnan(FCTD_GridData.density(:,kk)));
        if(numel(id)<4)
            continue;
        end
        dz=ddz(FCTD_GridData.depth(id));
        FCTD_GridData.N2(id,kk)=9.8/1024.*(dz*FCTD_GridData.density(id,kk));
    end
    
    [m,~]=size(FCTD_GridData.temperature);
    FCTD_GridData.lat=interp1(myFCTD.GPS.GPS_time,myFCTD.GPS.latitude,FCTD_GridData.time);
    FCTD_GridData.lat=repmat(FCTD_GridData.lat,m,1);
    FCTD_GridData.lon=interp1(myFCTD.GPS.GPS_time,myFCTD.GPS.longitude,FCTD_GridData.time);
    FCTD_GridData.lon=repmat(FCTD_GridData.lon,m,1);
    
    [m,~]=size(FCTD_GridData_up.temperature);
    FCTD_GridData_up.lat=interp1(myFCTD.GPS.GPS_time,myFCTD.GPS.latitude,FCTD_GridData_up.time);
    FCTD_GridData_up.lat=repmat(FCTD_GridData_up.lat,m,1);
    FCTD_GridData_up.lon=interp1(myFCTD.GPS.GPS_time,myFCTD.GPS.longitude,FCTD_GridData_up.time);
    FCTD_GridData_up.lon=repmat(FCTD_GridData_up.lon,m,1);
    
    
    FCTD = rmfield(myFCTD,'winch');
    fctd_fieldnames = fieldnames(FCTD);
    ind = FCTD.time < time_end & FCTD.time >= time_start;
    for i = 1:numel(fctd_fieldnames)
        if(~isempty(strfind(fctd_fieldnames{i},'header')))
            continue;
        elseif (~isempty(strfind(fctd_fieldnames{i},'GPS')))
            continue;
        end
        FCTD.(fctd_fieldnames{i}) = FCTD.(fctd_fieldnames{i})(ind,:);
    end
    
    fctd_fieldnames = fieldnames(FCTD.GPS);
    ind = FCTD.GPS.time < time_end & FCTD.GPS.time >= time_start;
    for i = 1:numel(fctd_fieldnames)
        FCTD.GPS.(fctd_fieldnames{i}) = FCTD.GPS.(fctd_fieldnames{i})(ind,:);
    end
    
    save(fullfile(matDir,'_hourly','combined',['FCTD_' new_fname '.mat']),'FCTD','-v7.3');
    
    FCTD=FCTD_GridData;
    
    fctd_fieldnames = fieldnames(FCTD);
    ind = FCTD.time < time_end & FCTD.time >= time_start;
    for i = 1:numel(fctd_fieldnames)
        if(~isempty(strfind(fctd_fieldnames{i},'depth')))
            continue;
        elseif (~isempty(strfind(fctd_fieldnames{i},'tGrid')))
            continue;
        end
        FCTD.(fctd_fieldnames{i}) = FCTD.(fctd_fieldnames{i})(:,ind);
    end
    
    fctd_fieldnames = fieldnames(FCTD.tGrid);
    ind = FCTD.tGrid.time < time_end & FCTD.tGrid.time >= time_start;
    for i = 1:numel(fctd_fieldnames)
        if(~isempty(strfind(fctd_fieldnames{i},'depth')))
            continue;
        end
        FCTD.tGrid.(fctd_fieldnames{i}) = FCTD.tGrid.(fctd_fieldnames{i})(:,ind);
    end
    
    save(fullfile(matDir,'_hourly','gridded',['FCTD_' new_fname '.mat']),'FCTD','-v7.3');
    
    FCTD=FCTD_GridData_up;
    
    fctd_fieldnames = fieldnames(FCTD);
    ind = FCTD.time < time_end & FCTD.time >= time_start;
    for i = 1:numel(fctd_fieldnames)
        if(~isempty(strfind(fctd_fieldnames{i},'depth')))
            continue;
        elseif (~isempty(strfind(fctd_fieldnames{i},'tGrid')))
            continue;
        end
        FCTD.(fctd_fieldnames{i}) = FCTD.(fctd_fieldnames{i})(:,ind);
    end
    
    fctd_fieldnames = fieldnames(FCTD.tGrid);
    ind = FCTD.tGrid.time < time_end & FCTD.tGrid.time >= time_start;
    for i = 1:numel(fctd_fieldnames)
        if(~isempty(strfind(fctd_fieldnames{i},'depth')))
            continue;
        end
        FCTD.tGrid.(fctd_fieldnames{i}) = FCTD.tGrid.(fctd_fieldnames{i})(:,ind);
    end
    
    save(fullfile(matDir,'_hourly','gridded','upcasts',['FCTD_' new_fname '.mat']),'FCTD','-v7.3');
    
    clear('myFCTD','FCTD_GridData','FCTD_GridData_up','FastCTD_MATfile_TimeIndex');
end
disp([datestr(now) ': Done ' mfilename '!']);
% clear FastCTD_MATfile_TimeIndex

% save([matDir '/FCTD_all.mat'],'FCTD_all','-v7.3')
% save([matDir '/FCTD_grid.mat'],'FCTD_grid','-v7.3')

%%
% load d1 to concatenate gridded downcasts, leave everything else alone for
% now DO NOT USE
%
% t=load('/Volumes/cruiseshare/data/FCTD/MAT/deployment_1/FCTD_grid_d1.mat');
% FCTD_grid=mergefields(t.FCTD_grid,FCTD_grid);
% FCTD_grid.depth=FCTD_grid.depth(:,1);
% save([matDir '/FCTD_grid.mat'],'FCTD_grid','-v7.3')

