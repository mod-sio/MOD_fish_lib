% RUN_epsiAuto_process_data.m
%
% This is a wrapper script to set up input values and run
% epsiAuto_process_data.m
%
% CREATE:
% input_struct, a structure containing the following fields
%
% Required:
%  .raw_dir = path to directory where the data are streaming in
%  .process_dir = path to directory where data will be copied and where
%              subdirectories raw, mat, and FCTDmat will be created
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .str_to_match = (default '*'), Name (or partial name) of first file in raw_dir to copy
%  .refresh_time_sec = (default 5*60), refresh period in seconds
%  .version       = (default 4), version of mod_som_read_epsi_files.m to use
% -------------------------------------------------------------------------
%
% USER INPUTS
% These will probably be the same for the whole cruise
clear input_struct
input_struct.Meta_Data_process_file = '/Volumes/Software_current_cruise/MOD_fish_lib/EPSILOMETER/Meta_Data_Process/MDP_motive_2024.txt';
input_struct.refresh_time_sec =  2*60;
input_struct.cruise_specifics = 'tfo_2024';
epsi_depth_array = 0:1600;
fctd_depth_array = 0:2100;

% Realtime or Simulator mode
data_mode = 'realtime'; %'realtime' or 'simulator'

% -------------------------------------------------------------------------



%% Find the fish name and the path to store processed data
% 1. If you're running this in realtime mode, the path to processed data will
% be in the most recent Setup file.
%
% 2. If you're running this is simulated mode, the most recent Setup file is
% probably not what you want or the path to it might not exist (because
% you're running this on a personal computer. In that case, ask for the
% directory where you want to store processed data.

% First, look for the instrument name and the survey name from the most
% recent .modraw file. Create the output directory based on the survey
% name. If it doesn't exist in the file, ask the user to define it. If the
% instrument name doesn't exist in the file, ask the user to define it.

% Then look for any (up to 10) .modraw files you might have missed that
% have the same



% Get the
switch data_mode

    case 'simulator'

        % Define raw directory for simulated data
        input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/Simulated_Data/Realtime_RAW/raw/';

        % Path to setup
        path2setup = '/Volumes/MOD HD/Users/Shared/FCTD_EPSI_DATA/Simulated_Data/Setup';

        process_dir_root = '/Users/Shared/EPSI_PROCESSING/Simulated_Data/Processed';

    case 'realtime'

        % Define raw directory for realtime data
        input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Realtime_RAW/raw/';

        % Path to setup file
        root_software='/Volumes/MOD HD/Users/Shared/Software_current_cruise/MOD_fish_lib/';
        path2setup=fullfile(root_software,'Acquisition/fctd_epsi_acq/build/fctd_epsi/Build/Products/Debug/Setup');

        process_dir_root = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed';
end

%% Read the Setup file
% Look for fish flag in Setup file
fid=fopen(path2setup,'r');
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);
newSetup_flag=contains(str,'CTD.fishflag=');
if newSetup_flag
    fishflag_str      = str(strfind(str,'CTD.fishflag=')+(0:100));
    fishflag_str      = fishflag_str(1:find(uint8(fishflag_str)==10,1,'first'));
    fishflag_name      = strsplit(fishflag_str,'=');
    fishflag_name      = fishflag_name{2}(2:end-2);
    instrument = fishflag_name;

else
    instrument = input('What fish are we using? [''epsi'',''fctd'']');

end

% Look for survey name in Setup file
newSurvey_flag=contains(str,'CTD.survey');
if newSetup_flag
    surveyflag_str      = str(strfind(str,'CTD.survey')+(0:100));
    surveyflag_str      = surveyflag_str(1:find(uint8(surveyflag_str)==10,1,'first'));
    surveyflag_name     = strsplit(surveyflag_str,'=');
    surveylag_name      = surveyflag_name{2}(1:end-1);
    survey_name = surveylag_name;

else
    survey_name = input('What is the survey name? [''yyyymmdd_d##_NAME'']');

end

% Define the name of the process directory based on either what you
% found in the Setup file or what you input
input_struct.process_dir = fullfile( ...
    process_dir_root, ...
    strrep(survey_name,'''','')); %This will create a directory with this name

%% Look for the files that match the survey name
N_files_to_search = 10; %Maximum number of files to search for CTD.survey
if exist(fullfile(input_struct.process_dir,'raw'),'dir')
    listfile_process_dir=dir(fullfile(input_struct.process_dir,'raw','*.modraw'));
    % If the processed data directory is not empty, that means you've
    % already started copying files, so 'str_to_match' is the last file in
    % the directory. Copy everything from the raw incoming directory
    % beginning with that file.
    if ~isempty(listfile_process_dir)
        input_struct.str_to_match = listfile_process_dir(end).name;
    else
        % If the processed data directory is empty, look for files in the
        % raw directory. Look for what you defined as 'survey_name' in the
        % last 10 files and start copying from the first one with that
        % survey name.
        

        listfile_raw_dir=dir(fullfile(input_struct.raw_dir,'*.modraw'));
        count=0;
        % quick look in the previous files if survey name exist.
        % In case blue matlab was launched few minutes after fctd_epsi
        
        while count<min([N_files_to_search,length(listfile_raw_dir)]) %Count up to the total number of files in the directory or N_files_to_search, whatever is smaller
            path2setup1=fullfile(...
                listfile_raw_dir((length(listfile_raw_dir)-count)).folder,...
                listfile_raw_dir((length(listfile_raw_dir)-count)).name);
            fid=fopen(path2setup1,'r');
            fseek(fid,0,1);
            frewind(fid);
            str = fread(fid,'*char')';
            fclose(fid);
            % If you find 'CTD.survey' in the .mod raw file, get the survey
            % name
            if ~isempty(strfind(str,'CTD.survey'))
                surveyflag1_str      = str(strfind(str,'CTD.survey')+(0:100));
                surveyflag1_str      = surveyflag1_str(1:find(uint8(surveyflag1_str)==10,1,'first'));
                surveyflag1_name     = strsplit(surveyflag1_str,'=');
                surveylag1_name      = surveyflag1_name{2}(1:end-1);
                survey1_name         = surveylag1_name;

                if contains(survey1_name,survey_name)
                    input_struct.str_to_match = ...
                        listfile_raw_dir(length(listfile_raw_dir)-count).name;
                    count=count+1;
                else
                    input_struct.str_to_match = ...
                        listfile_raw_dir(length(listfile_raw_dir)-count+1).name;
                    % setting count to 10 to break out from while loop
                    count=N_files_to_search; % The previous file is the first of the survey
                end % if survey1=survey
            else % If you don't find 'CTD.survey', count up.
                count=count+1;
            end %end if CTD.survey exists
        end % end of while loop.
    end
else
    listfile_raw_dir=dir(fullfile(input_struct.raw_dir,'*.modraw'));
    count=0;
    % quick look in the previous files if survey name exist.
    % In case RUN_Auto_process_data was launched few minutes after ./fctd_epsi
    while count<min([N_files_to_search,length(listfile_raw_dir)])
        path2setup1=fullfile(...
            listfile_raw_dir((length(listfile_raw_dir)-count)).folder,...
            listfile_raw_dir((length(listfile_raw_dir)-count)).name);
        fid=fopen(path2setup1,'r');
        fseek(fid,0,1);
        frewind(fid);
        str = fread(fid,'*char')';
        fclose(fid);% If you find 'CTD.survey' in the .mod raw file, get the survey
        % name
        if ~isempty(strfind(str,'CTD.survey'))
            surveyflag1_str      = str(strfind(str,'CTD.survey')+(0:100));
            surveyflag1_str      = surveyflag1_str(1:find(uint8(surveyflag1_str)==10,1,'first'));
            surveyflag1_name     = strsplit(surveyflag1_str,'=');
            surveylag1_name      = surveyflag1_name{2}(1:end-1);
            survey1_name         = surveylag1_name;

            if contains(survey1_name,survey_name)
                input_struct.str_to_match = ...
                    listfile_raw_dir(length(listfile_raw_dir)-count).name;
                count=count+1;
            else
                input_struct.str_to_match = ...
                    listfile_raw_dir(length(listfile_raw_dir)-count+1).name;
                % setting count to 10 to break out from while loop
                count=N_files_to_search; % The previous file is the first of the survey
            end % if survey1=survey
        else
            count=count+1;
        end %end if CTD.survey exists
    end % end of while loop.
end

% If you get to the end of that and haven't found any files that match
% CTD.survey (maybe you're using old .modraw files that don't have them),
% ask for the name of the first file to process
if ~isfield(input_struct,'str_to_match')
    input_struct.str_to_match = input('What is the first file in this deployment? [''EPSIyy_mm_dd_HHMMSS.modraw'']');
end


%% All options have been determined. Now get depth array, create output directory, and get ready to run the processing script
% Get the depth array based on instrument choice
switch instrument
    case {'epsi','EPSI'}
        input_struct.depth_array = epsi_depth_array;
    case {'fctd','FCTD'}
        input_struct.depth_array = fctd_depth_array;
end

% Make the process directory
if ~exist(fullfile(input_struct.process_dir,'raw'),'dir')
    mkdir(input_struct.process_dir);
end

% Set command window color
set_window_color('cyan')

%% Copy Setup file onto the end of Meta_Data_Process file to make a new file with datetime stamped
% Make a new file in the process directory called 'MD' with the
% current datetime
newfile_name = fullfile(input_struct.process_dir,...
    sprintf('MD_%04.0f_%02.0f_%02.0f_%02.0f%02.0f%02.0f.txt',...
    year(now),month(now),day(now),hour(now),minute(now),second(now)));

% Open Meta_Data_Process file for reading and get content
fid1 = fopen(input_struct.Meta_Data_process_file,'r');
content1 = fread(fid1,'*char')';
fclose(fid1);

% Open Setup file for reading and get content
fid2 = fopen(path2setup,'r');
content2 = fread(fid2,'*char')';
fclose(fid2);

% Open the new file for writing amd write content of both files
fidNew = fopen(newfile_name,'w');
fwrite(fidNew,content1,'char');
fwrite(fidNew,content2,'char');
fclose(fidNew);

input_struct.Meta_Data_process_file = newfile_name;


%% Run the processing script on a timer
switch instrument
    case {'epsi','EPSI'}
        epsiAuto_process_data
    case {'fctd','FCTD'}
        fctdAuto_process_data
end

