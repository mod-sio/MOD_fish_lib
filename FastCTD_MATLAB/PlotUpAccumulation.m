% ALB matDataDir is the FCTDmat folder of the current deployment/section
% Be carefull to update file_index_ofset and rot_count_offset
%
% rot_count_offset is the last rot_count_offset obtained at the previous
% deployment
%
% file_index_ofset is the 
%
%%

matDataDir = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/MAT_full_cruise/'; %Directory where FCTD*.mat are stored
rotDataDir = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/ROT_full_cruise/'; %Directory where rotation data for each FCTD*.mat file will be stored
%matFiles = dir([matDataDir 'FCTD*.mat']);
matFiles = dir([matDataDir 'EPSI*.mat']);
rotFiles = dir([rotDataDir 'EPSI*.mat']);


% if ~exist(rotDataDir,'dir')
%     mkdir(rotDataDir)
% end


% -----------------------------------------------------------------------------


% matDataFilesLoaded{end+1:numel(matFiles)} = cell(1,numel(matFiles)-numel(matDataFilesLoaded));
for i = 1+numel(rotFiles):numel(matFiles)
    %find out what files are new
    % ind = find(strcmpi(matFiles(i).name,{matDataFilesLoaded.filename}));
    % if ~isempty(ind)
    %     %search for time differences
    %     if(matFiles(i).datenum <= matDataFilesLoaded(ind).time)
    %         continue;
    %     end
    % end

    %load FCTD data
    %disp(matFiles(i).name);
    load([matDataDir '/' matFiles(i).name]);

    % %if the index file is empty, create new index
    % if isempty(matDataFilesLoaded)
    %     matDataFilesLoaded = struct('filename',matFiles(i).name,'time',matFiles(i).datenum);
    % else
    %     if ~isempty(ind)
    %         matDataFilesLoaded(ind) = struct('filename',matFiles(i).name,'time',matFiles(i).datenum);
    %     else
    %         matDataFilesLoaded(end+1) = struct('filename',matFiles(i).name,'time',matFiles(i).datenum);
    %     end
    % end

    % convert mat data file to rotation accumulation file
    if exist('vnav','var') && exist('ctd','var') && ...
       isstruct(vnav) && ~isempty(vnav.time_s)
    
        diff_not_neg = [0;diff(vnav.dnum)]>0;
        keep = ~isnan(vnav.dnum) & ~isinf(vnav.dnum) & diff_not_neg;
        rot.compass      = vnav.compass(keep,:);
        rot.gyro         = vnav.gyro(keep,:);
        rot.acceleration = vnav.acceleration(keep,:)./9.81;
        rot.time_s       = vnav.time_s(keep);
        rot.dnum         = vnav.dnum(keep);
        if ~isempty(ctd)
            rot.pressure     = interp1(ctd.dnum,ctd.P,rot.dnum);
        else
            rot.pressure     = rot.dnum.*nan;
        end

        time = rot.dnum;
        pts = numel(time);
%         multiplier = -1;
        multiplier = 1;

        %ALB why is this line here?
        % Phi and other are redefined later, so no need for this line.
        [Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis([0, 0,1]);

        Rot_Mat = @(p,t,s)[ cos(t)*cos(s), -cos(p)*sin(s) + sin(p)*sin(t)*cos(s),  sin(p)*sin(s) + cos(p)*sin(t)*cos(s);
            cos(t)*sin(s),  cos(p)*cos(s) + sin(p)*sin(t)*sin(s), -sin(p)*cos(s) + cos(p)*sin(t)*sin(s);
            -sin(t),         sin(p)*cos(t),                         cos(p)*cos(t)];

        acce = rot.acceleration; %(Rot_mat*(rot.acceleration'))';
        acce(:,3) = multiplier*acce(:,3);

        acc_xy_length = sqrt(sum(acce(:,1:2).^2,2));
        acc_length = sqrt(sum(acce.^2,2));

        comp = rot.compass;
        gyro = rot.gyro;
        acce = rot.acceleration;

        comp(:,3) = multiplier*comp(:,3);
        gyro(:,3) = multiplier*gyro(:,3);
        acce(:,3) = multiplier*acce(:,3);
        comp = comp';
        gyro = gyro';
        acce = acce';

        new_comp = NaN(size(comp))';
        new_gyro = NaN(size(comp))';
        new_acce = NaN(size(comp))';
        tic
        for k = 1:pts
            [Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis(acce(:,k)');
            new_comp(k,:) = Rot_mat*(comp(:,k));
            new_gyro(k,:) = Rot_mat*(gyro(:,k));
            new_acce(k,:) = Rot_mat*(acce(:,k));

        end
        toc
        pressure     = rot.pressure;
        acceleration = rot.acceleration;
        org_gyro     = gyro';
        gyro         = new_gyro;
        compass      = new_comp; %medfilt1(new_comp,5,[],1);
%         compass(end-3:end,:) = new_comp(end-3:end,:);
        COMP_mag = repmat(sqrt(sum(compass(:,1:2).*compass(:,1:2),2)),[1 3]);
        compass = compass./COMP_mag;
        compass = compass(:,1)+1i*compass(:,2);

        dt = diff(time)*24*3600;
        dt = [dt; nanmedian(dt)];

        tot_rot_gyro = cumsum(gyro(:,1).*dt);
        tot_rot_gyro(:,2) = cumsum(gyro(:,2).*dt);
        tot_rot_gyro(:,3) = cumsum(gyro(:,3).*dt);
        tot_rot_acc = phase(compass);
        compass = new_comp;

        save([rotDataDir '/' matFiles(i).name],'compass','tot_rot_gyro','tot_rot_acc','time','org_gyro','gyro','acceleration','pressure');
    end
    clear FCTD;
end
%disp('Done');

% file_index_ofset = 1; %index of file in 'matFiles' to start plotting
% rot_count_offset = 0; %number of total twists that were on the cable at the end of file_index_offset-1

%%
rotFiles = dir([rotDataDir 'EPSI*.mat']);
% rotFiles = dir([rotDataDir '/deployment_1/FCTD*.mat']);

if exist([rotDataDir 'Latest_rot_acc_count.mat'],'file')
load([rotDataDir 'Latest_rot_acc_count.mat'], ...
     'tot_rot_gyro', ...
     'tot_rot_acc',...
     'tot_time',...
     'tot_pressure',...
     'last_file_idx');
    file_index_ofset=last_file_idx;
    nonan_tot_rot_acc=tot_rot_acc(~isnan(tot_rot_acc));
    rot_count_offset=-nonan_tot_rot_acc(end)/pi/2;

else
    tot_rot_gyro     = [];
    tot_rot_acc      = [];
    file_index_ofset = 1;
    tot_time         = [];
    tot_pressure     = [];    
    rot_count_offset = 0;
end

% compass      = [];
% gyro         = [];
% org_gyro     = [];
% pressure     = [];
% time         = [];
% acceleration = [];
tic
for i = file_index_ofset:numel(rotFiles)
    %disp(rotFiles(i).name);
    rot = load([rotDataDir '/' rotFiles(i).name]);
    % Make sure there are no nans in the data so you can sum
    % rot.tot_rot_gyro(isnan(rot.tot_rot_gyro)) = 0;
    % rot.tot_rot_acc(isnan(rot.tot_rot_acc)) = 0;
    % last_rot_acc  = find(~isnan(rot_total.tot_rot_acc),1,'last');
    if ~isempty(tot_time)
        start_idx=find(rot.time>=tot_time(end),1,'first');
        if isempty(start_idx)
            start_idx=1;
        end
    else
        start_idx=1;
    end
    last_rot_acc  = find(~isnan(tot_rot_acc),1,'last');
    if isempty(last_rot_acc)
        last_rot_acc  = 1;
    end
    last_rot_gyro = find(~isnan(tot_rot_gyro),1,'last');
    if isempty(last_rot_gyro)
        last_rot_gyro  = 1;
    end

    % rot = load([rotDataDir '/deployment_1/' rotFiles(i).name]);
    % compass = [compass; rot.compass(start_idx:end,:)];
    % gyro = [gyro; rot.gyro(start_idx:end,:)];
    % org_gyro = [org_gyro; rot.org_gyro(start_idx:end,:)];
    % acceleration = [acceleration; rot.acceleration(start_idx:end,:)];
    % tot_rot_gyro = [tot_rot_gyro; tot_rot_gyro(end,:)+rot.tot_rot_gyro];
    % tot_rot_acc = [tot_rot_acc; tot_rot_acc(end,:)+rot.tot_rot_acc];
    if ~isempty(tot_rot_gyro)
        tot_rot_gyro = [tot_rot_gyro(1:end-1); tot_rot_gyro(last_rot_gyro)+rot.tot_rot_gyro(start_idx:end,3)];
        tot_rot_acc  = [tot_rot_acc(1:end-1); tot_rot_acc(last_rot_acc)+rot.tot_rot_acc(start_idx:end)];
    else
        tot_rot_gyro = rot.tot_rot_gyro(start_idx:end,3);
        tot_rot_acc  = rot.tot_rot_acc(start_idx:end);

    end
    tot_pressure = [tot_pressure(1:end-1); rot.pressure(start_idx:end)];
    tot_time = [tot_time(1:end-1); rot.time(start_idx:end)];    
    last_file_idx=i;
    clear rot;
end
% tot_rot_gyro = tot_rot_gyro(2:end,:);
% tot_rot_acc = tot_rot_acc(2:end,:);
% tot_time=[tot_time(1:end-1) ; time];
% tot_pressure=[tot_pressure(1:end-1) ; pressure];


save([rotDataDir 'Latest_rot_acc_count.mat'], ...
     'tot_rot_gyro', ...
     'tot_rot_acc',...
     'tot_time',...
     'tot_pressure',...
     'last_file_idx');

toc
disp('Done');
fprintf('Last numel(rotFiles)=%i\n',numel(rotFiles))


%%
% t_offset = datenum(2019,07,00);
%SN_figure(2,'w',1500,'h',900);
if isempty(time)
    disp('no new file... Patience my friend.')
    last_value =rot_count_offset;
else
dt = diff(tot_time)*24*3600;
dt = [dt(1); dt];
dp = diff(tot_pressure);
dp = [dp(1); dp];
fall_rate = dp(:)./dt(:);
% figure(2);
% set(2,'position',[0 0 1500 900])
h = [0,0,0,0];
clf;
% h(1) = plot(datetime(tot_time,'ConvertFrom','datenum'), -tot_rot_acc/pi/2,'linewidth',2,'color','r');
% hold on;
% if sum(fall_rate<0.001)>10
%     h(2) = plot(datetime(tot_time(fall_rate<0.01),'ConvertFrom','datenum'), -tot_rot_acc(fall_rate<0.01)/pi/2,'.','linewidth',2,'color','k');
% end
last_value = -tot_rot_acc/pi/2;
last_value=last_value(~isnan(last_value));

disp(['Most recent turn count: ' datestr(now),'   ' num2str(round(last_value(end)))])

h(3) = plot(datetime(tot_time,'ConvertFrom','datenum'), tot_rot_gyro/pi/2,'linewidth',2,'color','b');
if sum(fall_rate<0.001)>10
    h(4) = plot(datetime(tot_time(fall_rate<0.001),'ConvertFrom','datenum'), tot_rot_gyro(fall_rate<0.001)/pi/2,'.','linewidth',2,'color','c');
end
% set(gca,'XTick',datetime(time(1):12/24:time(end),'ConvertFrom','datenum'))
set(gca,'XTickLabelRotation',45)
hold off;
grid on;
xlabel('time');
ylabel('Number of rotations');
%title('FCTD: Rotation count on current line');

title(['MODfish: Rotation count = ' num2str(round(last_value(end))) '  _ _ _ ']);

str_legend={'Rot by acc [dn]','Rot by acc [up]','Rot by gyro [dn]','Rot by gyro [up]'};
hl = legend(h(h>0),str_legend(h>0),'Location','NorthWest');
set(hl,'Fontsize',20);
set(gca,'Fontsize',20);
ax(1) = gca;

%SN_setTextInterpreter('latex');
% SN_printfig([rotDataDir '/currentRotFig.pdf']);
%
% figure(3);
% set(2,'position',[0 0 1500 900])
% clf;
% plot(datetime(time,'ConvertFrom','datenum'), acceleration,'linewidth',2);
% hold on;
% plot(datetime(time(fall_rate<0.01),'ConvertFrom','datenum'), acceleration(fall_rate<0.01,:),'.','linewidth',2,'color','k');
%
% hold off;
% grid on;
% xlabel('time');
% ylabel('acceleration');
% % title('FCTD: Rotation count on current line');
%
%
% ax(2) = gca;
%
% linkaxes(ax,'x')
end

fprintf("Last rotation count %i\r\n",round(last_value(end)))
