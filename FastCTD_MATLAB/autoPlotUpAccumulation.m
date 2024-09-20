clear;
clc;
cmdWinDoc=com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
listeners = cmdWinDoc.getDocumentListeners;
%find text area part
jTextArea=listeners(5); %or listeners(3) or listeners (4) depending on matlab
%set colour of command window
jTextArea.setBackground(java.awt.Color(0.5,0.2,0.2));


f = figure(1);
f.Units = 'normalized';
set(f,'position',[0 0,0.6,0.6]);
clf;

FCTD_rot_timer = timer;
s = false;

FCTD_rot_timer.StartFcn = 'disp(''Plot FCTD rotation begins now!'');';
FCTD_rot_timer.TimerFcn = [...
    'if s, '...
    'stop(FCTD_rot_timer); '...
    'delete(FCTD_rot_timer); '...
    'else, '...
    'try,'...
    'PlotUpAccumulation;'...
    'catch err,'...
    'display_error_stack(err); '...
    'end;'...
    'disp([''Done at '' datestr(now)]);'...
    'end;'];
FCTD_rot_timer.Period = 60*5;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
FCTD_rot_timer.BusyMode = 'drop';
FCTD_rot_timer.Name = 'FCTD_rot_timer';
FCTD_rot_timer.Tag = 'FCTD_rot_timer';
FCTD_rot_timer.StopFcn = 'clear(''rawDir'',''rawDirAway''); disp([datestr(now) '': Stopped FCTD_rot_timer'']);';
FCTD_rot_timer.ExecutionMode = 'fixedSpacing';
% FCTD_rot_timer.ExecutionMode = 'singleShot';
FCTD_rot_timer.TasksToExecute = Inf;
FCTD_rot_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(FCTD_rot_timer);