%%% Script for 2P Project Shock Obs, Laser Exp, and Squeak Playback
%%% Never have animals in line of laser or on footgrid until triggers have
%%% been initialized
%%% Make Sure Laser ALCO2 and Avisoft are connected before Starting
%%% By Steven Voges with help from scripts made my Efe Soyman

%% Get Information And Shuffle RNG
programStart = clock;
rng(sum(programStart));
expinfo = inputdlg({'Observer Animal Number','Demonstrator Animal number',...
    'Scanbox Field 1', 'Scanbox Field 2', 'Scanbox Field 3', 'Filename','Number of Shock Blocks',...
    'Number Shock Conditions','Number of Laser Blocks','Number of Laser Conditions',...
    'Number of total Squeak Trials(Conditions x Blocks','baseline length'});


fileName = expinfo(6);
shockblocks = str2double(expinfo(7));
shockcond = str2double(expinfo(8));
laserblocks = str2double(expinfo(9));
lasercond = str2double(expinfo(10));
squeaktotal = str2double(expinfo(11));
baseline =  str2double(expinfo(12));
%% Initialize DIO24
% Define board and channel numbers
board = 0; % board number
shock = 0;  %shocker 1 (mouse)
shockctrl = 1; %shocker 2(control)
%laser = 2; %laser on
avisoft = 3; %trigger for avisoft


dinit(board); % initialize the board
dasbit(shock,1); % set grid 1 to 1 (reversed)
dasbit(shockctrl,1);
%dasbit(laser,0)
dasbit(avisoft,1);

%% Randomize ITI and Condition Order

shocktotal = shockblocks * shockcond;
lasertotal = laserblocks * lasercond;
shocktriallist = zeros(shocktotal,7);
lasertriallist = zeros(lasertotal,7);
squeaktriallist = zeros(squeaktotal,6);
shocktriallist(:,1) = trialmaker(shockblocks,shockcond);
lasertriallist(:,1) = trialmaker(laserblocks,lasercond);
itishock = 30+round(rand((shocktotal),1)*30);
itisqueak =  30+round(rand((squeaktotal),1)*30);
itilaser = 30+round(rand((lasertotal),1)*30);
squeakbase = 30;
laserbase=30;
%% Set-Up Motor Controll
%motors_initialize
result=move_absolute_all(0,7,6.5);

%% Set Y Coordinates

xcor = 7;
ycor = [8 6.5 5 8 6.5 5 8 6.5 5 8 6.5 5 8 6.5 5 8 6.5 5 8 6.5 8 6.5 5 8 6.5 5 8 6.5 5 8 6.5 5 8 6.5 5 8 6.5 5  ]; 

%% Create Video Object
vid = videoinput('gentl', 1, 'Mono8');
src = getselectedsource(vid);

vid.FramesPerTrigger = 1;

src.ExposureTimeAbs = 10000;

preview(vid);

src.TriggerMode = 'On';

vid.LoggingMode = 'disk';

diskLogger = VideoWriter(char(fileName), 'Grayscale AVI');

vid.DiskLogger = diskLogger;

diskLogger.FrameRate = 15.5;

% TriggerRepeat is zero based and is always one
% less than the number of triggers.
vid.TriggerRepeat = Inf;

triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');%% Video will Cut because it needs a hardware trigger




%% Set Up ALCO2 Laser Control

arduino=serial('COM3','BaudRate',115200);
arduino.OutputBufferSize = 1024;
arduino.Terminator = 'CR/LF';
%Set COM Settings
fopen(arduino);%Open COM
pause(1)
fprintf(arduino, 'ALCO2!'); %Open ALCO@
pause(1)
fprintf(arduino, 'ALCO2:Start!'); %Tell ALCO that you ready
pause(1)
fprintf(arduino,'ALCO2:Set_RemoteOFF!');
cla
figure
text(.5,.5,'Turn on Diode on Front Panel of Laser and Press Any Key to Continue','HorizontalAlignment','Center','VerticalAlignment','Middle')%% Set visible laser diode on with front panel
axis off

% Halt the program until keyboard press
waitforbuttonpress;
cla
fprintf(arduino,'ALCO2:Set_RemoteON!');% Switches off front panel and enables control of laser via ALCO2
pause(0.5)
fprintf(arduino,'ALCO2:Set_PWMDuration_uS,30000!');% Set Laser Duration to Lowest Setting
pause(0.5)
fprintf(arduino,'ALCO2:Set_PWMDutyCycle_Percent,30!');% Set Laser Power to Lowest Setting



%% Insert Observer and start video recording

start(vid);
% Display program status
cla
figure
text(.5,.5,'Check that grid is off and insert demonstrator into set up, press any key to continue','HorizontalAlignment','Center','VerticalAlignment','Middle')
axis off

% Halt the program until keyboard press
waitforbuttonpress;


%% Setup and start scanbox





sbudp = udp('192.87.11.27', 'RemotePort', 7000);% Parameters for connection to scanbox

fopen(sbudp);% Connect to scanbox

%fprintf(sbudp,sprintf('A%s',expinfo{3})); % set animal id field
%fprintf(sbudp,sprintf('U%s',expinfo{4})); % set field number to 010
%fprintf(sbudp,sprintf('E%s',expinfo{5})); % set experiment number to 002
figure
cla
text(.5,.5,'Is Scanbox Focused and Ready to Go? Press Any Key To Start Scanbox Recording','HorizontalAlignment','Center','VerticalAlignment','Middle')
axis off 



%Halt the program until keyboard press
waitforbuttonpress;

fprintf(sbudp,'G');    % Go!  Start sampling
%% Run shock trials
cla
text(.5,.5,sprintf('Baseline Before Shock: %d seconds',baseline),'HorizontalAlignment','Center','VerticalAlignment','Middle')
axis off
pause(20)
fprintf(sbudp,'L0');
pause(baseline - 40);
fprintf(sbudp,'L1');
pause(20)
for i = 1:shocktotal
    if shocktriallist(i)== 1 %shock
        
        cla
        text(.5,.5, sprintf('Trial #%d: Shock! iti: %d seconds',i,itishock(i)),'HorizontalAlignment','Center','VerticalAlignment','Middle')
        axis off
        dasbit(shock,0)
        shocktriallist(i,2:7) = clock;
        pause(2)
        dasbit(shock,1)
        pause(10);
        fprintf(sbudp,'L0');
        pause(itishock(i)-20);
        fprintf(sbudp,'L1');
        pause(10)
    elseif shocktriallist(i)== 2 % control
        
        cla
        text(.5,.5,sprintf('Trial #%d: Control! iti: %d seconds ',i,itishock(i)),'HorizontalAlignment','Center','VerticalAlignment','Middle')
        axis off
        
        dasbit(shockctrl,0)
        shocktriallist(i,2:7) = clock;
        pause(2)
        dasbit(shockctrl,1)
        
        
        pause(10);
        fprintf(sbudp,'L0');
        pause(itishock(i)-20);
        fprintf(sbudp,'L1');
        pause(10)
        
    end
end
pause(10)
fprintf(sbudp,'L0');
cla
text(.5,.5,'Shock Trials Over Prepare for Laser Trial and press any Key when ready','HorizontalAlignment','Center','VerticalAlignment','Middle')
axis off



%Halt the program until keyboard press
waitforbuttonpress;

%% Run Laser trials
counter = 1;
cla

text(.5,.5,sprintf('Baseline Before Laser: %d seconds',laserbase),'HorizontalAlignment','Center','VerticalAlignment','Middle')
axis off
fprintf(sbudp,'L1');
pause(10) % Pause for baseline
fprintf(sbudp,'L0');
pause(10);
fprintf(sbudp,'L1');
pause(10)

for i = 1:lasertotal
    if lasertriallist(i,1)==1 %ctrl laser
        
        cla
        text(.5,.5, sprintf('Trial #%d: Laser Ctrl! iti: %d seconds',i,itilaser(i)),'HorizontalAlignment','Center','VerticalAlignment','Middle')
        axis off
        fprintf(arduino,'ALCO2:Set_PWMDuration_uS,30000!');
        
        pause(0.1)
        fprintf(arduino,'ALCO2:Set_PWMDutyCycle_Percent,30!');
        
        pause(0.1)
        
        fprintf(arduino,'ALCO2:PulseOut!');
        lasertriallist(i,2:7) = clock;
        pause(10);
        fprintf(sbudp,'L0');
        pause(itilaser(i)-20);
        fprintf(sbudp,'L1');
        pause(10)
        
    elseif lasertriallist(i,1)==2
        
        cla
        text(.5,.5, sprintf('Trial #%d: Laser Low! iti: %d seconds',i,itilaser(i)),'HorizontalAlignment','Center','VerticalAlignment','Middle')
        axis off
        fprintf(arduino,'ALCO2:Set_PWMDuration_uS,50000!');
        pause(0.1)
        fprintf(arduino,'ALCO2:Set_PWMDutyCycle_Percent,55!');
        disp('mid')
        pause(0.1)
        cla
        text(.5,.5,'Is Mouse in Line of Laser press any key for pulse out','HorizontalAlignment','Center','VerticalAlignment','Middle')
        axis off
        waitforbuttonpress;
        
        
        %Halt the program until keyboard press
        
        fprintf(arduino,'ALCO2:PulseOut!');
        lasertriallist(i,2:7) = clock;
        
        pause(8);
        move_absolute_all(0,xcor,ycor(counter))
        pause(2)
        fprintf(sbudp,'L0');
        pause(itilaser(i)-20);
        fprintf(sbudp,'L1');
        pause(10)
        counter = counter+1;
    elseif lasertriallist(i,1)==3
        
        cla
        text(.5,.5, sprintf('Trial #%d: Laser High! iti: %d seconds',i,itilaser(i)),'HorizontalAlignment','Center','VerticalAlignment','Middle')
        axis off
        fprintf(arduino,'ALCO2:Set_PWMDuration_uS,60000!');
        pause(0.1)
        fprintf(arduino,'ALCO2:Set_PWMDutyCycle_Percent,60!');
        disp('high')
        pause(0.1)
        cla
        text(.5,.5,'Is Mouse in Line of Laser press any key for pulse out','HorizontalAlignment','Center','VerticalAlignment','Middle')
        axis off
        waitforbuttonpress;
        fprintf(arduino,'ALCO2:PulseOut!');
        lasertriallist(i,2:7) = clock;
        
        pause(8);
        move_absolute_all(0,xcor,ycor(counter))
        pause(2)
        
        fprintf(sbudp,'L0');
        pause(itilaser(i)-20);
        fprintf(sbudp,'L1');
        pause(10)
        counter = counter+1;

    end
end
fprintf(sbudp,'L0');
cla
text(.5,.5,'Laser Trials Over Prepare for Squeak Trial and press any Key when ready','HorizontalAlignment','Center','VerticalAlignment','Middle')
axis off



%Halt the program until keyboard press
waitforbuttonpress;


%% Run Squeak Trials - Playlist must be Randomized before hand

cla
text(.5,.5,sprintf('Baseline Before Squeak: %d seconds',squeakbase),'HorizontalAlignment','Center','VerticalAlignment','Middle')
axis off
fprintf(sbudp,'L1');
pause(10) % Pause for baseline
fprintf(sbudp,'L0');
pause(10);
fprintf(sbudp,'L1');
pause(10) % 30 second baseline
squeaktriallist = zeros(squeaktotal,6);
for i = 1:squeaktotal
    
    cla
    text(.5,.5,sprintf('Trial #%d: Squeak ! iti: %d seconds ',i,itisqueak(i)),'HorizontalAlignment','Center','VerticalAlignment','Middle')
    axis off
    dasbit(avisoft,0)
    squeaktriallist(i,1:6) = clock;
    pause(1)
    dasbit(avisoft,1)
    pause(10);
    fprintf(sbudp,'L0');
    pause(itisqueak(i)-20);
    fprintf(sbudp,'L1');
    pause(10)
end


%% Stop Scanbox and Camera and Save Log Files
fprintf(sbudp,'S');    % Stop sampling
stoppreview(vid);
stop(vid);
fclose(sbudp);
save(char(fileName), 'expinfo','itishock','shocktriallist','itilaser','lasertriallist','itisqueak','squeaktriallist')