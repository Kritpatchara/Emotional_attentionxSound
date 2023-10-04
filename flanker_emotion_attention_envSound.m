 
%% Real
%% Flanker task: Examining emotional attention effect across different env. sounds
%% Condition information
% Face cue:
% 1. Emotion x 3 (HA, AN, NE)
% 2. Location x 12
% 3. Timing x 2 (125/ 500 ms)

% Simulus (bowtie, diamond):
% 4. Congruency x 2 (con, incon)
% Note: congruency means 'aaaaa' and incongruency 'aabaa'

% Total number of condition: 3x12x2x2 = 144

% For each subject: 2 sound conditions (positive and negative)

% Total trials per subject: 144 x 2 = 288

% Number of face exemplars = 24 (each exemplar repeated 36 times; exemplars
% are not counterbalanced with other conditions)
% Note: exemplar selection ensures that each photo has >70% accuracy score,
% and that mean acc scores do not vary across emotions (HA, AN, NE/CA), and
% that emo level ratings do not differ across HA vs. AN categories
% actorID selected:  1  7  8  9 10 11 13 14 17 18 19 20 22 23 24 25 26 32 34 35 36 37 42 45
% (11F, 13M)

% Timing information
% 1. face cue = 0.075s
% 2. cue-stimulus-interval = 0.175s
% 3. stimulus = 1.55s
% 4. search duration = 2.05s (maximum)
% 5. feedback = 0.3s
% 6. inter-trial-interval = 1s

% estimated total time per trial = 0.75+0.175+1.55+2.05+0.3+1 = 5.15

%
%% close & clear everything
clear all;
close all;

%% get subject info & audio file info

% Show dialog box
prompt = {'SubjectID (e.g., 1234)','Audio file (pos=1, neg=2):'};
defAns = {'9999','1'};
box = inputdlg(prompt,'Enter Subject Info and Audio File...', 1, defAns);

% Record subnum and audiotype
p.subNum = box{1};
p.audioType = box{2};

%% Set a few other things
% Set EEG
p.EEG = 0;

% Set response deadline
p.rspdeadline = 3;

% Set current working directory
p.root = pwd;

% Set the run number to 1
p.runNum = 1;

% Connect to serial port if EEG is used 
if p.EEG == 1
    SerialPort = BioSemiSerialPort();
end

% Initiate a few other variables
p.s = round(sum(100*clock));
rand('state', p.s);


%% Get demographics info (handedness, gender, age)

% Show dialog box to collect gender, hand, age info
prompt = {'Gender (m/f)', 'Hand (l/r)', 'Age'};
defAns = {'x', 'x' , '0'};
box = inputdlg(prompt,'Enter Subject Information...', 1, defAns);

% Store answers in variables
p.gender = box{1};
p.hand = box{2};
p.age = str2double(box{3});


%% Create files and folders for storing data from the experiment

ListenChar(2);

% If a folder called 'data_envSound' does not already exist in the working directory (p.root), create one 
if ~exist([p.root, '/data_envSound/'], 'dir')
    mkdir([p.root, '/data_envSound/']);
end

% Create a folder for the particular subject
mkdir([p.root, '/data_envSound/s', num2str(p.subNum)]);

% Within the created folder, create fil names for saving data
fName = [p.root, '/data_envSound/s', num2str(p.subNum), '/flanker' '_sbj',  num2str(p.subNum), '_audType_', p.audioType, '.mat'];

% If fName already exists, abort the program
if exist(fName,'file')
    Screen('CloseAll');
    msgbox('File name already exists, please specify another', 'modal');
    ListenChar(0);
    return
end

% Create a folder for r, if not already exist
if ~exist([p.root, '/data_envSound_csv/'], 'dir')
    mkdir([p.root, '/data_envSound_csv/']);
end
%% set up screen - load the CLUT to get the correct gamma value, open a screen that is filled with grey

AssertOpenGL;

% Open onscreen window:
Screen('Preference', 'SkipSyncTests', 1);

% Determine which screen to use (if error happened related to 'Screen', may
% need to adjust here)
screen = max(Screen('Screens'));

% Get window characteristics
[win, scr_rect] = PsychImaging('OpenWindow', screen, 128, round([-1, -1, 1921, 1081])); % hard code debugging
[winWidth, winHeight] = Screen('WindowSize', win);

% Set priority levels
priorityLevel = MaxPriority(win);
Priority(priorityLevel);
flipT = Screen('GetFlipInterval', win);

% Specify a few variables related to display colors
black = 0;
white = 255;
red = [255 0 0];
blue = [0 0 255];
rd = 255;
bl = 255;
background = 128; % background color

% Specify variables related to screen center
p.xcenter = winWidth/2; 
p.ycenter = winHeight/2;

% Fill screen with background color:
Screen('FillRect', win, background);

% Initialize display and sync to timestamp:
[VBLTimestamp] = Screen(win, 'Flip');

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% load images from 'stim' folder

% Get the folder directories where all the stim images are contained
tgfolder = [p.root, '/stim/'];
facefolder = [p.root, '/faces/'];

% Load bow image (also adjusting all white background to grey (RGB=128))
bow = imread([tgfolder 'bowtie.png']);
bow(bow == 255) = 128;

% Load diamond image (also adjusting all white background to grey (RGB=128))
diamond = imread([tgfolder 'diamond.png']);
diamond(diamond == 255) = 128;

% Create textures for bow and diamond images
tgtext{2,1} = Screen('MakeTexture', win, bow);
tgtext{1,1} = Screen('MakeTexture', win, diamond);

% Load & create image texture for feedback of correct responses
right = imread([tgfolder 'right.png']);
right(right == 0) = 128;
right(right <= 30) = 0;
right = Screen('MakeTexture', win, right);

% Load & create image texture for feedback of wrong responses
wrong = imread([tgfolder 'wrong.png']);
wrong(wrong == 0) = 128;
wrong(wrong <= 30) = 0;
wrong = Screen('MakeTexture', win, wrong);

% Load & create image texture for feedback of slow responses
slow = imread([tgfolder 'slow.png']);
slow(slow == 0) = 128;
slow(slow <= 30) = 0;
slow = Screen('MakeTexture', win, slow);

% Load & create image texture for instruction screen
init = imread([tgfolder 'start.png']);
init = Screen('MakeTexture', win, init);

%% load sound

% Load audio files (file to load depending on p.audioType)
% ONLY FOR NOW: Use cafe ambience just as an example for positive sound 
if p.audioType == '1'
    [x, envSound] = audioread(strcat(p.root, '/sound/nature_sound.mp3'));
elseif p.audioType == '2'
    [x, envSound] = audioread(strcat(p.root, '/sound/traffic_sound.mp3'));
end


%% set screen locations for different images

% Set radial distance for stimuli (bowtie & diamond). Unit in pixels
p.Rec_stim = 390;

% Set radial distance for cues (face images)
p.Rec_cue = 290;

% Set radial distance for central fixation
p.Rfix = 5;

% Set radial size for stimuli (bowtie & diamond)
p.Rstim = 70;

% Set height and width for face images
p.Hstim = p.Rstim*1.8;
p.Wstim = p.Hstim * (301/451); % This is the calculated W/H ratio of the face image

% This is for feedback images
p.Rcue = 75;
p.Rring = 140;

% Translate the above information in to coordinate distance (x,y)
deg_loc = 0:2*pi/12:2*pi-pi/12;
p.LocX = p.xcenter+p.Rec_stim*(cos(deg_loc)); % stimulus
p.LocY = p.ycenter+p.Rec_stim*(sin(deg_loc));
p.LocXc = p.xcenter+p.Rec_cue*(cos(deg_loc)); % cue
p.LocYc = p.ycenter+p.Rec_cue*(sin(deg_loc));

%% set response keys

% Get all KbName
KbName('UnifyKeyNames');

% Set allowed key responses
respCode = [KbName('j') KbName('k')];

% Set abort key
abortKey = KbName('ESCAPE');

FlushEvents;

%% set up experimental conditions

% Specify the number of total trials in the experiment
p.numTrials = 144;

% Create matrix for cue location (1,2,3,...,12)
p.cueloc = repmat([...
    ones(1, 1)*1 ones(1, 1)*2 ones(1, 1)*3 ones(1, 1)*4 ...
    ones(1, 1)*5 ones(1, 1)*6 ones(1, 1)*7 ones(1, 1)*8 ...
    ones(1, 1)*9 ones(1, 1)*10 ones(1, 1)*11 ones(1, 1)*12], [1, 12]);

% Create matrix for facial emotion (HA=1, AN=2, NE=3)
p.emo_HA_AN_NE = repmat([ones(1, 12) ones(1, 12)*2 ones(1, 12)*3], [1, 4]);

% Create matrix for csi (125=1, 500=2)
p.csi_125_500 = repmat([ones(1, 36)*1 ones(1, 36)*2], [1, 2]);

% Create matrix for stimulus congruency (con=1, incon=2)
p.coninc = repmat([ones(1, 72) ones(1, 72)*2], [1, 1]);

%     % sanity check (plotting to see if everything looks ok..)
%     plot(p.cueloc); hold on; plot(p.emo_HA_AN_NE); hold on; plot(p.csi_125_500); hold on; plot(p.coninc);

% Create matrix for stimulus shapes (diamond, hourglass)
p.shape_diamondglass = [ones(1,72)*1 ones(1,72)*2];

% Randomize order of each variable by shuffled trial order
trialorder = Shuffle(1:p.numTrials);
p.cueloc(trialorder) = p.cueloc;
p.emo_HA_AN_NE(trialorder) = p.emo_HA_AN_NE;
p.csi_125_500(trialorder) = p.csi_125_500;
p.coninc(trialorder) = p.coninc;

% Randomize separately for this variable (Hope this is right?....)
trialorder = Shuffle(1:p.numTrials);
p.shape_diamondglass(trialorder) = p.shape_diamondglass;

% Initiate variable for face exemplars
p.exshuffle = Shuffle(repmat(1:24, [1,6]));

% Initiate variables for collecting demographic data
p.hand = repmat(p.hand,p.numTrials,1);
p.gender = repmat(p.gender,p.numTrials,1);
p.age = repmat(p.age,p.numTrials,1);


%% Create 'p.stimshape' variable containing shape and location information for stimuli 

% Initiate p.stimshape variable
p.stimshape = randi(2, [12, p.numTrials]);

% Create disloc variable 
disloc = [12 1:11; 2:12 1; ... % near flankers
    11:12 1:10 ; 3:12 1:2]; % add this line for far flankers

for t = 1:p.numTrials

    p.stimshape(p.cueloc(t), t) = p.shape_diamondglass(t);

    % Adjust stimulus location according to whether p.coninc is con or incon 
    if p.coninc(t) == 1
        
        % If congruent, then adjust the shapes of the flanking stimuli so that they are the same as target 
        p.stimshape(disloc(:, p.cueloc(t)), t) = p.shape_diamondglass(t);
   
    elseif p.coninc(t) == 2

        % If incongruent, adjust the shapes of flanking stimuli so that they are different from target
        if p.shape_diamondglass(t) == 1
            p.stimshape(disloc(:, p.cueloc(t)), t) = 2;
        elseif p.shape_diamondglass(t) == 2
            p.stimshape(disloc(:, p.cueloc(t)), t) = 1;
        end
    end

end

%% read cue images (face images) and create face texture

% Create labels for emotions and face exemplars (i.e., actorID)
emo_label = {'HA','AN','NE'};
actorIdSelected = {'01F', '07F', '08F', '09F', '10F', '11F', '13F', '14F',...
    '17F', '18F', '19F', '20M', '22M', '23M', '24M', '25M', '26M', '32M',...
    '34M', '35M', '36M', '37M', '42M', '45M'};

% Shuffle order of face exemplar
exshuffle = Shuffle(actorIdSelected); 

% For each face exemplar and for each emotion:
for emo = 1:3
    for ex = 1:24

        % Read all face images
        face = imread([facefolder actorIdSelected{ex} '_' emo_label{emo}(1:2) '_C.png']);

        % Create face texture
        facetext{emo, ex} = Screen('MakeTexture', win, face);

    end
end


%% set timing parameters

% Duration for face image (i.e., cue)
p.cuedur = 0.075;
p.cuedur_real = nan(1, p.numTrials);

% Duration for feedback
p.feedback = 0.3;
p.feeddur_real = nan(1, p.numTrials);

% Duration for stimuli (bowtie, diamond)
p.stimdur = 1.55;

% Duration for response
p.searchdur = p.stimdur+.5;
p.searchdur_real = nan(1, p.numTrials);

% Duration for inter-stimulus-interval (isi)
p.isi = Shuffle(linspace(0.4, 0.6, p.numTrials)); 
p.isi_real =  nan(1, p.numTrials);

% Duration for inter-trial-interval (iti)
p.iti = Shuffle(linspace(0.6, 1, p.numTrials)); 
p.iti_real =  nan(1, p.numTrials);

% Duration for subjects' respsones (accuracy, RT, response type, slow or
% not)
p.hit = zeros(1, p.numTrials);
p.hit_notslow = zeros(1, p.numTrials);
p.RT = nan(1, p.numTrials);
p.resp = nan(1, p.numTrials);
p.slow = ones(1, p.numTrials);


%% STARTING EXPERIMENT

% Create textures for instruction screen (init)
Screen('DrawTexture',win, init, [], [], []);
intime = num2str(p.rspdeadline);

% Create texture for central fixation)
DrawFormattedText(win, [intime], p.xcenter+380, p.ycenter-115, black);
Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15); 

% Flip: Present instruction on screen (also wait till any key is pressed before moving on)
Screen(win, 'Flip');
KbWait;

% Present fixation for 2s break
Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15); % fixation
[VBLTimestamp] = Screen(win, 'Flip');
WaitSecs(2);

% Play sound
loadsound1 = x;
loadsound2 = envSound;
sound(loadsound1, loadsound2);

% Begin loop (N loop = p.numTrials)
%% trial loop
for t = 1:p.numTrials

    %% Display face cue

    % Specify target location
    loc =  p.cueloc(t);

    % Create face texture (at specified location)
    Screen('DrawTexture', win, facetext{ p.emo_HA_AN_NE(t), p.exshuffle(t) }, ...
        [],[p.LocXc(loc)-p.Wstim p.LocYc(loc)-p.Hstim ...
        p.LocXc(loc)+p.Wstim p.LocYc(loc)+p.Hstim], 0);

    % Create central fixation texture
    Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15);

    % Flip: Present face & central fixation (for duration = p.cuedur)
    [VBLTimestamp] = Screen(win, 'Flip');
    cue_onset = VBLTimestamp;
    WaitSecs(p.cuedur);

    % Create break time (for the duration = csi)
    if (p.csi_125_500(t) == 1)
        csi = 0.050;
    elseif (p.csi_125_500(t) == 2)
        csi = 0.425;
    end

    % Flip: Present break screen (i.e., just central fixation)
    Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15); % fixation
    csi_onset = Screen(win, 'Flip');
    WaitSecs(csi);

    % Record actual csi (just to account for possible MATLAB time variability)
    p.cuedur_real(t) = csi_onset-cue_onset;

    % This is for EEG
    if p.EEG == 1
        SerialPort.sendTrigger(t);
    end

    %% Display shape stimuli & collect response from subs

    % Initiate a few variables for updating loop/ collecting data
    key_was_pressed = 0;
    f = 0;
    p.searchdur_real(t) = 0;

    % Loop the program in a big 'while' loop
    while p.searchdur_real(t) <= p.searchdur

        % Update 'f'
        f = f+1;

        % As long as the actual search duration is less than p.stimdur (we
        % continue to present the stimuli--dimonds and bowties)
        if p.searchdur_real(t) <= p.stimdur % 1.55 sec
            for loc = 1:12
                Screen('DrawTexture', win, tgtext{p.stimshape(loc, t), 1}, ...
                    [],[p.LocX(loc)-p.Rstim p.LocY(loc)-p.Rstim ...
                    p.LocX(loc)+p.Rstim p.LocY(loc)+p.Rstim], []);
            end
        end

        % Create fixation texture
        Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15);

        % Flip to present stimuli with fixation
        [VBLTimestamp] = Screen(win, 'Flip');

        % If actual search duration = 0, we collect stimOnset time and send
        % EEG trigger (if EEG ==1)
        if p.searchdur_real(t) == 0
            if p.EEG == 1
                SerialPort.sendTrigger(101);
            end
            p.stimonset(t) = VBLTimestamp;
        end

        % The below 'if' loop is for detecting key press. If no key is pressed (key_was_pressed == 0), do the below:
        if key_was_pressed == 0

            % Keep looking for response (we use the function 'KbCheck')
            [key_was_pressed, press_time, key_list] = KbCheck;

            % If key is pressed:
            if key_was_pressed

                % Check what key is pressed (using the 'find' function)
                key_code = find(key_list);

                % If the key pressed is abort key ('ESC'), then close the
                % screen and exit psychtoolbox
                if key_code == abortKey % abort key
                    ListenChar(0);
                    ShowCursor;
                    clear screen;
                    clear sound;
                    Screen('CloseAll')

                    % If they key pressed is a member of the allowed keys ('j'
                    % or 'k'), then:
                elseif ismember(key_code(1), respCode)

                    % Send trigger (if EEG)
                    if p.EEG == 1
                        SerialPort.sendTrigger(102);
                    end

                    % Clear respindex
                    clear respindex

                    % Find respindex of the key pressed
                    respindex = find(respCode == key_code(1));

                    % Record reaction time (RT)
                    p.RT(t) = press_time-p.stimonset(t);

                    % Record what key is pressed
                    p.resp(t) = respindex ;

                    % Record whether the resposne is accurate or not
                    % (p.hit)
                    p.hit(t) = p.resp(t) == p.shape_diamondglass(t);

                    % Record p.slow and p.hit_notslow
                    p.hit_notslow(t) = p.hit(t);
                    if p.RT(t) <= p.rspdeadline
                        p.slow(t) = 0;
                    else
                        p.hit_notslow(t) = 0;
                    end

                    % If the key pressed is neither the abort key or a member of the allowed keys:
                else

                    % Set key_was_pressed = 0
                    key_was_pressed = 0;
                end
            end
        end

        % Update the actual time used for resposne searching
        p.searchdur_real(t) = GetSecs - p.stimonset(t);

    end

    % After exiting the while loop, present break screen (central fixatio)
    % for the duration = p.feecback
    Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15); % fixation
    Screen(win, 'Flip');
    WaitSecs(p.feedback);

    %% Present feedback

    % Create feedback image texture depending on whether the response was
    % slow, right or wrong
    if p.slow(t) == 1
        Screen('DrawTexture', win, slow, [], [p.xcenter-p.Rcue p.ycenter-p.Rcue p.xcenter+p.Rcue p.ycenter+p.Rcue], []); % slow
    elseif p.hit(t) == 1
        Screen('DrawTexture', win, right, [], [p.xcenter-p.Rcue p.ycenter-p.Rcue p.xcenter+p.Rcue p.ycenter+p.Rcue], []); % right
    else
        Screen('DrawTexture', win, wrong, [], [p.xcenter-p.Rcue p.ycenter-p.Rcue p.xcenter+p.Rcue p.ycenter+p.Rcue], []); % wrong
    end

    % For EEG
    if p.EEG == 1
        SerialPort.sendTrigger(104);
    end

    % Present feedback (for the duration of p.feedback)
    Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15);
    Screen(win, 'Flip');
    WaitSecs(p.feedback);

    % Present fixation (for break between trials; for the duration of p.iti(t))
    Screen('FillOval', win, black, [p.xcenter-p.Rfix p.ycenter-p.Rfix p.xcenter+p.Rfix p.ycenter+p.Rfix], 15);
    [VBLTimestamp] = Screen(win, 'Flip');
    WaitSecs(p.iti(t));

end

%% save trial data from this block

% Save 'p'
save(fName, 'p');

% Return cursor
ShowCursor;
ListenChar(0);

% Close everything
Screen('CloseAll');
clear sound;

% Display histogram of RT
hist(p.RT)

% Display total score in MATLAB console
['total score: ' num2str(round(mean(p.hit_notslow)*100,2)) '%  (block:' num2str(p.runNum(1)) ')']


%% Write table (for R)
% T = table(repmat(p.subNum,[p.numTrials,1]), repmat(p.runNum,[p.numTrials,1]), ...
%     p.gender,p.hand,p.age,...
% 
% p.shape_diamondglass', p.cueloc', p.coninc', ...
%     p.emo_HA_AN_NE', p.csi_125_500', p.exshuffle',  ...
%     p.resp', p.RT', p.hit', ...
%     p.isi', p.iti', p.hit_notslow', p.slow', ...
%     p.stimonset');
% 
% T.Properties.VariableNames = {'subNum', 'runNum',...
%     'gender', 'hand', 'age',...
%     'shape_diamondglass', 'cueloc', 'coninc',...
%     'emo_HA_AN_NE', 'csi_125_500', 'exshuffle',...
%     'resp','RT','hit',...
%     'isi','iti','hit_notslow','slow',...
%     'stimonset'};
% 
% writetable(T, [p.root '/data_envSound_csv/' 'sbj',  num2str(p.subNum), '_audType', p.audioType, '.csv'])

T = table(repmat(p.subNum,[p.numTrials,1]), repmat(p.runNum,[p.numTrials,1]),...
    p.gender,p.hand,p.age,...
    p.shape_diamondglass', p.cueloc', p.coninc',...
    p.emo_HA_AN_NE', p.csi_125_500', p.exshuffle',...
    p.resp', p.RT', p.hit', ...
    p.isi', p.iti', p.hit_notslow', p.slow', ...
    p.stimonset');

T.Properties.VariableNames = {'subNum', 'runNum',...
    'gender', 'hand', 'age',...
    'shape_diamondglass', 'cueloc', 'coninc',...
    'emo_HA_AN_NE', 'csi_125_500', 'exshuffle',...
    'resp','RT','hit',...
    'isi','iti','hit_notslow','slow',...
    'stimonset'};

writetable(T, [p.root '/data_envSound_csv/' 'sbj',  num2str(p.subNum), '_audType', p.audioType, '.csv'])