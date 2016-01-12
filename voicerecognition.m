%% Project: Laryengal Vibrations
% By  Rohith K D, Deepu A D, Amaya Shaji, Treasa Annie Jose.

%% Main Function Voice Recognition
function []=voicerecognition()
% For clear screen
clc;
% min_distance is a variable used to directly set the minimum distance for word recognition.
min_distance=20;
%the message to be displayed.
disp('Project: Laryengal Vibrations');
disp('By  Rohith K D, Deepu A D, Amaya Shaji, Trease Annie Jose.');

%sampling frequency is set as 8000.
samplingfrequency=8000;
% For making Choice
ch=0;
poss=2;
while ch~=poss
    
%menu is displayed. The first parameter passed is the title, and the
%coressponding parameters passed are the butons to be displayed. The
%program have two buttons. One button for the main operation and other for
%exiting the system.

      ch=menu('Laryengal to sound conversion System',....
        '1: Test with other speech files','2: Exit');
   
%% 1: Laryengal vibraton to word conversion by letting user enter into database and then compare
    %finds which word the person want to talk.

    if ch==1
        disp('> 1: Test with laryengal vibration files')
% The secondary menu. It is for adding new vibration to database, testing
% the vibrations for appropriate words. It also have options for viewing
% database and deleting database.

        chos=0;
        possibility=5;
        while chos~=possibility,
            chos=menu('Word Identification System','Add a new vibration from your device',...
                'Test the vibration from microphone',...
                'Database Info','Delete database','Exit');

            %----------------------------------------------------------------------

%% 10.1 Add a new vibration from the device.


            if chos==1
                
                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');                   
                    classe = input('Enter the word  that will be used for recognition:','s');
                    message=('The following parameters will be used during recording:');
                    disp(message);
                    message=strcat('Bits per sample',num2str(samplingbits));
                    disp(message);
                    durata = input('Insert the duration of the recording (in seconds):');
                    if isempty(durata)
                        durata = 3;
                        disp( num2str(durata) );
                    end
                    micrecorder = audiorecorder(8000,samplingbits,1);
                    disp('Now, try speaking. The device will pickup the vibrations');
                    record(micrecorder,durata);
                    
                    while (isrecording(micrecorder)==1)
                        disp('Recording...');
                        pause(0.5);
                    end
                    disp('Recording stopped.');
                    y1 = getaudiodata(micrecorder);
                    % converting the wav file to get the integer value.
                    y = getaudiodata(micrecorder, 'uint8');

                    if size(y,2)==2
                        y=y(:,1);
                    end
                    y = double(y);
                    sound_number = sound_number+1;
                    data{sound_number,1} = y;
                    data{sound_number,2} = classe;
                    data{sound_number,3} = 'Microphone';
                    data{sound_number,4} = 'Microphone';
                    st=strcat('u',num2str(sound_number));
                    wavwrite(y1,8000,samplingbits,st)
                    save('sound_database.dat','data','sound_number','-append');
                    msgbox('Vibration added to database','Database result','help');
                    disp('Vibration added to database');

                else
                    classe = input('Enter the word that will be used for recognition:','s');
                    durata = input('Insert the duration of the recording (in seconds):');
                    if isempty(durata)
                        durata = 3;
                        disp( num2str(durata) );
                    end
                  
                    samplingbits = input('Insert the number of bits per sample (8 recommended):');
                    if isempty(samplingbits )
                        samplingbits = 8;
                        disp( num2str(samplingbits) );
                    end
                    micrecorder = audiorecorder(8000,samplingbits,1);
                    disp('Now, speak into microphone...');
                    record(micrecorder,durata);

                    while (isrecording(micrecorder)==1)
                        disp('Recording...');
                        pause(0.5);
                    end
                    disp('Recording stopped.');
                    y1 = getaudiodata(micrecorder);
                    y = getaudiodata(micrecorder, 'uint8');

                    if size(y,2)==2
                        y=y(:,1);
                    end
                    y = double(y);
                    sound_number = 1;
                    data{sound_number,1} = y;
                    data{sound_number,2} = classe;
                    data{sound_number,3} = 'Microphone';
                    data{sound_number,4} = 'Microphone';
                    st=strcat('u',num2str(sound_number));
                    wavwrite(y1,8000,samplingbits,st)
                    save('sound_database.dat','data','sound_number','samplingfrequency','samplingbits');
                    msgbox('Vibration added to database','Database result','help');
                    disp('Vibration added to database');
                end
            end

            %----------------------------------------------------------------------

%% 10.2 Word Recognition from Device

            if chos==2

                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');
                    Fs = 8000;
                    durata = input('Insert the duration of the recording (in seconds):');
                    if isempty(durata)
                        durata = 3;
                        disp( num2str(durata) );
                    end
                    micrecorder = audiorecorder(8000,samplingbits,1);
                    disp('Now, try to speak with device...');
                    record(micrecorder,durata);

                    while (isrecording(micrecorder)==1)
                        disp('Recording...');
                        pause(0.5);
                    end
                    disp('Recording stopped.');
                    y = getaudiodata(micrecorder);
                    st='v';
                    wavwrite(y,8000,samplingbits,st);
                    y = getaudiodata(micrecorder, 'uint8');
                    % if the input sound is not mono

                    if size(y,2)==2
                        y=y(:,1);
                    end
                    y = double(y);
                    %----- code for speaker recognition -------
                    disp('MFCC cofficients computation and VQ codebook training in progress...');
                    disp(' ');
                    % Number of centroids required
                    k =16;

                    for ii=1:sound_number
                        % Compute MFCC cofficients for each sound present in database
                        v = mfcc(data{ii,1}, Fs);
                        % Train VQ codebook
                        code{ii} = vqlbg(v, k);
                        disp('...');
                    end
                    disp('Completed.');
                    % Compute MFCC coefficients for input sound
                    v = mfcc(y,Fs);
                    % Current distance and sound ID initialization
                    distmin = Inf;
                    k1 = 0;

                    for ii=1:sound_number
                        d = disteu(v, code{ii});
                        dist = sum(min(d,[],2)) / size(d,1);
                        message=strcat('For word #',num2str(ii),' Dist : ',num2str(dist));
                        disp(message);
             
                        if dist < distmin
                            distmin = dist;
                            k1 = ii;
                        end
                    end

                    if distmin < min_distance
                        min_index = k1;
                        speech_id = data{min_index,2};
                        %-----------------------------------------
                        disp('Matching sound:');
                        message=strcat('File:',data{min_index,3});
                        disp(message);
                        message=strcat('Location:',data{min_index,4});
                        disp(message);
                        message = strcat('Recognized word ID: ',num2str(speech_id));
                        disp(message);
                        msgbox(message,'Matching result','help');
                        

                        %speech production output. Done using Microsoft
                        %SAPI 5.0
                      textIn = speech_id;
                        ha = actxserver('SAPI.SpVoice');
                        invoke(ha,'speak',textIn);



                        ch3=0;
                        while ch3~=3
                            ch3=menu('Matched result verification:','Recognized Sound','Recorded sound','Exit');

                            if ch3==1
                                st=strcat('u',num2str(speech_id));
                                [s fs nb]=wavread(st);
                                p=audioplayer(s,fs,nb);
                                play(p);
                            end

                            if ch3==2
                                [s fs nb]=wavread('v');
                                p=audioplayer(s,fs,nb);
                                play(p);
                            end
                        end

                    else
                        warndlg('No matching Result.',' Warning ')
                    end
                else
                    warndlg('Database is empty. No word detection is possible.',' Warning ')
                end
            end
            %----------------------------------------------------------------------

%% 10.3 Database Info

            if chos==3

                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');
                    message=strcat('Database has #',num2str(sound_number),'words:');
                    disp(message);
                    disp(' ');

                    for ii=1:sound_number
                        message=strcat('Location:',data{ii,3});
                        disp(message);
                        message=strcat('File:',data{ii,4});
                        disp(message);
                        message=strcat('Sound ID:',num2str(data{ii,2}));
                        disp(message);
                        disp('-');
                    end

                    ch32=0;
                    while ch32 ~=2
                        ch32=menu('Database Information','Database','Exit');

                        if ch32==1
                            st=strcat('Sound Database has : #',num2str(sound_number),'words. Enter a database number : #');
                            prompt = {st};
                            dlg_title = 'Database Information';
                            num_lines = 1;
                            def = {'1'};
                            options.Resize='on';
                            options.WindowStyle='normal';
                            options.Interpreter='tex';
                            an = inputdlg(prompt,dlg_title,num_lines,def);
                            an=cell2mat(an);
                            a=str2double(an);
                            
                            if (isempty(an))

                            else

                                if (a <= sound_number)
                                    st=strcat('u',num2str(an));
                                    [s fs nb]=wavread(st);
                                    p=audioplayer(s,fs,nb);
                                    play(p);
                                else
                                    warndlg('Invalid Word ','Warning');
                                end
                            end
                        end
                    end

                else
                    warndlg('Database is empty.',' Warning ')
                end
            end
            %----------------------------------------------------------------------

%% 10.4 Delete database

            if chos==4
                %clc;
                close all;

                if (exist('sound_database.dat','file')==2)
                    button = questdlg('Do you really want to remove the Database?');

                    if strcmp(button,'Yes')
                        load('sound_database.dat','-mat');

                        for ii=1:sound_number
                            st=strcat('u',num2str(ii),'.wav');
                            delete(st);
                        end

                        if (exist('v.wav','file')==2)
                            delete('v.wav');
                        end

                        delete('sound_database.dat');
                        msgbox('Database was succesfully removed from the current directory.','Database removed','help');
                    end

                else
                    warndlg('Database is empty.',' Warning ')
                end
            end
        end
    end
close all;
end

close all;
%--------------------------------------------------------------------------

%% blockFrames Function
% blockFrames: Puts the signal into frames
%
% Inputs: s contains the signal to analize
% fs is the sampling rate of the signal
% m is the distance between the beginnings of two frames
% n is the number of samples per frame
%
% Output: M3 is a matrix containing all the frames

function M3 = blockFrames(s, fs, m, n)
l = length(s);
nbFrame = floor((l - n) / m) + 1;
for i = 1:n
    for j = 1:nbFrame
        M(i, j) = s(((j - 1) * m) + i);
    end
end
h = hamming(n);
M2 = diag(h) * M;
for i = 1:nbFrame
    M3(:, i) = fft(M2(:, i));
end
%--------------------------------------------------------------------------

%% MFCC Function
% MFCC
%
% Inputs: s contains the signal to analize
% fs is the sampling rate of the signal
%
% Output: r contains the transformed signal

function r = mfcc(s, fs)
m = 100;
n = 256;
frame=blockFrames(s, fs, m, n);
m = melfb(20, n, fs);
n2 = 1 + floor(n / 2);
z = m * abs(frame(1:n2, :)).^2;
r = dct(log(z));

%--------------------------------------------------------------------------

%% VQLBG Vector quantization using the Linde-Buzo-Gray algorithm
% VQLBG Vector quantization using the Linde-Buzo-Gray algorithm
%
% Inputs: d contains training data vectors (one per column)
% k is number of centroids required
%
% Output: r contains the result VQ codebook (k columns, one for each  centroids)

function r = vqlbg(d,k)
e = .01;
r = mean(d, 2);
dpr = 10000;
for i = 1:log2(k)
    r = [r*(1+e), r*(1-e)];
    while (1 == 1)
        z = disteu(d, r);
        [m,ind] = min(z, [], 2);
        t = 0;
        for j = 1:2^i
            r(:, j) = mean(d(:, find(ind == j)), 2); 
            x = disteu(d(:, find(ind == j)), r(:, j)); 
            for q = 1:length(x)
                t = t + x(q);
            end
        end
        if (((dpr - t)/t) < e)
            break;
        else
            dpr = t;
        end
    end
end