function [rslt]= calculateTryFeatures(myimu, timeframe)

%% select time frame
Fs = length(myimu)/ (10 * 60);
% timeframe = 15; % in seconds
n = floor(timeframe*Fs);            %number of elements per timeframe
m = floor(size(myimu,1) / n);   %number of output elements

normVal = [1.43834995975120e-05,	2.12393816751647e-06,	4.25041556702755e-06,	3.48071618328034e-06];

%% Parameter
numaxis = 3;
if (width(myimu)>3)
myimu = myimu(:,1:numaxis);end

rslt = zeros(m, 11);
labelOcc = zeros(m,6);

for i = 1: m
    % get index of timeframe
    t1 = 1 + (i-1)*n;
    t2 = t1 +n-1;
    if i == m
        t2 = length(myimu); end
    
    data = myimu(t1:t2,:);

    %Number of fft elements
    N = pow2(ceil(log(n)/log(2))); %round to next power of 2
    
    f = Fs *(0:(N/2))/N;

    %% compute X AXIS FEATURE
    fftX = calcOneSidedSpectrum(data(:,1), N, false);
    fftX = sgolayfilt(fftX,2,101);

    % Get Power in relevant Frequency Band
    fborder = [23 27];
    buffer = fftX((f >=fborder(1)) & (f < fborder(2)));
    bufferRest = fftX((f < fborder(1)) | (f >= fborder(2)));
    Xpower25 = mean(buffer) / mean(bufferRest);     % mean of highest 10 percent
    Xpower25abs = rms(buffer);     % mean of highest 10 percent

%     fborder = [0 1000];
%     buffer = fftX((f >=fborder(1)) & (f < fborder(2)));
%     XpowerAll = rms(buffer); 

    %% compute Y AXIS FEATURE
    fftY = calcOneSidedSpectrum(data(:,2), N, false);
    fftY = sgolayfilt(fftY,2,101);

    % Get Power in relevant Frequency Band
    fborder = [62 80];
    buffer = fftY((f >=fborder(1)) & (f < fborder(2)));
    bufferRest = fftY((f < fborder(1)) | (f >= fborder(2)));
    Ypower70 = mean(buffer) / mean(bufferRest);     % mean of highest 10 percent
    Ypower70abs = rms(buffer);     % mean of highest 10 percent
%     fborder = [0 1000];
%     buffer = fftY((f >=fborder(1)) & (f < fborder(2)));
%     buffer = sort(buffer, 'descend');                          % Sort Descending
%     YpowerAll = rms(buffer); %(1:ceil(length(buffer)*0.1)));     % mean of highest 10 percent

        
    %% compute Z AXIS FEATURE
    fftZ = calcOneSidedSpectrum(data(:,3), N, false);
    fftZ = sgolayfilt(fftZ,2,101);
    % Get Power in relevant Frequency Band
    fborder = [80 120];
    buffer = fftZ((f >=fborder(1)) & (f < fborder(2)));
    bufferRest = fftZ((f < fborder(1)) | (f >= fborder(2)));
    %     buffer = sort(buffer, 'descend');                          % Sort Descending
    Zpower100 = mean(buffer) / mean(bufferRest);  
    Zpower100abs = rms(buffer);  

    fborder = [60 180];
    buffer = fftZ((f >=fborder(1)) & (f < fborder(2)));
    bufferRest = fftZ((f < fborder(1)) | (f >= fborder(2)));
    Zpower150 = mean(buffer) / mean(bufferRest);  
    Zpower150abs = rms(buffer);  

%     fborder = [0 1000];
%     buffer = fftZ((f >=fborder(1)) & (f < fborder(2)));
%     %     buffer = sort(buffer, 'descend');                          % Sort Descending
%     ZpowerAll = rms(buffer); %(1:ceil(length(buffer)*0.1)));     % mean of highest 10 percent
    normAllPower = 6e-5;
    allPowerX = rms(data(:,1),1) ./ 0.0067;
    allPowerY = rms(data(:,2),1) ./ 0.0067;
    allPowerZ = rms(data(:,3),1) ./ 0.0067;
    allPower = sum(rms(data(:,:),1)) ./ 0.0067;

    rslt(i,:) = [Xpower25,Xpower25abs, Ypower70,Ypower70abs, Zpower100,Zpower100abs, Zpower150,Zpower150abs, allPowerX, allPowerY, allPowerZ];
end
end




