function settings = initSettings()
%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
%
settings.msToProcess        = 1000;        %[ms]需要处理的毫秒数

% Number of channels to be used for signal processing
settings.numberOfChannels   = 1;    %通道数（即卫星个数）

settings.PRN = 18;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
%移动数据处理的开始点。能在数据记录中的任何一点开始信号处理。
%fseek函数移动文件的读取点
settings.skipNumberOfBytes     = 0;  %跳过的字节数

%% Raw signal file name and other parameter ========原始信号文件名和其它参数=======================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode
settings.nco_Length = 32;                   %方便计算


% Data type used to store one sample
settings.dataType           = 'int8';   %存储一个采样的数据类型

% Intermediate, sampling and code frequencies
settings.IF1                 = 9.548e6 %1.42e6 %4.123968e6;      %[Hz]   %L1中频
settings.samplingFreq       = 38.192e6 %5.714e6 %16.367667e6;     %[Hz] %采样频率比
settings.codeFreqBasis      = 1.023e6;      %[Hz]   %码元的基频
settings.IF2                 = 14.548e6 %1.42e6 %4.123968e6;      %[Hz]   %L2中频
% Define number of chips in a code period
settings.codeLength         = 1023;     %一个码元周期的“片”数
%每个CA码周期的采样数，整数倍不好38192
settings.samplesPerCode = round(settings.samplingFreq /(settings.codeFreqBasis/settings.codeLength));  %一个码元有多少个采样点

%% Acquisition settings ==============捕获设置=====================================
% Skips acquisition in the script postProcessing.m if set to 1
%如果设置为1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition    %所搜寻的卫星列表，可以排除一些卫星以加快捕获
%settings.acqSatelliteList   = 1:32;         %[PRN numbers]  %捕获的卫星列表
% Band around IF to search for satellite signal. Depends on max Doppler
%由最大多普勒频率决定
settings.acqSearchBand      = 10;           %[kHz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2.5;  %判决阈值

%% Tracking loops settings =============跟踪环路设置===================================
% Code tracking loop parameters     码跟踪环路参数
settings.FLL_flag = 1;                      %FLL标志，因为刚开始用FLL所以初始时FLL的标志为1
settings.PLL_flag = 0;                      %和上面的同理
settings.FLL_bandwidth = 4.2;               %FLL噪声带宽
settings.PLL_bandwidth = 10;                %PLL噪声带宽
settings.DDLL_bandwidth = 2;                %码环滤波噪声带宽
settings.cofe_FLL_auxi_DDLL  = 1/763;       %载波辅助系数,码率/载波频率
settings.dllCorrelatorSpacing = 0.5;



%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off
                                            % 1 - On

                                            
%% Constants ==============================================================

settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 0;       %[ms] Initial sign. travel time


%% 

settings.dup_freq = 0;                     %多普勒频率
settings.noise_std = 1;
% setting.length = (1:10000);
% setting.length_no = 10000;
settings.sample_t = 1/settings.samplingFreq; %采样时间
settings.K = 1;                             %环路增益
settings.transfer_coef = (2^settings.nco_Length)/settings.samplingFreq;  %频率字转换系数，同时，将采样的频率还在其中
settings.middle_freq_nco1 = settings.IF1*settings.transfer_coef;%中频1对应的频率字
settings.middle_freq_nco2 = settings.IF2*settings.transfer_coef;%中频对应的频率字
settings.Ncoh = (settings.samplingFreq / settings.codeFreqBasis )*settings.codeLength;%一个积分清除时间内的采样点数
settings.Tcoh = settings.Ncoh *settings.sample_t;               %积分清除时间
settings.dot_length = [1:settings.Ncoh];                        %一个积分清除时间内的采样点数
settings.code_word = settings.codeFreqBasis * settings.transfer_coef;%码环控制字
settings.fd_code = settings.dup_freq*(1/763)*settings.transfer_coef;%添加在信号源的码上的多普勒，体现在码NCO上
%setting.fd_code = setting.dup_freq*settings.cofe_FLL_auxi_DDLL*setting.transfer_coef;%添加在信号源的码上的多普勒，体现在码NCO上
settings.e_code_original_phase = 0;         %nco扩频码早码的初相位
settings.modulate_code_bias_phsae = 0;      %调制时，本地B1码的初相位
settings.signal_phase = 0;
settings.local_phase = 0;
