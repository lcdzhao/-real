function settings = initSettings()
%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
%
settings.msToProcess        = 1000;        %[ms]��Ҫ����ĺ�����

% Number of channels to be used for signal processing
settings.numberOfChannels   = 1;    %ͨ�����������Ǹ�����

settings.PRN = 18;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
%�ƶ����ݴ���Ŀ�ʼ�㡣�������ݼ�¼�е��κ�һ�㿪ʼ�źŴ���
%fseek�����ƶ��ļ��Ķ�ȡ��
settings.skipNumberOfBytes     = 0;  %�������ֽ���

%% Raw signal file name and other parameter ========ԭʼ�ź��ļ�������������=======================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode
settings.nco_Length = 32;                   %�������


% Data type used to store one sample
settings.dataType           = 'int8';   %�洢һ����������������

% Intermediate, sampling and code frequencies
settings.IF1                 = 9.548e6 %1.42e6 %4.123968e6;      %[Hz]   %L1��Ƶ
settings.samplingFreq       = 38.192e6 %5.714e6 %16.367667e6;     %[Hz] %����Ƶ�ʱ�
settings.codeFreqBasis      = 1.023e6;      %[Hz]   %��Ԫ�Ļ�Ƶ
settings.IF2                 = 14.548e6 %1.42e6 %4.123968e6;      %[Hz]   %L2��Ƶ
% Define number of chips in a code period
settings.codeLength         = 1023;     %һ����Ԫ���ڵġ�Ƭ����
%ÿ��CA�����ڵĲ�����������������38192
settings.samplesPerCode = round(settings.samplingFreq /(settings.codeFreqBasis/settings.codeLength));  %һ����Ԫ�ж��ٸ�������

%% Acquisition settings ==============��������=====================================
% Skips acquisition in the script postProcessing.m if set to 1
%�������Ϊ1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition    %����Ѱ�������б������ų�һЩ�����Լӿ첶��
%settings.acqSatelliteList   = 1:32;         %[PRN numbers]  %����������б�
% Band around IF to search for satellite signal. Depends on max Doppler
%����������Ƶ�ʾ���
settings.acqSearchBand      = 10;           %[kHz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2.5;  %�о���ֵ

%% Tracking loops settings =============���ٻ�·����===================================
% Code tracking loop parameters     ����ٻ�·����
settings.FLL_flag = 1;                      %FLL��־����Ϊ�տ�ʼ��FLL���Գ�ʼʱFLL�ı�־Ϊ1
settings.PLL_flag = 0;                      %�������ͬ��
settings.FLL_bandwidth = 4.2;               %FLL��������
settings.PLL_bandwidth = 10;                %PLL��������
settings.DDLL_bandwidth = 2;                %�뻷�˲���������
settings.cofe_FLL_auxi_DDLL  = 1/763;       %�ز�����ϵ��,����/�ز�Ƶ��
settings.dllCorrelatorSpacing = 0.5;



%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off
                                            % 1 - On

                                            
%% Constants ==============================================================

settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 0;       %[ms] Initial sign. travel time


%% 

settings.dup_freq = 0;                     %������Ƶ��
settings.noise_std = 1;
% setting.length = (1:10000);
% setting.length_no = 10000;
settings.sample_t = 1/settings.samplingFreq; %����ʱ��
settings.K = 1;                             %��·����
settings.transfer_coef = (2^settings.nco_Length)/settings.samplingFreq;  %Ƶ����ת��ϵ����ͬʱ����������Ƶ�ʻ�������
settings.middle_freq_nco1 = settings.IF1*settings.transfer_coef;%��Ƶ1��Ӧ��Ƶ����
settings.middle_freq_nco2 = settings.IF2*settings.transfer_coef;%��Ƶ��Ӧ��Ƶ����
settings.Ncoh = (settings.samplingFreq / settings.codeFreqBasis )*settings.codeLength;%һ���������ʱ���ڵĲ�������
settings.Tcoh = settings.Ncoh *settings.sample_t;               %�������ʱ��
settings.dot_length = [1:settings.Ncoh];                        %һ���������ʱ���ڵĲ�������
settings.code_word = settings.codeFreqBasis * settings.transfer_coef;%�뻷������
settings.fd_code = settings.dup_freq*(1/763)*settings.transfer_coef;%������ź�Դ�����ϵĶ����գ���������NCO��
%setting.fd_code = setting.dup_freq*settings.cofe_FLL_auxi_DDLL*setting.transfer_coef;%������ź�Դ�����ϵĶ����գ���������NCO��
settings.e_code_original_phase = 0;         %nco��Ƶ������ĳ���λ
settings.modulate_code_bias_phsae = 0;      %����ʱ������B1��ĳ���λ
settings.signal_phase = 0;
settings.local_phase = 0;
