function settings = initSettings()
%% 处理设置 ====================================================

settings.msToProcess        = 230;        %[ms]需要处理的毫秒数

settings.PRN = 18;                        %伪码的编号，不同编号，伪码序列不同


%% 基本参数设置 ===========================================================

settings.ncoLength = 32;                        %字长，方便计算，提高精度

settings.IF1                 = 250e6;            %[Hz]   L1中频
settings.freqDiff           = 40e6;              %[HZ]   频差
settings.samplingFreq       = 613.8e6;           %[Hz]   采样频率比
settings.codeFreqBasis      = 1.023e6;           %[Hz]   码元的基频
settings.IF2                 = settings.IF1 + settings.freqDiff;   %[Hz]    L2中频


settings.codeLength         = 1023;              %一个伪码周期的码元数
settings.sampleT = 1/settings.samplingFreq;      %采样时间


%一个伪码码元有多少个采样点
settings.samplesPerCode = round(settings.samplingFreq /(settings.codeFreqBasis/settings.codeLength));  

%%  捕获设置==========================================================

settings.acqSearchBand      = 10;           %[kHz]  捕获时的带宽

settings.acqThreshold       = 2.5;  %判决阈值

%% 跟踪环路设置=======================================================
settings.K = 1;                             %环路增益

% DDLL参数设置
settings.DDLLBandwidth = 2;                %码环滤波噪声带宽
settings.cofeFLLAuxiDDLL  = 1/763;         %载波辅助系数,码率/载波频率
settings.dllCorrelatorLength = 18;         %相关计算时的码距

% FLL 参数设置
settings.FLLFlag = 1;                      %FLL标志，因为刚开始用FLL所以初始时FLL的标志为1
settings.FLLBandwidth = 4.2;               %FLL噪声带宽

% PLL 参数设置
settings.PLLFlag = 0;                      %和上面的同理
settings.PLLBandwidth = 10;                %PLL噪声带宽




                                            
%% 常数设置 ===========================================================

settings.c                  = 299792458;                                       % 光速, [m/s]
settings.CA_Period          = (1/settings.codeFreqBasis)*settings.codeLength;  % 每个CA码的周期

%% 多普勒频移设置及噪声设置 =====================================================
settings.dupFreq = 0;                            %多普勒频率
settings.noiseStd = 1;                           %噪声幅度






%% 经过字转化后的各个参数 =====================================================
settings.transferCoef = (2^settings.ncoLength)/settings.samplingFreq;      %频率字转换系数，同时，将采样的频率还在其中
settings.middleFreqNco1 = settings.IF1*settings.transferCoef;              %中频1对应的频率字
settings.middleFreqNco2 = settings.IF2*settings.transferCoef;              %中频对应的频率字
settings.codeSplitSpace = 1;                                               %积分清除时间占有一个伪码周期的几分之几
settings.Ncoh =  settings.codeSplitSpace * ...
    (settings.samplingFreq / settings.codeFreqBasis )*settings.codeLength; %一个积分清除时间内的采样点数
settings.Tcoh = settings.Ncoh *settings.sampleT;                           %积分清除时间
settings.dotLength = [1:settings.Ncoh];                                    %一个积分清除时间内的采样点数
settings.codeWord = settings.codeFreqBasis * settings.transferCoef;        %码环控制字
settings.fdCode = settings.dupFreq*(1/763)*settings.transferCoef;          %添加在信号源的码上的多普勒，体现在码NCO上


%% 初始值设置

settings.modulateCodeBiasPhsae = 0;                   %发送端CA码的初相位
settings.signalPhase = 0;                             %发送端载波的初相位
settings.eCodeOriginalPhase = 0;                      %接收端CA码早码的初相位
settings.localPhase = 0;                              %接收端载波的初相位
