function settings = initSettings()
%% �������� ====================================================

settings.msToProcess        = 230;        %[ms]��Ҫ����ĺ�����

settings.PRN = 18;                        %α��ı�ţ���ͬ��ţ�α�����в�ͬ


%% ������������ ===========================================================

settings.ncoLength = 32;                        %�ֳ���������㣬��߾���

settings.IF1                 = 250e6;            %[Hz]   L1��Ƶ
settings.freqDiff           = 40e6;              %[HZ]   Ƶ��
settings.samplingFreq       = 613.8e6;           %[Hz]   ����Ƶ�ʱ�
settings.codeFreqBasis      = 1.023e6;           %[Hz]   ��Ԫ�Ļ�Ƶ
settings.IF2                 = settings.IF1 + settings.freqDiff;   %[Hz]    L2��Ƶ


settings.codeLength         = 1023;              %һ��α�����ڵ���Ԫ��
settings.sampleT = 1/settings.samplingFreq;      %����ʱ��


%һ��α����Ԫ�ж��ٸ�������
settings.samplesPerCode = round(settings.samplingFreq /(settings.codeFreqBasis/settings.codeLength));  

%%  ��������==========================================================

settings.acqSearchBand      = 10;           %[kHz]  ����ʱ�Ĵ���

settings.acqThreshold       = 2.5;  %�о���ֵ

%% ���ٻ�·����=======================================================
settings.K = 1;                             %��·����

% DDLL��������
settings.DDLLBandwidth = 2;                %�뻷�˲���������
settings.cofeFLLAuxiDDLL  = 1/763;         %�ز�����ϵ��,����/�ز�Ƶ��
settings.dllCorrelatorLength = 18;         %��ؼ���ʱ�����

% FLL ��������
settings.FLLFlag = 1;                      %FLL��־����Ϊ�տ�ʼ��FLL���Գ�ʼʱFLL�ı�־Ϊ1
settings.FLLBandwidth = 4.2;               %FLL��������

% PLL ��������
settings.PLLFlag = 0;                      %�������ͬ��
settings.PLLBandwidth = 10;                %PLL��������




                                            
%% �������� ===========================================================

settings.c                  = 299792458;                                       % ����, [m/s]
settings.CA_Period          = (1/settings.codeFreqBasis)*settings.codeLength;  % ÿ��CA�������

%% ������Ƶ�����ü��������� =====================================================
settings.dupFreq = 0;                            %������Ƶ��
settings.noiseStd = 1;                           %��������






%% ������ת����ĸ������� =====================================================
settings.transferCoef = (2^settings.ncoLength)/settings.samplingFreq;      %Ƶ����ת��ϵ����ͬʱ����������Ƶ�ʻ�������
settings.middleFreqNco1 = settings.IF1*settings.transferCoef;              %��Ƶ1��Ӧ��Ƶ����
settings.middleFreqNco2 = settings.IF2*settings.transferCoef;              %��Ƶ��Ӧ��Ƶ����
settings.codeSplitSpace = 1;                                               %�������ʱ��ռ��һ��α�����ڵļ���֮��
settings.Ncoh =  settings.codeSplitSpace * ...
    (settings.samplingFreq / settings.codeFreqBasis )*settings.codeLength; %һ���������ʱ���ڵĲ�������
settings.Tcoh = settings.Ncoh *settings.sampleT;                           %�������ʱ��
settings.dotLength = [1:settings.Ncoh];                                    %һ���������ʱ���ڵĲ�������
settings.codeWord = settings.codeFreqBasis * settings.transferCoef;        %�뻷������
settings.fdCode = settings.dupFreq*(1/763)*settings.transferCoef;          %������ź�Դ�����ϵĶ����գ���������NCO��


%% ��ʼֵ����

settings.modulateCodeBiasPhsae = 0;                   %���Ͷ�CA��ĳ���λ
settings.signalPhase = 0;                             %���Ͷ��ز��ĳ���λ
settings.eCodeOriginalPhase = 0;                      %���ն�CA������ĳ���λ
settings.localPhase = 0;                              %���ն��ز��ĳ���λ
