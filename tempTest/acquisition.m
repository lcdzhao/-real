function acqResult = acquisition(longSignal, settings)
%% ��ʼ�� ============================================================

%--- ��ÿ������CA��Ĳ��������� ---------------------
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

%--- ������1ms������ʸ�� ----------------------------
signal1 = longSignal(1 : samplesPerCode);
signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);

%--- �ź� -------------------------------------------
signal0DC = longSignal - mean(longSignal);      

%--- ��������� -------------------------------------
ts = 1 / settings.samplingFreq;                         

%--- ��λ�� -----------------------------------------
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;   

%--- Ƶ�ʴ��������Ĵ��� -----------------------------
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;   

%--- ����CA�� ---------------------------------------
caCodeTable = makeCaTable(0,settings.PRN,...
    settings.codeLength,settings.codeFreqBasis ,settings.samplingFreq,settings);

%--- ��ʼ������� -----------------------------------
results     = zeros(numberOfFrqBins, samplesPerCode);

%--- ��ʼ���ز�Ƶ��������Ƶ�� -----------------------
frqBins     = zeros(1, numberOfFrqBins);



%% ����ź� ===========================================================   
%--- CA���DFT------------------------------------------
caCodeFreqDom = conj(fft(caCodeTable));


%--- ������Ƶ�ν���������� ----------------------------
for frqBinIndex = 1:numberOfFrqBins

    %--- �����ز�Ƶ������0.5KHz������ ----------------
    frqBins(frqBinIndex) = settings.IF1 - ...
                           (settings.acqSearchBand/2) * 1000 + ...
                           0.5e3 * (frqBinIndex - 1);

    %--- �������ص�sin��cos ---------------------------
    sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
    cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

    %--- ��� -----------------------------------------
    I1      = sinCarr .* signal1;
    Q1      = cosCarr .* signal1;
    I2      = sinCarr .* signal2;
    Q2      = cosCarr .* signal2;

    %--- �������ź�ת��ΪƵ�� -------------------------
    IQfreqDom1 = fft(I1 + j*Q1);
    IQfreqDom2 = fft(I2 + j*Q2);

    %--- Ƶ��˷���ʱ����أ�--------------------------
    convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
    convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;

    %--- ִ�з���DFT���洢��ؽ�� --------------------
    acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
    acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;

    %--- ����ĸ�MSEC���и���Ĺ��ʲ�����--------------
    %--- ����ϵ�һ���͵ڶ���MSEC��������������λ����--
    if (max(acqRes1) > max(acqRes2))
        results(frqBinIndex, :) = acqRes1;
    else
        results(frqBinIndex, :) = acqRes2;
    end

end % frqBinIndex = 1:numberOfFrqBins

%% �ڽ���в�����ط�ֵ ===============================================
% �ҵ���߷壬��������ڶ��߷���бȽ�
% �ڶ�����ֵ����߷�ֵ�ľ��벻С��1����Ƭ

%--- �ҳ���ط���ز�Ƶ�� -------------------------------------
[peakSize frequencyBinIndex] = max(max(results, [], 2));

%--- Ѱ����ͬ��ط������ -------------------------------------
[peakSize codePhase] = max(max(results));

%--- ����1����Ƭ��Χ��C/A����λ�ų���Χ -----------------------
samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
excludeRangeIndex1 = codePhase - samplesPerCodeChip;
excludeRangeIndex2 = codePhase + samplesPerCodeChip;

%--- �����Χ��������߽磬�����C/A������λ�ų���Χ ----------
if excludeRangeIndex1 < 2
    codePhaseRange = excludeRangeIndex2 : ...
                     (samplesPerCode + excludeRangeIndex1);

elseif excludeRangeIndex2 >= samplesPerCode
    codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                     excludeRangeIndex1;
else
    codePhaseRange = [1:excludeRangeIndex1, ...
                      excludeRangeIndex2 : samplesPerCode];
end

%--- ��ͬһƵ�����ҵ��ڶ��������ط� ------------------------
secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

%--- �洢��� ------------------------------------------------
acqResult.peakMetric = peakSize/secondPeakSize;

% ������������ֵ��������ź� 
if (peakSize/secondPeakSize) > settings.acqThreshold
     acqResult.carrFreq  = frqBins(frequencyBinIndex);
     acqResult.codePhase = codePhase;
end

