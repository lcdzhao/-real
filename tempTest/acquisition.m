function acqResult = acquisition(longSignal, settings)
%% 初始化 ============================================================

%--- 求每个周期CA码的采样点总数 ---------------------
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

%--- 产生两1ms的数据矢量 ----------------------------
signal1 = longSignal(1 : samplesPerCode);
signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);

%--- 信号 -------------------------------------------
signal0DC = longSignal - mean(longSignal);      

%--- 求采样周期 -------------------------------------
ts = 1 / settings.samplingFreq;                         

%--- 相位点 -----------------------------------------
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;   

%--- 频率带宽搜索的次数 -----------------------------
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;   

%--- 产生CA码 ---------------------------------------
caCodeTable = makeCaTable(0,settings.PRN,...
    settings.codeLength,settings.codeFreqBasis ,settings.samplingFreq,settings);

%--- 初始化结果集 -----------------------------------
results     = zeros(numberOfFrqBins, samplesPerCode);

%--- 初始化载波频率搜索的频率 -----------------------
frqBins     = zeros(1, numberOfFrqBins);



%% 相关信号 ===========================================================   
%--- CA码的DFT------------------------------------------
caCodeFreqDom = conj(fft(caCodeTable));


%--- 对整个频段进行相关运算 ----------------------------
for frqBinIndex = 1:numberOfFrqBins

    %--- 生成载波频率网格（0.5KHz步进） ----------------
    frqBins(frqBinIndex) = settings.IF1 - ...
                           (settings.acqSearchBand/2) * 1000 + ...
                           0.5e3 * (frqBinIndex - 1);

    %--- 产生本地的sin和cos ---------------------------
    sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
    cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

    %--- 点乘 -----------------------------------------
    I1      = sinCarr .* signal1;
    Q1      = cosCarr .* signal1;
    I2      = sinCarr .* signal2;
    Q2      = cosCarr .* signal2;

    %--- 将基带信号转换为频域 -------------------------
    IQfreqDom1 = fft(I1 + j*Q1);
    IQfreqDom2 = fft(I2 + j*Q2);

    %--- 频域乘法（时域相关）--------------------------
    convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
    convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;

    %--- 执行反向DFT并存储相关结果 --------------------
    acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
    acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;

    %--- 检查哪个MSEC具有更大的功率并保存--------------
    %--- 将混合第一个和第二个MSEC，但将纠正数据位问题--
    if (max(acqRes1) > max(acqRes2))
        results(frqBinIndex, :) = acqRes1;
    else
        results(frqBinIndex, :) = acqRes2;
    end

end % frqBinIndex = 1:numberOfFrqBins

%% 在结果中查找相关峰值 ===============================================
% 找到最高峰，并将其与第二高峰进行比较
% 第二个峰值与最高峰值的距离不小于1个码片

%--- 找出相关峰和载波频率 -------------------------------------
[peakSize frequencyBinIndex] = max(max(results, [], 2));

%--- 寻找相同相关峰的码相 -------------------------------------
[peakSize codePhase] = max(max(results));

%--- 查找1个码片范围的C/A码相位排除范围 -----------------------
samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
excludeRangeIndex1 = codePhase - samplesPerCodeChip;
excludeRangeIndex2 = codePhase + samplesPerCodeChip;

%--- 如果范围包括数组边界，则更正C/A代码相位排除范围 ----------
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

%--- 在同一频率中找到第二个最高相关峰 ------------------------
secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

%--- 存储结果 ------------------------------------------------
acqResult.peakMetric = peakSize/secondPeakSize;

% 如果结果高于阈值，则存在信号 
if (peakSize/secondPeakSize) > settings.acqThreshold
     acqResult.carrFreq  = frqBins(frequencyBinIndex);
     acqResult.codePhase = codePhase;
end

