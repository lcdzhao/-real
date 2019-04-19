function trackResult = tracking(fid, acqResult, settings,data)
%% 初始化 =============================================================

loopPara = loopCanshuCalculate(settings);           %计算环路滤波器参数
Tcoh = settings.Tcoh;                               %积分清零时间
blksize = settings.codeSplitSpace * ...
    settings.samplesPerCode;                        %每次处理数据

%--- 鉴相时使用 ------------------------------------------------
nIQ = 2;                            
n = 3;

%--- 产生CA码 --------------------------------------------------
codeTable = cacode(settings.PRN,settings);          %调用函数，产生伪随机码

%--- fll dll pll的初始化 -------------------------------------------
fllNcoAdder = 0;                                    %fll  NCO加法器
carrierNcoSum = 0;                                  %给这个值再乘2*pi就是相位了
pllNcoAdder = 0;                                    %pll  NCO加法器
loopCount = 0;                                      %循环次数
codeNcoSum = 0;                                     %码相位和
codeNcoAdder = 0;                                   %ddll  NCO加法器
outputFll(2:3) = 0;                                 %fll鉴频器的输出
outputFilterFll(1:3) = 0;                           %fll 环路滤波器的输出
outputFilterPll(1:3) = 0;                           %pll 环路滤波器的输出
outputPll(2:3) = 0;                                 %pll 鉴频器的输出
outputFilterDdll(1:3) = 0;                          %ddll 环路滤波器的输出

trackResult.carrNcoPhases = ...
    zeros(1,settings.msToProcess);                  %每次积分清零前载波的nco相位,未经过转换
trackResult.codeNcoPhases = ...
    zeros(1,settings.msToProcess);                  %每次积分清零前CA码的nco相位,未经过转换
trackResult.carrFreq = ...
    zeros(1,settings.msToProcess);                  %每次跟踪的载波频率
trackResult.fllDiscr = ...
     zeros(1,settings.msToProcess);                 %FLL 鉴相输出
trackResult.fllDiscrFilt = ...                      
     zeros(1,settings.msToProcess);                 %FLL 环路滤波器输出
trackResult.pllDiscr = ...
    zeros(1,settings.msToProcess);                  %PLL 鉴相输出
trackResult.pllDiscrFilter = ...
    zeros(1,settings.msToProcess);                  %PLL 环路滤波器的输出


%--- 接收端早码初始化 -----------------------------------------------------
global earlyCodeNco;      %是接收端的，因为早码的话，即时码和晚码都可以从其中而来
earlyCodeNco = ((1 - (acqResult.codePhase-settings.dllCorrelatorLength)/settings.samplesPerCode)...
    *settings.codeLength) * 2^settings.ncoLength;
earlyCodeNco = mod(earlyCodeNco,2^settings.ncoLength*settings.codeLength);
localEarlyCodeLast = localEarlycodeInitial(settings,codeTable); %产生本地超前码，接收端使用，因为早码的话，即时码和晚码都可以从其中而来


%--- 跟踪设置 ----------------------------------------------------------
trackResult.trackFlag = 0;                  %跟踪成功标志位
startCountPhase = -210;                     %进行相位收集的次数
carrStartPhaseSum = 0;                      %载波相位起始点收集的次数
codeStartPhaseSum = 0;                      %CA码相位起始点收集的次数



%% 开始跟踪,并记录结果 ================================================
for loopNum = 1 : settings.msToProcess
    %--- 计算载波这一ms的初始相位并保存 -------------------------------------
    carrNcoPhase = mod(carrierNcoSum,2^settings.ncoLength);      
    trackResult.carrNcoPhases(loopNum) = carrNcoPhase;
    
    %--- 计算CA码这一ms的初始相位并保存 -------------------------------------
    trackResult.codeNcoPhases(loopNum) = ...
        ((((earlyCodeNco)/settings.codeLength)*settings.samplesPerCode...
         -(settings.dllCorrelatorLength)*2^settings.ncoLength) /settings.samplesPerCode)*settings.codeLength...
         - mod(loopNum,1/settings.codeSplitSpace) * settings.codeLength * 2^settings.ncoLength ;       
    
    %--- 计算载波这一ms的频率，并保存 --------------------------------------- 
    trackResult.carrFreq(loopNum) = ...
         (settings.middleFreqNco1 + fllNcoAdder + pllNcoAdder)/settings.transferCoef;
    
    %--- 保存该次循环是否进行了PLL锁定 --------------------------------------
    trackResult.flag(loopNum) = settings.PLLFlag;         
    
    receiveSignal = data(fid:fid + blksize - 1);    %读取接收数据
 
    %--- 计算初始相位和，这里后面改成只计算误差小的相位，丢弃误差大的--------
    if 1 == settings.PLLFlag
        startCountPhase = startCountPhase + 1;
        if startCountPhase >= 1
            carrStartPhaseSum = carrStartPhaseSum + trackResult.carrNcoPhases(loopNum);
            codeStartPhaseSum = codeStartPhaseSum + trackResult.codeNcoPhases(loopNum);
        end
    else     
        if startCountPhase >= -5
            startCountPhase = -10;
            carrStartPhaseSum = 0;
            codeStartPhaseSum = 0;
        end
    end
    
    
    %--- 产生本地再生载波 ---------------------------------------------
    for demondNum = 1:settings.Ncoh 
        localCos(demondNum) = cos(2*pi*carrierNcoSum/2^settings.ncoLength);
        localSin(demondNum) = -sin(2*pi*carrierNcoSum/2^settings.ncoLength);
        carrierNcoSum = carrierNcoSum + settings.middleFreqNco1 + fllNcoAdder + pllNcoAdder ;%本地再生载波NCO
    end
    
    
    %--- 本地再生码环NCO,和上面的carrierNcoSum作用一样 --------------------
    codeNcoSum = codeNcoAdder + settings.codeWord ...              
        + fllNcoAdder*settings.fdCode;                

    %--- 产生本地超前，即时，滞后码 --------------------------------------
    [localEarlyCode,localPromptCode,localLateCode,settings.localPhase]=localcodeGenerate(localEarlyCodeLast,codeNcoSum,codeTable,settings);
    localEarlyCodeLast = localEarlyCode;
    
    %--- 载波解调 --------------------------------------------------------   
    IDemonCarrier = localCos.*receiveSignal;
    QDemonCarrier = localSin.*receiveSignal;

    
    %--- 信号解扩并积分清除 ----------------------------------------------
    I_E_final = sum(IDemonCarrier.*localEarlyCode);
    Q_E_final = sum(QDemonCarrier.*localEarlyCode);
    I_P_final(nIQ) = sum(IDemonCarrier.*localPromptCode);
    Q_P_final(nIQ) = sum(QDemonCarrier.*localPromptCode);
    I_L_final = sum(IDemonCarrier.*localLateCode);
    Q_L_final = sum(QDemonCarrier.*localLateCode);
    
    
    
    if  1 == loopNum
        I_P_final(nIQ - 1) = I_P_final(nIQ);
        Q_P_final(nIQ - 1) = Q_P_final(nIQ);
    else
        if settings.PLLFlag == 0
            %--- FLL四象限反正切鉴频器,进行环路滤波，并保存结果 ------------------------
            dotFll = I_P_final(nIQ - 1) * I_P_final(nIQ) + Q_P_final(nIQ - 1) * Q_P_final(nIQ);
            crossFll = I_P_final(nIQ - 1) * Q_P_final(nIQ) - I_P_final(nIQ) * Q_P_final(nIQ - 1);       
            outputFll(n) = atan2(crossFll,dotFll)/(Tcoh*2*pi); 
            trackResult.fllDiscr(loopNum) = outputFll(n);
            outputFilterFll(n) = (loopPara.cofeone_FLL * outputFll(n)) + (loopPara.cofetwo_FLL * outputFll(n - 1)) + (2 * outputFilterFll(n - 1)) - outputFilterFll(n - 2);
            trackResult.fllDiscrFilt(loopNum) = outputFilterFll(n);
            fllNcoAdder = outputFilterFll(n) * settings.transferCoef ;  %频率字转换      
            outputFll(n - 1)=outputFll(n);
            outputFilterFll(n - 2)=outputFilterFll(n - 1);
            outputFilterFll(n - 1)=outputFilterFll(n);
        end
        if settings.PLLFlag == 1
             %--- FLL鉴相，进行环路滤波，并保存结果 ------------------------
            outputPll(n) = atan2(Q_P_final(nIQ),I_P_final(nIQ)); 
            outputFilterPll(n) = loopPara.cofeone_PLL*outputPll(n) + loopPara.cofetwo_PLL*outputPll(n-1)+loopPara.cofethree_PLL*outputPll(n-2)+2*outputFilterPll(n-1)-outputFilterPll(n-2);
            trackResult.pllDiscr(loopNum) = outputPll(n);
            trackResult.pllDiscrFilter(loopNum) = outputFilterPll(n);
            pllNcoAdder = (outputFilterPll(n)/(2*pi)) * settings.transferCoef;  %频率字转换
            outputPll(n-2) = outputPll(n-1);
            outputPll(n-1) = outputPll(n);
            outputFilterPll(n-2) = outputFilterPll(n-1);
            outputFilterPll(n-1) = outputFilterPll(n);
        end
        %--- 记录上次的值 ---------------------------------------------------
        I_P_final(nIQ - 1) = I_P_final(nIQ);
        Q_P_final(nIQ - 1) = Q_P_final(nIQ);
        
        if 0 == settings.PLLFlag  && abs(outputFll(n))<10  %锁频环工作状态下，信号与本地频差小于10时
            loopCount = loopCount + 1;
            if  loopCount>20            
                   settings.PLLFlag = 1;
            end
        elseif  1 == settings.PLLFlag && abs(outputFll(n))>30      %在锁相环工作状态下，锁频环所鉴出的信号与本地频差大于30时
            loopCount = loopCount-1;
            if  0 == loopCount
                settings.PLLFlag = 0;
            end
        end
    end
    
    %--- 码环鉴别，环路滤波，并保存结果 ---------------------------------------
    outputDdll(n) = ((I_E_final - I_L_final)*I_P_final(nIQ) + (Q_E_final - Q_L_final)*Q_P_final(nIQ) )/((I_P_final(nIQ)^2 + Q_P_final(nIQ)^2)*2);  % DDLL_discri_1      
    trackResult.dllDiscr(loopNum) = outputDdll(n);
    %码环滤波器（二阶）
    outputFilterDdll(n) = outputFilterDdll(n -1) + (loopPara.cofeone_DDLL*outputDdll(n)) + loopPara.cofetwo_DDLL*outputDdll(n - 1);
    trackResult.dllDiscrFilter(loopNum) = outputFilterDdll(n);
    % 转换成频率控制字
    codeNcoAdder = outputFilterDdll(n) * settings.transferCoef ; %频率字转换
    outputDdll(n - 1)=outputDdll(n);
    outputFilterDdll(n - 1) = outputFilterDdll(n);
    
    
end

%--- 计算相位 ---------------------------------------------------------------
trackResult.carrPhase = carrStartPhaseSum / startCountPhase;
trackResult.codePhase = codeStartPhaseSum / startCountPhase;
if  startCountPhase > 0
    trackResult.trackFlag = 1;
end