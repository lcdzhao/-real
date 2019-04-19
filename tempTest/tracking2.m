function trackResult2 = tracking2(fid, acqResult, settings,data)
% 这里和tracking1相同，因此不再细加注释
%% 初始化 ----------------------------------------------------------
loopPara = loopCanshuCalculate(settings);%计算环路滤波器参数

%codeTable = cacode(settings.PRN);          %调用函数，产生伪随机码

fllNcoAdder = 0;                   %fll  NCO加法器，应该在滤波器中使用
carrierNcoSum = 0;                 %给这个值再乘2*pi就是相位了
pllNcoAdder = 0;                   %pll  NCO加法器
loopCount = 0;                      %循环次数

nIQ = 2;                            
n = 3;
outputFll(2:3) = 0;                %fll鉴频器的输出
outputFilterFll(1:3) = 0;          %fll 环路滤波器的输出
outputFilterPll(1:3) = 0;          %pll 环路滤波器的输出
outputPll(2:3) = 0;                %pll 鉴频器的输出

trackResult2.pllDiscrFilter = zeros(1,settings.msToProcess);



Tcoh = settings.Tcoh;          %积分清零时间


trackResult2.carrNcoPhases = zeros(1,settings.msToProcess);       %每次积分清零前载波的nco相位,未经过转换
trackResult2.carrFreq = zeros(1,settings.msToProcess);
trackResult2.trackFlag = 0;                  %捕获成功标志位
blksize = round(settings.samplingFreq * settings.Tcoh);
startCountPhase = -200;
carrStartPhaseSum = 0;
for loopNum = 1 : settings.msToProcess
    carrNcoPhase = mod(carrierNcoSum,2^settings.ncoLength);  
    trackResult2.carrNcoPhases(loopNum) = carrNcoPhase;
    trackResult2.carrFreq(loopNum) = ...
         (settings.middleFreqNco2 + fllNcoAdder + pllNcoAdder)/settings.transferCoef;
    trackResult2.flag(loopNum) = settings.PLLFlag;         %标识该次循环有没有进行PLL锁定
   
    %读取接收数据，因为数据是循环的，所以这里是这样   
    receiveSignal = data(fid:fid + blksize - 1);
 %   fid = fid + blksize ;

    if 1 == settings.PLLFlag
        startCountPhase = startCountPhase + 1;
        if startCountPhase >= 1
            carrStartPhaseSum = carrStartPhaseSum + trackResult2.carrNcoPhases(loopNum);
        end
    else     
        if startCountPhase >= -5
            startCountPhase = -10;
            carrStartPhaseSum = 0;

        end
    end
    
    %产生本地再生载波
    for demondNum = 1:settings.Ncoh 
        localCos(demondNum) = cos(2*pi*carrierNcoSum/2^settings.ncoLength);
        localSin(demondNum) = -sin(2*pi*carrierNcoSum/2^settings.ncoLength);
        carrierNcoSum = carrierNcoSum + settings.middleFreqNco2 + fllNcoAdder + pllNcoAdder ;%本地再生载波NCO
    end
    IDemonCarrier = localCos.*receiveSignal;
    QDemonCarrier = localSin.*receiveSignal;
    
    %信号解扩并积分清除

    I_P_final(nIQ) = sum(IDemonCarrier);
    Q_P_final(nIQ) = sum(QDemonCarrier);

    
    
    if  1 == loopNum
        I_P_final(nIQ - 1) = I_P_final(nIQ);
        Q_P_final(nIQ - 1) = Q_P_final(nIQ);
    else
        if settings.PLLFlag == 0
            %         四象限反正切鉴频器
            dotFll = I_P_final(nIQ - 1) * I_P_final(nIQ) + Q_P_final(nIQ - 1) * Q_P_final(nIQ);
            crossFll = I_P_final(nIQ - 1) * Q_P_final(nIQ) - I_P_final(nIQ) * Q_P_final(nIQ - 1);
            outputFll(n) = atan2(crossFll,dotFll)/(Tcoh*2*pi); 
            trackResult2.FllDiscr(loopNum) = outputFll(n);

            outputFilterFll(n) = (loopPara.cofeone_FLL * outputFll(n)) + (loopPara.cofetwo_FLL * outputFll(n - 1)) + (2 * outputFilterFll(n - 1)) - outputFilterFll(n - 2);
            trackResult2.fllDiscrFilt(loopNum) = outputFilterFll(n);

            fllNcoAdder = outputFilterFll(n) * settings.transferCoef ;  %频率字转换      
            outputFll(n - 1)=outputFll(n);
            outputFilterFll(n - 2)=outputFilterFll(n - 1);
            outputFilterFll(n - 1)=outputFilterFll(n);
        end
        if settings.PLLFlag == 1
            %锁相环鉴相器
            outputPll(n) = atan2(Q_P_final(nIQ),I_P_final(nIQ)); 
            outputFilterPll(n) = loopPara.cofeone_PLL*outputPll(n) + loopPara.cofetwo_PLL*outputPll(n-1)+loopPara.cofethree_PLL*outputPll(n-2)+2*outputFilterPll(n-1)-outputFilterPll(n-2);
            trackResult2.pllDiscr(loopNum) = outputPll(n);
            trackResult2.pllDiscrFilter(loopNum) = outputFilterPll(n);
            pllNcoAdder = (outputFilterPll(n)/(2*pi)) * settings.transferCoef;  %频率字转换
            outputPll(n-2) = outputPll(n-1);
            outputPll(n-1) = outputPll(n);
            outputFilterPll(n-2) = outputFilterPll(n-1);
            outputFilterPll(n-1) = outputFilterPll(n);
        end
        
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
 %码环鉴别器
end

trackResult2.carrPhase = carrStartPhaseSum / startCountPhase;

if  startCountPhase > 0
    trackResult2.trackFlag = 1;
end