function trackResult2 = tracking2(fid, acqResult, settings,data)
% �����tracking1��ͬ����˲���ϸ��ע��
%% ��ʼ�� ----------------------------------------------------------
loopPara = loopCanshuCalculate(settings);%���㻷·�˲�������

%codeTable = cacode(settings.PRN);          %���ú���������α�����

fllNcoAdder = 0;                   %fll  NCO�ӷ�����Ӧ�����˲�����ʹ��
carrierNcoSum = 0;                 %�����ֵ�ٳ�2*pi������λ��
pllNcoAdder = 0;                   %pll  NCO�ӷ���
loopCount = 0;                      %ѭ������

nIQ = 2;                            
n = 3;
outputFll(2:3) = 0;                %fll��Ƶ�������
outputFilterFll(1:3) = 0;          %fll ��·�˲��������
outputFilterPll(1:3) = 0;          %pll ��·�˲��������
outputPll(2:3) = 0;                %pll ��Ƶ�������

trackResult2.pllDiscrFilter = zeros(1,settings.msToProcess);



Tcoh = settings.Tcoh;          %��������ʱ��


trackResult2.carrNcoPhases = zeros(1,settings.msToProcess);       %ÿ�λ�������ǰ�ز���nco��λ,δ����ת��
trackResult2.carrFreq = zeros(1,settings.msToProcess);
trackResult2.trackFlag = 0;                  %����ɹ���־λ
blksize = round(settings.samplingFreq * settings.Tcoh);
startCountPhase = -200;
carrStartPhaseSum = 0;
for loopNum = 1 : settings.msToProcess
    carrNcoPhase = mod(carrierNcoSum,2^settings.ncoLength);  
    trackResult2.carrNcoPhases(loopNum) = carrNcoPhase;
    trackResult2.carrFreq(loopNum) = ...
         (settings.middleFreqNco2 + fllNcoAdder + pllNcoAdder)/settings.transferCoef;
    trackResult2.flag(loopNum) = settings.PLLFlag;         %��ʶ�ô�ѭ����û�н���PLL����
   
    %��ȡ�������ݣ���Ϊ������ѭ���ģ���������������   
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
    
    %�������������ز�
    for demondNum = 1:settings.Ncoh 
        localCos(demondNum) = cos(2*pi*carrierNcoSum/2^settings.ncoLength);
        localSin(demondNum) = -sin(2*pi*carrierNcoSum/2^settings.ncoLength);
        carrierNcoSum = carrierNcoSum + settings.middleFreqNco2 + fllNcoAdder + pllNcoAdder ;%���������ز�NCO
    end
    IDemonCarrier = localCos.*receiveSignal;
    QDemonCarrier = localSin.*receiveSignal;
    
    %�źŽ������������

    I_P_final(nIQ) = sum(IDemonCarrier);
    Q_P_final(nIQ) = sum(QDemonCarrier);

    
    
    if  1 == loopNum
        I_P_final(nIQ - 1) = I_P_final(nIQ);
        Q_P_final(nIQ - 1) = Q_P_final(nIQ);
    else
        if settings.PLLFlag == 0
            %         �����޷����м�Ƶ��
            dotFll = I_P_final(nIQ - 1) * I_P_final(nIQ) + Q_P_final(nIQ - 1) * Q_P_final(nIQ);
            crossFll = I_P_final(nIQ - 1) * Q_P_final(nIQ) - I_P_final(nIQ) * Q_P_final(nIQ - 1);
            outputFll(n) = atan2(crossFll,dotFll)/(Tcoh*2*pi); 
            trackResult2.FllDiscr(loopNum) = outputFll(n);

            outputFilterFll(n) = (loopPara.cofeone_FLL * outputFll(n)) + (loopPara.cofetwo_FLL * outputFll(n - 1)) + (2 * outputFilterFll(n - 1)) - outputFilterFll(n - 2);
            trackResult2.fllDiscrFilt(loopNum) = outputFilterFll(n);

            fllNcoAdder = outputFilterFll(n) * settings.transferCoef ;  %Ƶ����ת��      
            outputFll(n - 1)=outputFll(n);
            outputFilterFll(n - 2)=outputFilterFll(n - 1);
            outputFilterFll(n - 1)=outputFilterFll(n);
        end
        if settings.PLLFlag == 1
            %���໷������
            outputPll(n) = atan2(Q_P_final(nIQ),I_P_final(nIQ)); 
            outputFilterPll(n) = loopPara.cofeone_PLL*outputPll(n) + loopPara.cofetwo_PLL*outputPll(n-1)+loopPara.cofethree_PLL*outputPll(n-2)+2*outputFilterPll(n-1)-outputFilterPll(n-2);
            trackResult2.pllDiscr(loopNum) = outputPll(n);
            trackResult2.pllDiscrFilter(loopNum) = outputFilterPll(n);
            pllNcoAdder = (outputFilterPll(n)/(2*pi)) * settings.transferCoef;  %Ƶ����ת��
            outputPll(n-2) = outputPll(n-1);
            outputPll(n-1) = outputPll(n);
            outputFilterPll(n-2) = outputFilterPll(n-1);
            outputFilterPll(n-1) = outputFilterPll(n);
        end
        
        I_P_final(nIQ - 1) = I_P_final(nIQ);
        Q_P_final(nIQ - 1) = Q_P_final(nIQ);
       if 0 == settings.PLLFlag  && abs(outputFll(n))<10  %��Ƶ������״̬�£��ź��뱾��Ƶ��С��10ʱ
            loopCount = loopCount + 1;
            if  loopCount>20            
                   settings.PLLFlag = 1;
            end
       elseif  1 == settings.PLLFlag && abs(outputFll(n))>30      %�����໷����״̬�£���Ƶ�����������ź��뱾��Ƶ�����30ʱ
            loopCount = loopCount-1;
            if  0 == loopCount
                settings.PLLFlag = 0;
            end
       end
    end
 %�뻷������
end

trackResult2.carrPhase = carrStartPhaseSum / startCountPhase;

if  startCountPhase > 0
    trackResult2.trackFlag = 1;
end