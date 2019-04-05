function trackResult = tracking(fid, acqResult, settings,data)
%% test
loop_para = loop_canshu_calculate(settings);%���㻷·�˲�������

code_table = cacode(settings.PRN);          %���ú���������α�����

fll_nco_adder = 0;                   %fll  NCO�ӷ�����Ӧ�����˲�����ʹ��
carrier_nco_sum = 0;                 %�����ֵ�ٳ�2*pi������λ��
pll_nco_adder = 0;                   %pll  NCO�ӷ�����Ӧ�����˲�����ʹ��
loop_count = 0;                      %ѭ������
code_nco_sum = 0;                    %�����һ������λ����Ķ���
code_nco_adder = 0;                  %ddll  NCO�ӷ�����Ӧ�����˲�����ʹ��
n_IQ = 2;                            
n = 3;
output_fll(2:3) = 0;                 %����������fll��Ƶ�������
output_filter_fll(1:3) = 0;          %fll ��·�˲��������
output_filter_pll(1:3) = 0;          %pll ��·�˲��������
output_pll(2:3) = 0;                 %pll ��Ƶ�������
output_filter_ddll(1:3) = 0;         %ddll ��·�˲��������
pll_after_filter = 0;

Tcoh = settings.Tcoh;          %��������ʱ��
global early_code_nco;      %�ǽ��ն˵ģ���Ϊ����Ļ�����ʱ������붼���Դ����ж���
early_code_nco = ((1 - (acqResult.codePhase-3)/settings.samplesPerCode)...
    *settings.codeLength) * 2^settings.nco_Length;
early_code_nco = mod(early_code_nco,2^settings.nco_Length*1023);
local_early_code_last = local_earlycode_initial(settings,code_table); %�������س�ǰ�룬���ն�ʹ�ã���Ϊ����Ļ�����ʱ������붼���Դ����ж���

trackResult.carrier_nco_phase = zeros(1,settings.msToProcess);       %ÿ�λ�������ǰ�ز���nco��λ
trackResult.code_nco_phase = zeros(1,settings.msToProcess);          %ÿ�λ�������ǰB1���nco��λ
trackResult.carrierFreq = zeros(1,settings.msToProcess);
blksize = settings.samplesPerCode;

for loop_num = 1 : settings.msToProcess
    carrier_nco_phase = mod(carrier_nco_sum/2^settings.nco_Length,1) * 2 * pi;  
    if carrier_nco_phase > pi
        trackResult.carrier_nco_phase(loop_num) = ...
            2*pi - carrier_nco_phase;
    else
        trackResult.carrier_nco_phase(loop_num) = carrier_nco_phase;
    end
    trackResult.code_nco_phase(loop_num) = ...
        ((((early_code_nco/2^settings.nco_Length)/settings.codeLength)*settings.samplesPerCode...
         -3) /settings.samplesPerCode)*settings.codeLength;       
    trackResult.carrierFreq(loop_num) = ...
         (settings.middle_freq_nco1 + fll_nco_adder + pll_nco_adder)/settings.transfer_coef;
    trackResult.flag(loop_num) = settings.PLL_flag;         %��ʶ�ô�ѭ����û�н���PLL����
    %fd_plot(loop_num) = settings.dup_freq;      %��¼�ôεĶ�����Ƶ�Ƶ���ʵֵ
    


    % Read in the appropriate number of samples to process this
    % interation 
    receive_signal = data(fid:fid + blksize - 1);
    fid = fid + blksize ;

          
    %code_nco_sum = code_nco_adder + settings.code_word + fll_nco_adder*(1/763);  %���������뻷NCO     
  
    
    
    
    
    
    
    
    
    %[signal_modulate_code,settings.signal_phase] = signalcode(settings,code_table);  %�����������CB1�룬��ʼ��λΪsettings.modulate_code_bias_phsae



    
%     receive_signal =  original_signal;  %��Ƶ����ź�
    %�������������ز�
    for demond_num = 1:settings.Ncoh 
        local_cos(demond_num) = cos(2*pi*carrier_nco_sum/2^settings.nco_Length);
        local_sin(demond_num) = -sin(2*pi*carrier_nco_sum/2^settings.nco_Length);
        carrier_nco_sum = carrier_nco_sum + settings.middle_freq_nco1 + fll_nco_adder + pll_nco_adder ;%���������ز�NCO
        
        %carrier_nco_sum�����һ������λ�й�ϵ�Ķ���
    end
    
    code_nco_sum = code_nco_adder + settings.code_word ...              %���������뻷NCO,�������carrier_nco_sum����һ��,2048Ϊ��׼
        + fll_nco_adder*settings.cofe_FLL_auxi_DDLL;                
  
      %�������س�ǰ����ʱ���ͺ���
    [local_early_code,local_prompt_code,local_late_code,settings.local_phase]=localcode_generate(local_early_code_last,code_nco_sum,code_table,settings);
    local_early_code_last = local_early_code;
    %�ز����    
    I_demon_carrier = local_cos.*receive_signal;
    Q_demon_carrier = local_sin.*receive_signal;
%     save_I_demon_carrier = [save_I_demon_carrier I_demon_carrier];
%     save_Q_demon_carrier = [save_Q_demon_carrier Q_demon_carrier];
    
    %�źŽ������������
    I_E_final = sum(I_demon_carrier.*local_early_code);
    Q_E_final = sum(Q_demon_carrier.*local_early_code);
    I_P_final(n_IQ) = sum(I_demon_carrier.*local_prompt_code);
    Q_P_final(n_IQ) = sum(Q_demon_carrier.*local_prompt_code);
    I_L_final = sum(I_demon_carrier.*local_late_code);
    Q_L_final = sum(Q_demon_carrier.*local_late_code);
    
    
%     I_P_final(n_IQ) = sum(I_demon_carrier);
%     Q_P_final(n_IQ) = sum(Q_demon_carrier);
    
    
    if  1 == loop_num
        I_P_final(n_IQ - 1) = I_P_final(n_IQ);
        Q_P_final(n_IQ - 1) = Q_P_final(n_IQ);
    else
% %         �����޷����м�Ƶ��
        dot_fll = I_P_final(n_IQ - 1) * I_P_final(n_IQ) + Q_P_final(n_IQ - 1) * Q_P_final(n_IQ);
        cross_fll = I_P_final(n_IQ - 1) * Q_P_final(n_IQ) - I_P_final(n_IQ) * Q_P_final(n_IQ - 1);
        output_fll(n) = atan2(cross_fll,dot_fll)/(Tcoh*2*pi); 
        result_discriminator_Fll(loop_num) = output_fll(n);
        
        output_filter_fll(n) = (loop_para.cofeone_FLL * output_fll(n)) + (loop_para.cofetwo_FLL * output_fll(n - 1)) + (2 * output_filter_fll(n - 1)) - output_filter_fll(n - 2);
        fll_after_filter(loop_num) = output_filter_fll(n);
        
        fll_nco_adder = output_filter_fll(n) * settings.transfer_coef ;  %Ƶ����ת��      
        output_fll(n - 1)=output_fll(n);
        output_filter_fll(n - 2)=output_filter_fll(n - 1);
        output_filter_fll(n - 1)=output_filter_fll(n);
        
         if settings.PLL_flag == 1
            %���໷������
            output_pll(n) = atan2(Q_P_final(n_IQ),I_P_final(n_IQ)); 
            output_filter_pll(n) = loop_para.cofeone_PLL*output_pll(n) + loop_para.cofetwo_PLL*output_pll(n-1)+loop_para.cofethree_PLL*output_pll(n-2)+2*output_filter_pll(n-1)-output_filter_pll(n-2);
            result_discriminator_Pll(loop_num) = output_pll(n);
            pll_after_filter(loop_num) = output_filter_pll(n);
            pll_nco_adder = (output_filter_pll(n)/(2*pi)) * settings.transfer_coef;  %Ƶ����ת��
            
%             output_pll(1:2) = output_pll(2:3);
%             output_filter_pll(1:2) = output_filter_pll(2:3);
            output_pll(n-2) = output_pll(n-1);
            output_pll(n-1) = output_pll(n);
            output_filter_pll(n-2) = output_filter_pll(n-1);
            output_filter_pll(n-1) = output_filter_pll(n);
         end
        
        I_P_final(n_IQ - 1) = I_P_final(n_IQ);
        Q_P_final(n_IQ - 1) = Q_P_final(n_IQ);
       if 0 == settings.PLL_flag  && abs(output_fll(n))<10  %��Ƶ������״̬�£��ź��뱾��Ƶ��С��10ʱ
            loop_count = loop_count + 1;
            if  loop_count>200            
                   settings.PLL_flag = 1;
            end
       elseif  1 == settings.PLL_flag && abs(output_fll(n))>30      %�����໷����״̬�£���Ƶ�����������ź��뱾��Ƶ�����30ʱ
            loop_count = loop_count-1;
            if  0 == loop_count
                settings.PLL_flag = 0;
            end
       end
    end
 %�뻷������
    output_ddll(n) = ((I_E_final - I_L_final)*I_P_final(n_IQ) + (Q_E_final - Q_L_final)*Q_P_final(n_IQ) )/((I_P_final(n_IQ)^2 + Q_P_final(n_IQ)^2)*2);  % DDLL_discri_1
       
    result_ddll(loop_num) = output_ddll(n);
    %�뻷�˲��������ף�
    output_filter_ddll(n) = output_filter_ddll(n -1) + (loop_para.cofeone_DDLL*output_ddll(n)) + loop_para.cofetwo_DDLL*output_ddll(n - 1);
    result_DDLL_filter(loop_num) = output_filter_ddll(n);
    % ת����Ƶ�ʿ�����
    code_nco_adder = output_filter_ddll(n) * settings.transfer_coef ; %Ƶ����ת��
%     Code_NCO=0;
%     C(loop_num)=code_nco_adder;
    %�滻
    output_ddll(n - 1)=output_ddll(n);
    output_filter_ddll(n - 1) = output_filter_ddll(n);
    code_phase_discrim(loop_num) = settings.signal_phase - settings.local_phase ;
end

