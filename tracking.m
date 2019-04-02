function trackResult = tracking(fid, acqResult, settings,data)
% Performs code and carrier tracking for all channels.
%�����������ź�ʵ�� ����ٺ��ز�����
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record. �źż�¼���ļ���ʶ��
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRun.m from
%                       acquisition results).
%                       ��Ҫ���ٵ��������ǵ��ز�Ƶ�ʺ�����λ����preRun.m���ݲ������ṩ��
%       settings        - receiver settings. ���ջ�������
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.
%       ���ٽ�� - ���ٽ�����ṹ�����飩
%                             ������ in-phase prompt outputs
%                                         ������չ�����ʼλ��
%                                         �Ӹ��ٻ��õ��Ĺ۲�����

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================
%   ��ʼ������ṹ��
% Channel status    ͨ��״̬
trackResult.status         = '-';      % No tracked signal, or lost lock   �޸����źŻ�ʧ��

% The absolute sample in the record of the C/A code start:  C/A����ʼ��¼�ľ��Բ���
trackResult.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the C/A code: C/A���Ƶ��
trackResult.codeFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave: ���ٵ��ز�Ƶ��
trackResult.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):  I��
trackResult.I_P            = zeros(1, settings.msToProcess);
trackResult.I_E            = zeros(1, settings.msToProcess);
trackResult.I_L            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase): Q��
trackResult.Q_E            = zeros(1, settings.msToProcess);
trackResult.Q_P            = zeros(1, settings.msToProcess);
trackResult.Q_L            = zeros(1, settings.msToProcess);

% Loop discriminators ��·����
trackResult.dllDiscr       = inf(1, settings.msToProcess);
trackResult.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResult.pllDiscr       = inf(1, settings.msToProcess);
trackResult.pllDiscrFilt   = inf(1, settings.msToProcess);

% ÿһ���������ز���C/A��ĳ�ʼ��λ
trackResult.carrStartPhase      = inf(1, settings.msToProcess);
trackResult.codeStartPhase      = inf(1, settings.msToProcess);

%% Initialize tracking variables ==========================================
codePeriods = settings.msToProcess;

K_step = settings.K_step;  %FLLת��PLLʱ��

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;   %dll��ؼ��

% Summation interval �ܵ�ʱ����,��������ʱ��
PDIcode = 0.001;
 
% Calculate filter coefficient values   �����˲���ϵ��ֵ
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);
%FLL������ʼ��
past_I_accum = 0;
past_Q_accum = 0;
step_cnt = 0;

hwb = waitbar(0,'Tracking...');

%% Start processing channels ==============================================
% Get a vector with the C/A code sampled 1x/chip
caCode = cacode(settings.PRN);
% Then make it possible to do early and late versions
caCode = [caCode(1023) caCode caCode(1)];
%--- Perform various initializations ------------------------------

% define initial code frequency basis of NCO
codeFreq      = settings.codeFreqBasis;
codeFreqBasis = settings.codeFreqBasis;
% define residual code phase (in chips) ����ʣ�����׶�
remCodePhase  = settings.codeLength * ( 1 - acqResult.codePhase/settings.samplesPerCode);
% define carrier frequency which is used over whole tracking period
carrFreq      = acqResult.carrFreq;
carrFreqBasis = acqResult.carrFreq;
% ����ʣ���ز���λ,ָ����ÿһ����Ƭ֮�󣬽��������ز���λ�Ƕ��١�
remCarrPhase  = 0.0;        

%code tracking loop parameters
oldCodeNco   = 0.0;
oldCodeError = 0.0;

%carrier/Costas loop parameters
oldCarrNco   = 0.0;
oldCarrError = 0.0;

%=== Process the number of specified code periods =================
for loopCnt =  1:codePeriods

%% GUI update -------------------------------------------------------------
    if (rem(loopCnt, 50) == 0)
        try
            waitbar(loopCnt/codePeriods, ...
                    hwb, ...
                    ['Completed ',int2str(loopCnt), ...
                      ' of ', int2str(codePeriods), ' msec']);                       
        catch
            % The progress bar was closed. It is used as a signal
            % to stop, "cancel" processing. Exit.
            disp('Progress bar closed, exiting...');
            return
        end
    end
    
    
%% Read next block of data ------------------------------------------------            
    %�����������в��ҡ��顱��������ڵĴ�С�����ݴ���Ƶ�ʣ��ɱ䣩�Ͳ���Ƶ�ʣ��̶���������λ���ԡ�
    trackResult.codeStartPhase(loopCnt) = remCodePhase;
    trackResult.carrStartPhase(loopCnt) = remCarrPhase;
    
    codePhaseStep = codeFreq / settings.samplingFreq;

    blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);

    % Read in the appropriate number of samples to process this
    % interation 
    rawSignal = data(fid:fid + blksize - 1);
    fid = fid + blksize ;

%% Set up all the code phase tracking information -------------------------
    % Define index into early code vector
    %����������Ϊ���ڴ�������
    tcode       = (remCodePhase-earlyLateSpc) : ...
                  codePhaseStep : ...
                  ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
    tcode2      = ceil(tcode) + 1;
    earlyCode   = caCode(tcode2);

    % Define index into late code vector
    tcode       = (remCodePhase+earlyLateSpc) : ...
                  codePhaseStep : ...
                  ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
    tcode2      = ceil(tcode) + 1;
    lateCode    = caCode(tcode2);

    % Define index into prompt code vector
    tcode       = remCodePhase : ...
                  codePhaseStep : ...
                  ((blksize-1)*codePhaseStep+remCodePhase);
    tcode2      = ceil(tcode) + 1;
    promptCode  = caCode(tcode2);

    remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;

%% Generate the carrier frequency to mix the signal to baseband -----------
    time    = (0:blksize) ./ settings.samplingFreq;

    % Get the argument to sin/cos functions
    trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
    remCarrPhase = rem(trigarg(blksize+1), (2 * pi));

    % Finally compute the signal to mix the collected data to bandband
    % �������źţ����ɼ�����������Ƶ����ϡ�
    carrCos = cos(trigarg(1:blksize));
    carrSin = sin(trigarg(1:blksize));

%% Generate the six standard accumulated values ---------------------------
    % First mix to baseband
    qBasebandSignal = carrCos .* rawSignal;
    iBasebandSignal = carrSin .* rawSignal;

    % Now get early, late, and prompt values for each
    I_E = sum(earlyCode  .* iBasebandSignal);
    Q_E = sum(earlyCode  .* qBasebandSignal);
    I_P = sum(promptCode .* iBasebandSignal);
    Q_P = sum(promptCode .* qBasebandSignal);
    I_L = sum(lateCode   .* iBasebandSignal);
    Q_L = sum(lateCode   .* qBasebandSignal);

%% Find FLL error or  PLL error and update carrier NCO ----------------------------------
    if step_cnt<=K_step 
        %FLL������
        dot = past_I_accum * I_P   +   past_Q_accum * Q_P;
        cross = past_I_accum * Q_P  -    past_Q_accum * I_P;

        %�����޷����м���
        carrError = atan2(cross,dot) / (2*pi*PDIcode);
        past_I_accum = I_P;
        past_Q_accum = Q_P;   
        carrNco = oldCarrNco + carrError;

    else
        %PLL������
        carrError = atan(Q_P / I_P) / (2.0 * pi);
        carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
        (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);

    end 

    % Implement carrier loop filter and generate NCO command ִ���ز���·�˲���������NCO����


    oldCarrNco   = carrNco;
    oldCarrError = carrError;
    if step_cnt == K_step 
        oldCarrError = atan(past_Q_accum / past_I_accum) / (2.0 * pi);
    end
    step_cnt = step_cnt + 1;
    %CarrNcos = [CarrNcos  carrNco];
    % Modify carrier freq based on NCO command ����NCO�����޸��ز�Ƶ��
    carrFreq = carrFreqBasis + carrNco;
    trackResult.carrFreq(loopCnt) = carrFreq;

%% Find DLL error and update code NCO -------------------------------------
    codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
        (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

    % Implement code loop filter and generate NCO command
    codeNco = oldCodeNco + (tau2code/tau1code) * ...
        (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
    oldCodeNco   = codeNco;
    oldCodeError = codeError;
    %codeNcos = [codeNcos codeNco];
    % Modify code freq based on NCO command
    codeFreq = codeFreqBasis - codeNco;

    trackResult.codeFreq(loopCnt) = codeFreq;

%% Record various measures to show in postprocessing ----------------------
    % Record sample number (based on 8bit samples)��¼��Ʒ��ţ�����8λ��Ʒ��
    trackResult.absoluteSample(loopCnt) = fid - blksize;

    trackResult.dllDiscr(loopCnt)       = codeError;
    trackResult.dllDiscrFilt(loopCnt)   = codeNco;
    trackResult.pllDiscr(loopCnt)       = carrError;
    trackResult.pllDiscrFilt(loopCnt)   = carrNco;

    trackResult.I_E(loopCnt) = I_E;
    trackResult.I_P(loopCnt) = I_P;
    trackResult.I_L(loopCnt) = I_L;
    trackResult.Q_E(loopCnt) = Q_E;
    trackResult.Q_P(loopCnt) = Q_P;
    trackResult.Q_L(loopCnt) = Q_L;

end % for loopCnt
close(hwb)
