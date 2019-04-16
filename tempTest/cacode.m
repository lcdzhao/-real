function CAcode = cacode(PRN,settings)
% generateCAcode.m generates one of the 32 GPS satellite C/A codes.
%
% CAcode = generateCAcode(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       CAcode      - a vector containing the desired C/A code sequence 
%                   (chips).  

%--- Make the code shift array. The shift depends on the PRN number -------
% The g2s vector holds the appropriate shift of the g2 code to generate
% the C/A code (ex. for SV#19 - use a G2 shift of g2s(19) = 471)
g2s = [ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,...
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

%--- Pick right shift for the given PRN number ----------------------------
g2shift = g2s(PRN);

%--- Generate G1 code -----------------------------------------------------

%--- Initialize g1 output to speed up the function ---
g1 = zeros(1, settings.codeLength);
%--- Load shift register ---
reg = -1*ones(1, 10);

%--- Generate all G1 signal chips based on the G1 feedback polynomial -----
for i=1:settings.codeLength
    g1(i)       = reg(10);
    saveBit     = reg(3)*reg(10);
    reg(2:10)   = reg(1:9);
    reg(1)      = saveBit;
end

%--- Generate G2 code -----------------------------------------------------

%--- Initialize g2 output to speed up the function ---
g2 = zeros(1, settings.codeLength);
%--- Load shift register ---
reg = -1*ones(1, 10);

%--- Generate all G2 signal chips based on the G2 feedback polynomial -----
for i=1:settings.codeLength
    g2(i)       = reg(10);
    saveBit     = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
    reg(2:10)   = reg(1:9);
    reg(1)      = saveBit;
end

%--- Shift G2 code --------------------------------------------------------
%The idea: g2 = concatenate[ g2_right_part, g2_left_part ];
g2 = [g2(settings.codeLength-g2shift+1 : settings.codeLength), g2(1 : settings.codeLength-g2shift)];

%--- Form single sample C/A code by multiplying G1 and G2 -----------------
CAcode = -(g1 .* g2);
