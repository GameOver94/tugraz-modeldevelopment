%% Function tau
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  Estimate liquid-liquid equilibrium 
%           of ternary system - n-hexane(1)/methanol(2)/chloroform(3)
%           or Reactive LLE - of A+B->C

% Calculation binary interaction parameter for NRTL

function tij = tau(T,i)

% 1 - n- hexane
% 2 - methanol
% 3 - chloroform

switch i
    case 13 %parameter 13
        A=0;
        B=-154.502;
        E=0;
        F=0;        
    case 31 %parameter 31
        A=0;
        B=360.657;
        E=0;
        F=0;        
    case 12 %parameter 21
        A=-3.6551;
        B=1507.15;
        E=0;
        F=0;       
    case 21 %parameter 12
        A=-1.1544;
        B=734.514;
        E=0;
        F=0;        
    case 23 %parameter 23
        A=0;
        B= -71.9029;
        E=0;
        F=0;        
    case 32 %parameter 32
        A=0;
        B= 690.066;
        E=0;
        F=0;
    otherwise
        error('no valid argument')
end

tij = A+B./T+E.*log(T)+F.*T;
end