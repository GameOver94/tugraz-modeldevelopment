% Calculation binary interaction parameter for NRTL

% 1 - Hexane
% 2 - Methanol
% 3 - Cloroform

function tij = tau(T,i)

switch i
    case 13 %parameter 13
        A=0;%NRTL(2,j);
        B=-154.502;%NRTL(4,j); %3407.1*4.2/8.314;
        E=0;
        F=0;
        
    case 31 %parameter 31
        A=0;%NRTL(3,j);
        B=360.657;%NRTL(5,j); %1662*4.2/8.314;
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