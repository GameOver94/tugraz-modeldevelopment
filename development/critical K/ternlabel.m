% TERNLABEL label ternary phase diagram
%   TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') labels a ternary phase diagram created using TERNPLOT
%   
%   H = TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') returns handles to the text objects created.
%   with the labels provided.  TeX escape codes are accepted.
%
%   See also TERNPLOT

% Author: Carl Sandrock 20020827

% To Do

% Modifications

% Modifiers



function h = ternlabel(A, B, C)

    r(1) = text(0.5, -0.05, A, 'horizontalalignment', 'center','Color','blue','FontSize',12);
    r(2) = text(1-0.4*sin(deg2rad(30)), 0.5, B,'rotation', -60,'horizontalalignment', 'center','Color','red','FontSize',12);
    r(3) = text(0.4*sin(deg2rad(30)), 0.5, C, 'rotation', 60,'horizontalalignment', 'center','Color','green','FontSize',12);

if nargout > 0
    h = r;
end