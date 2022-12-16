%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: 
% getDzPID
% Usage:
% Finds the part of the dynamics of a PID controller
% that is associated with the zeros chosen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dynamics = getDzPID(z1_neg, z2_neg)
    s = tf('s');
    dynamics = (1/(z1_neg * z2_neg)) * ((s-z1_neg)*(s-z2_neg));
end