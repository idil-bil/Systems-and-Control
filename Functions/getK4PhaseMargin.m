%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: 
% getK4PhaseMargin
% Usage:
% Finds the master gain K needed for 
% a desired Phase Margin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = getK4PhaseMargin(Kref, fullDynamics, openLoop, targetPhaseMargin)
    K = Kref;
    while 1
        [~, testPM] = margin(K * fullDynamics * openLoop);
        if abs(testPM) <= targetPhaseMargin
            break
        end
        K = K + 0.01;
    end
end