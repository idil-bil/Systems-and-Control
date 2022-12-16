%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: 
% maxPmComplexZeros
% Usage:
% Finds the complex conjugate zeros that result in the
% maximum phase margin for a given set of 
% cross-over frequency, Kref, Dp, OL Xfer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dynamics, zeros, maxPm, valid, invalid_matrix] = maxPmComplexZeros(crossoverFreq, Kref, poleDynamicsXfer, openLoopXfer)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VARIABLE INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = tf('s');
    % flags to move in certain directions
    axis_sel = 1; % 1: real is set, jw changes ; 0: jw is set, real changes
    vertical = 1; % 1: go upwards; 0: go downwards
    horizontal = 1; % 1:
    % store relevant Pms
    maxPm = 0;
    prevPm = 0;
    checkPm1 = 0;
    checkPm2 = 0;
    % store zeros
    z1_test = intmin;
    z2_test = intmin;
    z1_res = intmin;
    z2_res = intmin;
    % store the Dynamics xfer
    D_res = s;
    % adjust the resolution
    resolution = 0.0005; 
    % start at the cross-over freq
    sigma = -crossoverFreq; 
    % start at real-axis
    omega = resolution;
    % valid flag to check the selected point at the end
    valid = 1;
    % invalid matrix
    invalid_matrix = [0 0 0; 0 0 0; 0 0 0];
    % counting variable to move diagonally when gets too high
    noMoveCnt = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ITERATION BLOCK / NEWTON'S METHOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iteration = 1:1000 % total number of direction traversing iterations
        while 1
            % two constraints
            if(sigma + resolution >= 0)
                break;
            end
            if(omega - resolution < 0)
                break;
            end
            % zero declarations
            z1_test = sigma + omega*1i;
            z2_test = sigma - omega*1i;
            % Dz computation
            test_Dz = getDzPID(z1_test, z2_test);
            % Pm computation
            [~, testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);

            % if local max is found, i.e. Pm starts dropping, take action
            if ((testPm) < (prevPm))

                % increment noMoveCnt by 1 each time this happens
                noMoveCnt = noMoveCnt + 1;

                % restore the maxPm point sigma and omega
                if (axis_sel == 1)
                    if(vertical == 1)
                        omega = omega - resolution;
                    else
                        omega = omega + resolution;
                    end
                else
                    if(horizontal == 1)
                        sigma = sigma - resolution;
                    else
                        sigma = sigma + resolution;
                    end
                end

                % change axis of direction
                if (axis_sel == 1)
                    axis_sel = 0;
                else
                    axis_sel = 1;
                end

                % check adjacent points to decide which direction is the best to move in
                if (axis_sel == 1)
                    z1_test = sigma + (omega+resolution)*1i;
                    z2_test = sigma - (omega+resolution)*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, checkPm1] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
                    z1_test = sigma + (omega-resolution)*1i;
                    z2_test = sigma - (omega-resolution)*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, checkPm2] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
                    if(checkPm1>checkPm2)
                        vertical = 1;
                    else 
                        vertical = 0;
                    end
                else
                    z1_test = (sigma+resolution) + omega*1i;
                    z2_test = (sigma+resolution) - omega*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, checkPm1] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
                    z1_test = (sigma-resolution) + omega*1i;
                    z2_test = (sigma-resolution) - omega*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, checkPm2] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
                    if(checkPm1>checkPm2)
                        horizontal = 1;
                    else 
                        horizontal = 0;
                    end
                end
                break; % finish moving in one direction, iterate in another direction
            else  % update the values with higher Pm point
                maxPm = testPm;
                z1_res = z1_test;
                z2_res = z2_test;
                D_res = test_Dz * poleDynamicsXfer;
                prevPm = testPm;
                noMoveCnt = 0;
            end

            % depending on the direction flags, update the omega and sigma points for a new set of zeros
            if (axis_sel == 1)
                if(vertical == 1)
                    omega = omega + resolution;
                else
                    omega = omega - resolution;
                end
            else
                if(horizontal == 1)
                    sigma = sigma + resolution;
                else
                    sigma = sigma - resolution;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RETURN VALUES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dynamics = D_res;
    zeros = [z1_res, z2_res];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAX PM VALIDITY CHECK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z1 = z1_res-resolution;
    z2 = z2_res-resolution;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [0 0 0; 1 0 0; 0 0 0];
    end
    z1 = z1_res+resolution;
    z2 = z2_res+resolution;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [0 0 0; 0 0 1; 0 0 0];
    end
    z1 = double(z1_res)-resolution*1i;
    z2 = double(z2_res)+resolution*1i;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [0 0 0; 0 0 0; 0 1 0];
    end
    z1 = double(z1_res)+resolution*1i;
    z2 = double(z2_res)-resolution*1i;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [0 1 0; 0 0 0; 0 0 0];
    end
    z1 = double(z1_res)-resolution-resolution*1i;
    z2 = double(z2_res)-resolution+resolution*1i;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [0 0 0; 0 0 0; 1 0 0];
    end
    z1 = double(z1_res)-resolution+resolution*1i;
    z2 = double(z2_res)-resolution-resolution*1i;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [1 0 0; 0 0 0; 0 0 0];
    end
    z1 = double(z1_res)+resolution-resolution*1i;
    z2 = double(z2_res)+resolution+resolution*1i;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [0 0 0; 0 0 0; 0 0 1];
    end
    z1 = double(z1_res)+resolution+resolution*1i;
    z2 = double(z2_res)+resolution-resolution*1i;
    test_Dz = getDzPID(z1,z2);
    [~,testPm] = margin(Kref * poleDynamicsXfer * test_Dz * openLoopXfer);
    if(maxPm < testPm) 
        valid = 0; 
        invalid_matrix = invalid_matrix + [0 0 1; 0 0 0; 0 0 0];
        z1_res = z1;
        z2_res = z2;
    end
end
