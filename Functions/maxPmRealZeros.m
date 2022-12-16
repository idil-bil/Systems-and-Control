function [dynamics, zeros, max_PM] = maxPmRealZeros(crossoverMagnitude, Kref, poleDynamics, openLoop)
    s = tf('s');
    z1_res = -crossoverMagnitude;
    z2_res = -crossoverMagnitude;
%     test_Dz = (1/(z1_res * z2_res)) * ((s-z1_res)*(s-z2_res));
%     [~, max_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
    max_PM = intmin;
    resolution = 0.005;

    for iter = 1:100
        while 1
            z1_test = z1_res + resolution;
            z2_test = z2_res - resolution;
            if z1_test >= 0
                break
            end
            test_Dz = (1/(z1_test * z2_test)) * ((s-z1_test)*(s-z2_test));
            [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
            if test_PM > max_PM
                max_PM = test_PM;
                z1_res = z1_test;
                z2_res = z2_test;
                D_res = test_Dz * poleDynamics;
            else
                break
            end
        end

        while 1
            z1_test = z1_res - resolution;
            z2_test = z2_res - resolution;
            test_Dz = (1/(z1_test * z2_test)) * ((s-z1_test)*(s-z2_test));
            [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
            if test_PM > max_PM
                max_PM = test_PM;
                z1_res = z1_test;
                z2_res = z2_test;
                D_res = test_Dz * poleDynamics;
            else
                break
            end
        end
    
        while 1
            z1_test = z1_res - resolution;
            z2_test = z2_res + resolution;
            if (z1_test) < (z2_test)
                break
            end
            test_Dz = (1/(z1_test * z2_test)) * ((s-z1_test)*(s-z2_test));
            [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
            if test_PM > max_PM
                max_PM = test_PM;
                z1_res = z1_test;
                z2_res = z2_test;
                D_res = test_Dz * poleDynamics;
            else
                break
            end
        end
    
        while 1
            z1_test = z1_res + resolution;
            z2_test = z2_res + resolution;
            if z1_test >= 0
                break
            end
            test_Dz = (1/(z1_test * z2_test)) * ((s-z1_test)*(s-z2_test));
            [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
            if test_PM > max_PM
                max_PM = test_PM;
                z1_res = z1_test;
                z2_res = z2_test;
                D_res = test_Dz * poleDynamics;
            else
                break
            end
        end
    end
    dynamics = D_res;
    zeros = [z1_res, z2_res];
end