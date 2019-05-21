function [ detectedCycleSlip ] = CycleSlipDetector_VelocityTrend( iGNSS, L1C_R, D1C_R, satPRNL1, satLOSL1, L2C_R, D2C_R, satPRNL2, satLOSL2, interval_R, cycleSlipDetectionThreshold )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if iGNSS>2
    % Cycle Slip Detector
    recordCycleSlipL1_ourRTK = [];
    if length(satLOSL1)>2
        currentPhiL1_R = L1C_R(satLOSL1, iGNSS);
        previousDopplerL1_R = -D1C_R(satLOSL1, iGNSS-1);
        currentDopplerL1_R = -D1C_R(satLOSL1, iGNSS-1);
        previousPhiL1_R = L1C_R(satLOSL1, iGNSS-1);
        currentPhiL1_B = L1C_B(satLOSL1, iGNSS);
        currentDopplerL1_B = -D1C_B(satLOSL1, iGNSS-1);
        previousDopplerL1_B = -D1C_B(satLOSL1, iGNSS-1);
        previousPhiL1_B = L1C_B(satLOSL1, iGNSS-1);
        velocityTrendTestL1 = [];
        estimatedPhiL1_R = [];
        for iSat = 1:length(satLOSL1)
            estimatedPhiL1_R(iSat) = previousPhiL1_R(iSat) + 1/2 * interval_R * ( previousDopplerL1_R(iSat) + currentDopplerL1_R(iSat) );
            velocityTrendTestL1(iSat) = abs( estimatedPhiL1_R(iSat) - currentPhiL1_R(iSat) );
            if velocityTrendTestL1(iSat) > cycleSlipDetectionThreshold
                detectedCycleSlip = [ detectedCycleSlip; satPRNL1(iSat) ];
                recordCycleSlipL1_ourRTK(iGNSS,satPRNL1(iSat)) = 0;
            end
        end
    end
    if length(satLOSL2)>2
        detectedCycleSlip = [];
        currentPhiL2_R = L2C_R(satLOSL2, iGNSS);
        previousDopplerL2_R = -D2C_R(satLOSL2, iGNSS-1);
        currentDopplerL2_R = -D2C_R(satLOSL2, iGNSS-1);
        previousPhiL2_R = L2C_R(satLOSL2, iGNSS-1);
        currentPhiL2_B = L2C_B(satLOSL2, iGNSS);
        currentDopplerL2_B = -D2C_B(satLOSL2, iGNSS-1);
        previousDopplerL2_B = -D2C_B(satLOSL2, iGNSS-1);
        previousPhiL2_B = L2C_B(satLOSL2, iGNSS-1);
        velocityTrendTestL2 = [];
        estimatedPhiL2_R = [];
        for iSat = 1:length(satLOSL2)
            estimatedPhiL2_R(iSat) = previousPhiL2_R(iSat) + 1/2 * interval_R * ( previousDopplerL2_R(iSat) + currentDopplerL2_R(iSat) );
            velocityTrendTestL2(iSat) = abs( estimatedPhiL2_R(iSat) - currentPhiL2_R(iSat) );
            if velocityTrendTestL2(iSat) > cycleSlipDetectionThreshold
                detectedCycleSlip = [ detectedCycleSlip; satPRNL2(iSat) + 100 ];
                recordCycleSlipL2(iGNSS,satPRNL2(iSat)) = 0;
            end
        end
    end
    
end

end

