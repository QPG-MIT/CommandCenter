function stitchedClosedDAQ(managers)
%     GG = 10.^(-2:.5:1.5);
% %     GG = 10.^(1.5:1:1.5);
%     
%     C = Sources.Cobolt_PB.instance;
%     
%     C.arm()
%     C.on()
    dwl = .01;
    WL = 406.550:dwl:406.850;
    
    S = Sources.msquared.SolsTiS.instance;
    
    for wl = WL
        wl
        S.WavelengthLock(true);
        S.TuneSetpoint(wl);
        pause(.5)
        S.GetPercent
        
        while abs(S.GetPercent - 50) > 5
            S.TuneSetpoint(wl + 2*dwl);
            S.TuneSetpoint(wl + (rand-.5)*dwl/40);
            pause(.5)
            S.GetPercent
        end

        managers.Experiment.run()
    end
% g
%     E = Experiments.WidefieldSlowScan.Closed.instance;
%     E.repump_always_on = false;
%     
%     for g = GG(end:-1:1)
%         C.CW_power = g;
% 
%         managers.Experiment.run()
%     end
end