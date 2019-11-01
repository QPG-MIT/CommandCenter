classdef GreenCalib < Experiments.PulseSequenceSweep.PulseSequenceSweep_invisible
    %GreenCalib Calibrates the powers and pulse durations for 532 NV
    %polarization and state readout

    properties(SetObservable,AbortSet)
        exciteLaser = Modules.Source.empty(1,0); % Allow selection of source
        APDline = 3;
        repumpTime = 1; 
        countTime_us = 'linspace(0.02,0.2,10)';
        delayTime = 1;
    end
    properties
        countTime = linspace(0.1,0.5,41); %will be in us
    end
    
    properties(Constant)
        nCounterBins = 2; %number of APD bins for this pulse sequence
        vars = {'countTime'}; %names of variables to be swept
    end
    methods(Static)
        obj = instance()
    end
    methods(Access=private)
        function obj = GreenCalib()
            obj.prefs = [obj.prefs,{'exciteLaser','APDline','repumpTime','countTime_us','delayTime'}]; %additional preferences not in superclass
            obj.loadPrefs;
        end
    end

    methods
        pulseSeq = BuildPulseSequence(obj,repInd,countInd) %Defined in separate file
        
        function PreRun(obj,~,~,ax)
            %prepare axes for plotting
            hold(ax,'on');
            plotH = plot(ax,obj.countTime,squeeze(obj.data.sumCounts(1,:,1)),'color','b');
            ylabel(ax,'Normalized PL');
            xlabel(ax,'Delay Time \tau (\mus)');
            hold(ax,'off');
            set(ax,'xlimmode','auto','ylimmode','auto','ytickmode','auto')
        end
        
        function UpdateRun(obj,~,~,ax,~,~)
            if obj.averages > 1
                averagedData = squeeze(nanmean(obj.data.sumCounts,1));
                meanError = squeeze(sqrt(nanmean(obj.data.stdCounts.^2,1)));
            else
                averagedData = squeeze(obj.data.sumCounts);
                meanError = squeeze(sqrt(nanmean(obj.data.stdCounts.^2,1)));
            end
            
            %grab handles to data from axes plotted in PreRun
            ax.UserData.plots(1).YData = averagedData(:,1);
            plot(ax,obj.countTime,(averagedData(:,2)./averagedData(:,1)),'color','b')
            drawnow;
        end
        
        function set.countTime_us(obj,val)
            %Note that order matters here; setting tauTimes first is
            %important in case of error
            tempvals = eval(val);
            obj.countTime = tempvals;
            obj.countTime_us = val;
        end
        
%         function set.repumpTime_us(obj,val)
%             %Note that order matters here; setting tauTimes first is
%             %important in case of error
%             tempvals = eval(val);
%             obj.repumpTime = tempvals;
%             obj.repumpTime_us = val;
%         end
    end
end
