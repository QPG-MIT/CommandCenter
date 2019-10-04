classdef GreenCalib < Experiments.PulseSequenceSweep.PulseSequenceSweep_invisible
    %GreenCalib Calibrates the powers and pulse durations for 532 NV
    %polarization and state readout

    properties(SetObservable,AbortSet)
        exciteLaser = Modules.Source.empty(1,0); % Allow selection of source
        APDline = 3;
        repumpTime_us = 'linspace(0.1,1,101)'; 
        countTime_us = 'linspace(0.1,0.5,41)';
        delayTime = 100;
    end
    
    properties
        repumpTime = linspace(0.1,1,101);
        countTime = linspace(0.1,0.5,41); %will be in us
    end
    
    properties(Constant)
        nCounterBins = 2; %number of APD bins for this pulse sequence
        vars = {'repumpTime','countTime'}; %names of variables to be swept
    end
    methods(Static)
        obj = instance()
    end
    methods(Access=private)
        function obj = GreenCalib()
            obj.prefs = [obj.prefs,{'exciteLaser','APDline','repumpTime','countTime','delayTime'}]; %additional preferences not in superclass
            obj.loadPrefs;
        end
    end

    methods
        pulseSeq = BuildPulseSequence(obj,repInd,countInd) %Defined in separate file
        
        function PreRun(obj,~,~,ax)
            %prepare axes for plotting
            hold(ax,'on');
            %plot data bin 1
            plotH = plot(ax,obj.tauTimes,obj.data.sumCounts(:,1,1),'color','b');
            %plot data bin 1 errors
            plotH(2) = plot(ax,obj.tauTimes,obj.data.sumCounts(:,1,1)+obj.data.stdCounts(:,1,1),'color',[1 .5 0],'LineStyle','--'); %upper bound
            plotH(3) = plot(ax,obj.tauTimes,obj.data.sumCounts(:,1,1)-obj.data.stdCounts(:,1,1),'color',[1 .5 0],'LineStyle','--'); %lower bound
            %plot data bin 2
            plotH(4) = plot(ax,obj.tauTimes,obj.data.sumCounts(:,2,1),'color','b');
            %plot data bin 2 errors
            plotH(5) = plot(ax,obj.tauTimes,obj.data.sumCounts(:,2,1)+obj.data.stdCounts(:,2,1),'color',[1 .5 0],'LineStyle','--'); %upper bound
            plotH(6) = plot(ax,obj.tauTimes,obj.data.sumCounts(:,2,1)-obj.data.stdCounts(:,2,1),'color',[1 .5 0],'LineStyle','--'); %lower bound
            ax.UserData.plots = plotH;
            ylabel(ax,'Normalized PL');
            xlabel(ax,'Delay Time \tau (\mus)');
            hold(ax,'off');
            set(ax,'xlimmode','auto','ylimmode','auto','ytickmode','auto')
        end
        
        function UpdateRun(obj,~,~,ax,~,~)
            if obj.averages > 1
                averagedData = squeeze(nanmean(obj.data.sumCounts,3));
                meanError = squeeze(nanmean(obj.data.stdCounts,3));
            else
                averagedData = obj.data.sumCounts;
                meanError = obj.data.stdCounts;
            end
            
            %grab handles to data from axes plotted in PreRun
            ax.UserData.plots(1).YData = averagedData(:,1);
            ax.UserData.plots(2).YData = averagedData(:,1) + meanError(:,1);
            ax.UserData.plots(3).YData = averagedData(:,1) - meanError(:,1);
            ax.UserData.plots(4).YData = averagedData(:,2);
            ax.UserData.plots(5).YData = averagedData(:,2) + meanError(:,2);
            ax.UserData.plots(6).YData = averagedData(:,2) - meanError(:,2);
            drawnow;
        end
        
        function set.countTime_us(obj,val)
            %Note that order matters here; setting tauTimes first is
            %important in case of error
            tempvals = eval(val);
            obj.tauTimes = tempvals;
            obj.tauTimes_us = val;
        end
        
        function set.repumpTime_us(obj,val)
            %Note that order matters here; setting tauTimes first is
            %important in case of error
            tempvals = eval(val);
            obj.tauTimes = tempvals;
            obj.tauTimes_us = val;
        end
    end
end
