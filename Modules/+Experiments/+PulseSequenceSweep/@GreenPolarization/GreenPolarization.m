classdef GreenPolarization < Experiments.PulseSequenceSweep.PulseSequenceSweep_invisible
    %GreenPolarization measures the time dependence of the PLE signal

    properties(SetObservable,AbortSet)
        greenLaser = Modules.Source.empty(1,0); % Allow selection of source
        APDline = 3;
        delayTime = 10;
        offset = 1;
        countDur = 0.02;
        onTime = 100;
    end
    
    properties
        placeHolderVariable = 1; %all APD bins are acquired in one shot, no variable is swept
        counterTimes = 0;
    end
    
    properties(Constant)
        nCounterBins = 20; %number of APD bins for this pulse sequence (with more than 20 the PB errors)
        vars = {'placeHolderVariable'}; %names of variables to be swept
    end
    methods(Static)
        obj = instance()
    end
    methods(Access=private)
        function obj = GreenPolarization()
            obj.prefs = [obj.prefs,{'APDline','delayTime','offset',...
            'countDur','onTime','greenLaser'}]; %additional preferences not in superclass
            obj.loadPrefs;
        end
    end

    methods
        pulseSeq = BuildPulseSequence(obj,Ind) %Defined in separate file
        
        function PreRun(obj,~,~,ax)
            %prepare axes for plotting
            hold(ax,'on');
            %plot data bin 1
            plotH = plot(ax,obj.counterTimes,squeeze(obj.data.sumCounts(1,:)),'color','b');
            ylabel(ax,'Counts');
            xlabel(ax,'counter time (\mus)');
            hold(ax,'off');
            set(ax,'xlimmode','auto','ylimmode','auto','ytickmode','auto')
        end
        
        function UpdateRun(obj,~,~,ax,~,~)
            if obj.averages > 1
                averagedData = squeeze(nanmean(obj.data.sumCounts,1));
                
            else
                averagedData = squeeze(obj.data.sumCounts);
               
            end
            
             plot(ax,linspace(0,obj.onTime + obj.offset,obj.nCounterBins) - obj.offset,averagedData,'color','b')
             drawnow;
        end
        
        function set.onTime(obj,val)
            %Note that order matters here; setting tauTimes first is
            %important in case of error
            obj.onTime = val;
            obj.counterTimes = linspace(0,obj.onTime + obj.offset,obj.nCounterBins);
        end
    end
end