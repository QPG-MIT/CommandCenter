times = 10:2:20;

for i=1:length(times)
    exp.resTime_us = times(i);
    managers.Experiment.run;
    if managers.Experiment.aborted
        return
    end
end
