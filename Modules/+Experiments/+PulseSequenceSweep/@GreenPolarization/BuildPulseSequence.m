function s = BuildPulseSequence(obj,~)
%BuildPulseSequence Builds pulse sequence for performing Optical Spin
%Polarization measurements

counterTimes = linspace(0,obj.onTime + obj.offset,obj.nCounterBins);

s = sequence('GreenPolarization');
laserChannel = channel('repump','color','g','hardware',obj.greenLaser.PBline-1);
APDchannel = channel('APDgate','color','b','hardware',obj.APDline-1,'counter','APD1');
s.channelOrder = [laserChannel, APDchannel];

count_s = node(s.StartNode,laserChannel,'units','us','delta', obj.delayTime);
node(count_s,laserChannel,'units','us','delta',obj.onTime);

for ii = 1:length(counterTimes)
    counterStart = node(count_s,APDchannel,'units','us','delta',counterTimes(ii) - obj.offset);
    node(counterStart,APDchannel,'units','us','delta',obj.countDur);
end

end