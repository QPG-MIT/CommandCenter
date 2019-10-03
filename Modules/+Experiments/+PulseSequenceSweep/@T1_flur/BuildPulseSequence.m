function s = BuildPulseSequence(obj,tauIndex)
%BuildPulseSequence Builds pulse sequence for performing all-optical T1
%characterization given the index (tauIndex) in tauTimes

s = sequence('T1_flur');
exciteChannel = channel('repump','color','g','hardware',1);
% exciteChannel = channel('repump','color','g','hardware',obj.exciteLaser.PBline-1);
APDchannel = channel('APDgate','color','b','hardware',obj.APDline,'counter','APD1');
s.channelOrder = [exciteChannel, APDchannel];
g = node(s.StartNode,exciteChannel,'delta',0);
g = node(g,exciteChannel,'units','us','delta',obj.repumpTime_us);

r = node(g,exciteChannel,'units','us','delta',obj.tauTimes(tauIndex));
node(r,APDchannel,'delta',0);
r = node(r,exciteChannel,'units','us','delta',obj.countTime);
node(r,APDchannel,'delta',0);

end

