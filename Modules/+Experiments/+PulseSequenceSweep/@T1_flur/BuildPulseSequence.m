function s = BuildPulseSequence(obj,tauIndex)
%BuildPulseSequence Builds pulse sequence for performing all-optical T1
%characterization given the index (tauIndex) in tauTimes

s = sequence('T1_flur');
exciteChannel = channel('repump','color','g','hardware',1);
% exciteChannel = channel('repump','color','g','hardware',obj.exciteLaser.PBline-1);
APDchannel = channel('APDgate','color','b','hardware',obj.APDline,'counter','APD1');
s.channelOrder = [exciteChannel, APDchannel];

t0 = node(s.StartNode,exciteChannel,'delta',0);
node(t0, APDchannel, 'delta', 0);
t_RepEnd = node(t0,exciteChannel,'units','us','delta',obj.repumpTime_us);
node(t_RepEnd, APDchannel, 'delta', 0);

t_Read = node(t_RepEnd, exciteChannel, 'units', 'us', 'delta', obj.tauTimes(tauIndex));
node(t_Read,APDchannel,'delta',0);
t_ReadEnd = node(t_Read ,exciteChannel,'units','us','delta',obj.countTime);
node(t_ReadEnd,APDchannel,'delta',0);

end

