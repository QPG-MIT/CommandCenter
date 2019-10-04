function s = BuildPulseSequence(obj,repInd,countInd)
%BuildPulseSequence Builds pulse sequence for performing all-optical T1
%characterization given the index (tauIndex) in tauTimes

s = sequence('GreenCalib');
exciteChannel = channel('repump','color','g','hardware',1);
% exciteChannel = channel('repump','color','g','hardware',obj.exciteLaser.PBline-1);
APDchannel = channel('APDgate','color','b','hardware',obj.APDline,'counter','APD1');
s.channelOrder = [exciteChannel, APDchannel];
tb_rep = node(s.StartNode,exciteChannel,'delta',0);
te_rep = node(tb_rep,exciteChannel,'units','us','delta',obj.repumpTime(repInd)+obj.countTime(countInd));

node(te_rep, APDchannel,'units','us','delta',-obj.countTime(countInd));
te_count = node(te_rep, APDchannel, 'units','us', 'delta', 0);

tb_count = node(te_count,exciteChannel,'units','us','delta',obj.delayTime);
node(tb_count, APDchannel,'units','us','delta',0);
te_count = node(tb_count, exciteChannel,'units','us', 'delta',obj.countTime(countInd));
node(te_count, APDchannel, 'units','us', 'delta', 0);

end

