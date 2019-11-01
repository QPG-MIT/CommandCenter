function s = BuildPulseSequence(obj,countInd)
%BuildPulseSequence Builds pulse sequence for performing all-optical T1
%characterization given the index (tauIndex) in tauTimes

s = sequence('GreenCalib');
%exciteChannel = channel('repump','color','g','hardware',1);
exciteChannel = channel('repump','color','g','hardware',obj.exciteLaser.PBline-1);
APDchannel = channel('APDgate','color','b','hardware',obj.APDline-1,'counter','APD1');
s.channelOrder = [exciteChannel, APDchannel];
tb_rep = node(s.StartNode,exciteChannel,'units','us','delta',0);
te_rep = node(tb_rep,exciteChannel,'units','us','delta',obj.repumpTime);

ts_count = node(te_rep,exciteChannel,'units','us','delta',10);
te_count = node(ts_count,exciteChannel,'units','us','delta',obj.repumpTime);
node(ts_count, APDchannel,'units','us','delta',0);
node(ts_count, APDchannel, 'units','us', 'delta', obj.countTime(countInd));

tb_count = node(te_count,exciteChannel,'units','us','delta',obj.delayTime);
node(tb_count, exciteChannel,'units','us', 'delta',obj.repumpTime);
node(tb_count, APDchannel,'units','us','delta',0);
node(tb_count, APDchannel, 'units','us', 'delta', obj.countTime(countInd));

end

