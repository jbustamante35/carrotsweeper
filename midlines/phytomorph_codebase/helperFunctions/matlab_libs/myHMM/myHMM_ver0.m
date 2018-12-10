clear classes
%%
n1 = hmm_node('state1');
n2 = hmm_node('state2');
n1.attachNode(n2,.5);
n1.attachNode(n1,.5);
n2.attachNode(n1,.7);
n2.attachNode(n2,.3);
nn = n1.getNextNode();
n1_d = myProb(10,5);
n2_d = myProb(0,5);
n1.attachDistribution(n1_d);
n2.attachDistribution(n2_d);
O = [];
for e = 1:10000
    nn = nn.getNextNode();
    O(e,:) = nn.observe();
end
hist(O,100)
%%
n1 = hmm_node('state1');
n2 = hmm_node('state2');
n3 = hmm_node('state3');
n1.attachNode(n2,.5);
n1.attachNode(n1,.5);
n2.attachNode(n1,.7);
n2.attachNode(n2,.2);
n2.attachNode(n3,.1);
n3.attachNode(n3,.9);
n3.attachNode(n2,.1);
nn = n1.getNextNode();
n1_d = myProb(10,5);
n2_d = myProb(0,5);
n3_d = myProb(30,5);
n1.attachDistribution(n1_d);
n2.attachDistribution(n2_d);
n3.attachDistribution(n3_d);
hmm = my_hmm();
hmm.addNode(n1);
hmm.addNode(n2);
hmm.addNode(n3);
hmm.state = 1;
hmm.dn = .1;
O = hmm.getObservations(10);
P = hmm.forwardCalc(O,[.5 .25 .25]);
Po = hmm.forwardCalc(-10*ones(1,10),[.5 .25 .25]);
tm = hmm.buildTransitionMatrix();
%s = hmm.Viterbi(O);
s = hmm.Viterbi(30*ones(1,10));
%s = hmm.Viterbi([0 10 0 30])
