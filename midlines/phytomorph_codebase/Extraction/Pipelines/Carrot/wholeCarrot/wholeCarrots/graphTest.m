clear all
clear classes
close all
%%
close all
X = randi([-50 50],[1000 1]);
Y = randi([-50 50],[1000 1]);

G = phG([X Y]);
% take top N points from table and add as nodes
N = 100;
for e = 1:N
    n = phN(e);
    G.addNode(n);
end
G.displayNodes();
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%
% add random single edges
% take random M from N
M = 20;
for e = 1:M
    i1 = randi(N,1);
    i2 = randi(N,1);
    ed = phE([i1 i2]);
    G.addEdge(ed);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%
% add edges of length pl
[G] = generatePaths(G,3,[15],N);
G.displayEdges();
G.displayEndPoints();

fidx = findnEndPoints(G);

startPoint = sparse(size(X,1),1);
startPoint(fidx(1)) = 1;
[startPoint] = G.traceToEnd(startPoint,5);
str = find(startPoint(:,1));
stp = find(startPoint(:,end));
figure;
hold all
G.displayEdges();
G.displayNodes('g*',str);
G.displayNodes('k*',stp);



