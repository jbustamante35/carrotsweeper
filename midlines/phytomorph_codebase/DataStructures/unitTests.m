%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit testing for tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%
    % constructor tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % create a basic tensor with co-low signature
    t1 = tensor(rand(3,4,5));
    % create a basic tensor with high signature
    t2 = tensor(rand(3,1,2));
    % create a basic tensor with high signature
    t3 = tensor(rand(3,1));
    t4 = tensor(rand(1,3));
    % create a tensor with defined signature
    t5 = tensor(rand(3,4,5));

    %%%%%%%%%%%%%%%%%%%%%%%%
    % dot product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % test the dot product along the dims
    r = grt.dotProduct(t1,t3,[[1 1];[1 1]]);
    r = grt.dotProduct(t1,t2,[[1 1];[1 1]]);
    r = grt.dotProduct(t1,t5,[[1 1];[3 3]]);
    
    r = grt.apply([t1;t1;t1],[t5;t5;t5],@plus);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % map test with plus
    %%%%%%%%%%%%%%%%%%%%%%%%
    op = @(A,B)plus;
    
    [r] =  map(objA,objB,d,op)
    
    
    f1 = tensor([tensor([t1;t1;t1]),tensor([t1;t1;t1]),tensor([t1;t1;t1])]);
    f2 = tensor([tensor([t1;t1;t1]),tensor([t1;t1;t1]),tensor([t1;t1;t1])]);
    f3 = tensor([tensor([t1;t1;t1]),tensor([t1;t1;t1]),tensor([t1;t1;t1])]);
    fT = [[f1;f2;f3],[f1;f2;f3]];
    
    g1 = tensor([tensor([t1;t1;t1]),tensor([t1;t1;t1]),tensor([t1;t1;t1])]);
    g2 = tensor([tensor([t1;t1;t1]),tensor([t1;t1;t1]),tensor([t1;t1;t1])]);
    g3 = tensor([tensor([t1;t1;t1]),tensor([t1;t1;t1]),tensor([t1;t1;t1])]);
    gT = [[g1;g2;g3],[g1;g2;g3]];
    
    
    r = grt.apply(fT,gT,@plus);
    
    
    r = grt.dotProduct([t1;t1;t1],[t5;t5;t5],[[1 1];[3 3]]);
%%
%% TEST:1
clear all
syms f g y1 y2
syms a b real


A = rand(2,2);
B = rand(2,2);

A = rand(2,3,5);
B = rand(3,2,6);


along = [1 1];
along = [[1 2];[2 1]];

AS = sym(['A'],size(A));
BS = sym(['B'],size(B));
%t(1).d = [2 5];
%t(2).d = [4 2 5];
t(1).d = [2];
t(2).d = [2];

t(1).d = [2];
t(2).d = [2 3];

t(1).d = [2 3];
t(2).d = [2 3];

t(1).d = [2 3 4];
t(2).d = [2 3];

t(1).d = [3 2 4];
t(2).d = [2 3];

t(1).d = size(A);
t(2).d = size(B);


d = [t.d];
L = [1 cumprod(d)];
iia = 1:prod(t(1).d);
iib = 1:prod(t(2).d);

iib = iib + iia(end);
basis = [];
mIN = 1:prod(d);
IN =  sym(['IN'],[1 prod(d)],'real');
initC = 0;
cnt = 1;
for ten = 1:numel(t)
    % spin up basis for the ten-th tensor
    for order = 1:numel(t(ten).d)
        % assign basis for the o-th order
        basis(ten).f{order} = sym(['t' num2str(ten) num2str(order) '_'],[1 t(ten).d(order)],'real');
        basis(ten).g{order} = sym(['t' num2str(ten) num2str(order) '_'],[1 t(ten).d(order)],'real');
        basis(ten).e{order} = rand(1,t(ten).d(order));
        iVec = (1:t(ten).d(order)) + initC;
        basis(ten).m{order} = IN(iVec);
        basis(ten).n{order} = iVec;
        initC = iVec(end);
        if order == 1
            lN = 1;
        else
            lN = cumprod(t(ten).d(1:order-1));
        end
        basis(ten).new{order} = (t(ten).d(order)^-1)*lN;
        cnt= 1;
    end
end
%% TEST:2 dot(plus) product along 
%along = [1 1];
%along = [[1 1];[2 2]];
%along = [[2 1];[1 2]];
for e = 1:size(along,1)
    % numeric values
    basis(1).e{along(e,1)} = basis(2).e{along(e,2)};
    % symbolic values
    basis(1).f{along(e,1)} = basis(2).f{along(e,2)};
end
%% TEST:3 accounting for result size
ra = ones(size(t(1).d));
rb = ones(size(t(2).d));
rr = zeros(1,size(along,1));
for e = 1:size(along,1)
    ra(along(e,1)) = 0;
    rb(along(e,2)) = 0;
    if t(1).d(along(e,1)) == t(2).d(along(e,2))
        rr(e) = t(1).d(along(e,1));
    else
        fprintf(['Dims to not add up for instruction./n']);
    end
end

ra = ra.*t(1).d;
rb = rb.*t(2).d;
ra(ra==0) = [];
rb(rb==0) = [];
rt = [ra rr rb];

pa = 0:numel(ra);
pb = 0:numel(rb);
pr = 0:numel(rr);

% funny integration
pr = pr + pa(end);
pb = pb + pr(end);



pr(1) = [];
pa(1) = [];
pb(1) = [];
pt = [pr pa pb];
%% TEST:4 law
x = sym(['x'],[1 numel(d)]);
%f(x(:)) = x(1) == x(2);
%f(x(:)) = x(1) == x(4);
%f(x(:)) = x(1) == x(2);
%f(x(:)) = x(1) == x(5) | x(2) == x(4);
%f(x(:)) = isequaln(x(1),x(2));
x0 = numel(t(1).d);
IDX = nchoosek(1:numel(d),2);
e = 1;
f(x(:)) = kroneckerDelta(x(IDX(e,1)) - x(IDX(e,2)));
for e = 2:size(IDX,1)
    %f(x(:)) = x(t(1).d(along(e,1))) == x(x0 + t(2).d(along(e,2)));
    f(x(:)) = [f , kroneckerDelta(x(IDX(e,1)) - x(IDX(e,2)))];
    %f(x(:)) = f*kroneckerDelta(x(IDX(e,1)) - x(IDX(e,2)));
end
Z = zeros(1,size(IDX,1));
for e = 1:size(along,1)
    query = along(e,:);
    query(2) = query(2) + x0;
    fidx = find(all(bsxfun(@eq,IDX,query),2) | all(bsxfun(@eq,IDX,fliplr(query)),2));
    %fidx = find(all(bsxfun(@eq,IDX,query),2));
    Z(fidx) = 1;
end
f = prod(f.^Z);
q = f;
%% TEST:5 subs on f
IN_withPairing(x(:)) = prod(x);
IN_with_OUT_Pairing(x(:)) = prod(x);
F = f;
slot = numel(x);
% due to fill order
tenIndex = 1:numel(basis);
tenIndex = flipdim(tenIndex,2);
for ten = tenIndex
    orderIndex = flipdim(1:numel(basis(ten).e),2);
    for order = orderIndex
        f = subs(f,x(slot),basis(ten).e{order}(:));
        F = subs(F,x(slot),basis(ten).m{order}(:));
        IN_withPairing = subs(IN_withPairing,x(slot),basis(ten).f{order}(:));
        IN_with_OUT_Pairing = subs(IN_with_OUT_Pairing,x(slot),basis(ten).g{order}(:));
        
        slot = slot - 1;
    end
end
g = formula(IN_with_OUT_Pairing);
h = formula(IN_withPairing);
simplify(f)
sum(simplify(f))
numel(formula(f))
%% TEST:6 add in operator
op = @(x,y)plus(x,y);
f = formula(f);
r = zeros(numel(f),1);
vA = A(:);
vB = B(:);
cnt = 1;
syms rs;
r = [];

%IN = [AS(:);BS(:)];
for bi = 1:numel(vB)
    for ai = 1:numel(vA)
        if f(cnt)
            rs(cnt) = op(IN(iia(ai)),IN(iib(bi)));
            rs2(cnt) = op(AS(ai),BS(bi));
        end
        cnt = cnt + 1;
    end
end


tic
cnt = 1;
for bi = 1:numel(vB)
    for ai = 1:numel(vA)
        if f(cnt)
            r(cnt) = op(A(ai),B(bi));
        end
        cnt = cnt + 1;
    end
end
toc


rs(r==0) = [];

r(r==0) = [];
r = reshape(r,rt);
%rs = reshape(rs,rt);
func = matlabFunction(rs);
func = func2str(func);
clip = strfind(func,')');
func(1:clip(1)) = [];
for e = fliplr(1:numel(IN))
    toR = ['X(' num2str(e) ')'];
    func = strrep(func,char(IN(e)),toR);
end
func = str2func(['@(X)' func]);

inV = [A(:);B(:)];
tic;here = func(inV);toc
tic;what = inV + inV;toc
inV = inV';
tic;hey = t0(inV,inV);toc
%%
x = rand(1,6600);
y = rand(1,6600);
tic;z1 = t0(x,y);toc
tic;z2 = x+y;toc
tic;z3 = t0(inV,inV);toc
%%
rs(IN(:)) = rs;
NIN = [A(:);B(:)];

func = matlabFunction(rs);
func(NIN(:)')
%rs(FS(:)) = permute(rs,pt);
%%

%%
for e = 1:numel(d)
    v = sym(['e' num2str(e) '_'],[1 d(e)],'real');
    f = subs(f,x(e), g(1:d(e)));
    f = subs(f,x(e),v);
end
%% remember to keep around the idea of the result tensor
for e = 1:numel(d)
    v = sym(['e' num2str(e) '_'],[1 d(e)],'real');
    f = subs(f,x(e), g(1:d(e)));
    f = subs(f,x(e),v);
end
%%
clear all
a = sym(['a'],[3 1]);    
b = sym(['b'],[1 4]);
c = a*b;
%% 
clear all

d = [2 3 4];
r = 1;
for e = 1:numel(d)
    v = sym(['e' num2str(e)],[1 d(e)],'real');
    %v = sym([num2str(e)],[1 d(e)]);
    r = r(:)*v;
end
%%
clear all
syms f g y1 y2
syms a b real
a = [2 3];
b = [2];
sv = {'a' 'b'};
d = [a b];
x = sym(['x'],[1 numel(d)]);



%f(x(:)) = sin(sum(x));
%f(x(:)) = kroneckerDelta(x(1)-x(2));
%f(x(:)) = x(1) == x(4) | x(2) == x(5) | x(3) == x(6);
%f(x(:)) = x(1) == x(3) | x(2) == x(4);
f(x(:)) = x(1) == x(2);


g = f;
%f(x(:)) = g(x(:));
r = 1;

for e = 1:numel(d)
    
    
    v = sym(['e' num2str(e) '_'],[1 d(e)],'real');
    
    func = @(a,b)a*b';
    
    tmp(y1,y2) = kroneckerDelta(y1-y2);
    %tmp(y1) = (subs(tmp,y2,1:d(e)).*v)*ones(d(e),1);
    tmp(y1) = (subs(tmp,y2,1:d(e)));
    g = func(tmp,v);
    %tmp(1:d(e))
    %f = subs(f,x(e), g(1:d(e)));
    %f = subs(f,x(e),tmp);
    
    %?
    %f = subs(f,x(e),v);
    
    
    %v = sym(['e'],[1 d(e)],'real');
    
    
    
    %v = sym(['e'],[1 d(e)]);
    %{
    % try to scub off component therefore use repmat
    v = sym(['e' num2str(e)]);
    v = repmat(v,[1 d(e)]);
    %}
    % best so far
    f = subs(f,x(e),v);
    %r = r(:)*v;
    %r = [r;v(:)];
end
%%
clear all
syms f x e g
e(x) = sin(x);
v = formula(subs(e,x,1:4))
%array = sym(['e'],[1 4],'real');
%%
syms f x y
a = sym(['e'],[1 3],'real');
f(x,y) = kroneckerDelta(x-y);
f(y) = (subs(f,x,1:3).*a)*ones(3,1);
%%
clear all
syms f x y a
a = sym(['e'],[3 1],'real');
f(x,y) = kroneckerDelta(x-y);
f(y,a) = (subs(f,x,1:numel(a)).*a');
%%
clear all
syms f x y a
a = sym(['e'],[3 1],'real');
v = sym(['e'],[3 1],'real');
f(a,v) = a'*v;
%%
f([y(:);a(:)]) = (subs(f,x,1:numel(a)).*a);
%% whats this, whats this
clear all
syms f x y a
v = sym(['e'],[3 1],'real');
f(x,y) = kroneckerDelta(x-y);
%f(y) = subs(f,x,1:3);
f = subs(f,x,1:3);
func = @(a,b)a*b;
func(formula(f),v)
%%
w = sym(['e'],[1 3],'real');
f(1,w)
%%
b = [v(1);v(2)];
%%
v = sym(['e'],[3 1],'real');
x = sym(['x'],[3 1],'real');
f(x) = eye(3)*x
%%
a = sym(['a' num2str(e)],[1 2],'real');
b = sym(['b' num2str(e)],[1 2],'real');
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    % tensor product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    s = grt.tensorProduct(t1,t3);
    q = grt.tensorProduct(t1,t2);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % permute tests  note that the signature needs to change too
    %%%%%%%%%%%%%%%%%%%%%%%%
    q.permute([1 2 4 3 5]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % folding tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    fT = tensor(rand(3,4,5,6,7,8));
    fT.fold({3 [4 5] [2 1 6]});
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % plus test
    %%%%%%%%%%%%%%%%%%%%%%%%
    a = tensor(rand(3,4,5,6,7,8));
    b = tensor(rand(3,4,5,6,7,8));
    c = a + b;
    d = a + 1;
    e = 1 + b;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q0 = tensor(rand(3,4,5));
    q0.fold({[1:q0.order()]});
    q0 = tensor(rand(3,4,5));
    1<q0
    q1 = tensor(rand(3,4,5));
    q1>1

    q0 = tensor(rand(3,4,5));
    q1 = tensor(rand(3,4,5));
    r = 1<q0|q1>1;

    q0 = tensor(rand(60,1));
    q1 = tensor(rand(60,1));
    Q = tensor(rand(60,60));
    r = 1<q0|Q|q1>1;



    q0 = zeros(11,1);
    q0(3) = 1;
    q0 = tensor(q0);
    q1 = zeros(11,1);
    q1(3) = 1;
    q1 = tensor(q1);
    Q = tensor(diag(-5:5));
    r = 1<q0|Q|q1>1;
    
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['************************************************\n']);
fprintf(['*********** Start:unit test: class - tE *********\n']);
fprintf(['*** Purpose: label a type co/vector along a dim.\n']);
fprintf(['************************************************\n']);
%%%%%%%%%%%%%%%%%%%%%%
% constructor testing
%%%%%%%%%%%%%%%%%%%%%%
t1 = tE([1],10);
t2 = tE([0]);
t3 = tE([1],3);
t4 = tE([1 0 1],20);
t5 = tE([1 1],5);
t6 = tE([1 0],2);
T = [t1 t2 t3];
t2 == t3;
t4.toString()
t4 = t5.directSum(t6);
t6 = T.fold();
t8.toString()
fprintf(['************************************************\n']);
fprintf(['*********** End:unit test: class - tE *********\n']);
fprintf(['************************************************\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create 
t1 = tE([1]);
t2 = tE([1]);
t3 = tE([1]);

t4 = tE([0]);
t5 = tE([0]);



T1 = TE([t1 t4]);
T2 = TE([t2 t4]);
T3 = TE([t3 t4]);

T4 = TE([t1 t5]);
T5 = TE([t2 t5]);
T6 = TE([t3 t5]);

T = [[T1;T2;T3],[T4;T5;T6]];

%T.permute([2 1]);
%T = T.fold({[1 2]});



q1 = [T1;T2;T3];


q1 = TE([1]);
q2 = TE([1]);
qT = [q1 q2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% smooth image
I = tensor(imread('/home/nate/Downloads/20170329n02_04.tif'));
interpFunc = @(P,X)ba_interp2(X{1},P{1},P{2});
mI = smoothImage(interpFunc,str2func('smoothImage'),I);
X = [];
[X(:,:,2),X(:,:,1)] = ndgrid(linspace(1,300,300),linspace(1,300,300));
squareDomain = smoothDomain('','',X);
sI = mI.evaluate(squareDomain);
sI.plot()
mI.plot();
%%

% make the inner part of the grid
[POL(:,:,1),POL(:,:,2)] = ndgrid(linspace(0,6,10),linspace(-pi,pi,100));

EE = [[1,0];[0,1]];
clear x
x = sym(['I' '_%d'],[2 1],'real');
T = symbolicTensor(EE,[x(1)*cos(x(2));x(1)*sin(x(2))]);
Y = T(POL)
%% test smooth domain and symbolic tensor
clear p
p = sym(['p' '_%d'],[2 1],'real');
ROT = [[cos(p(1));sin(p(1))],[-sin(p(1));cos(p(1))]];
R = symbolicTensor(ROT);
clear POL
[POL(:,:,1),POL(:,:,2)] = ndgrid(linspace(0,6,10),linspace(-pi,pi,100));
testDomain = smoothDomain('','',POL);
paraDomain = smoothDomain('','',pi/8);
Y = R(testDomain,paraDomain);

%% test symbolic tensor compose
clear M rho theta
syms M iM rho theta x y real

M(rho,theta) = [rho*cos(theta);rho*sin(theta)];
iM(x,y) = [(x.^2 + y.^2).^-.5;atan2(y,x)];
T = symbolicTensor(M,[1 2],iM);
[POL(:,:,1),POL(:,:,2)] = ndgrid(linspace(0,6,10),linspace(-pi,pi,100));
polarDomain = geometry2D('','',POL);
Y = T(polarDomain);
Y.plot();
%% play with polar cart
syms alpha M real
x = sym(['x'],[2 1],'real');
M(alpha) = [[cos(alpha);sin(alpha)],[-sin(alpha);cos(alpha)]];
M(0)
r(x) = x'*M(0)*x;
%%
W = sym(['W'],[2 2 2 2],'real');
S = symbolicTensor(W,[2 2 2 2]);
S.tryStuff()
%%
clear classes
syms alpha M
M(alpha) = [[cos(alpha);sin(alpha)],[-sin(alpha);cos(alpha)]];
x = sym(['x'],[2 1],'real');
y = sym(['y'],[1 2],'real');
s([alpha(:);x(:);y(:)]) = y*M*x
t([alpha(:);x(:);y(:)]) = finverse(s,x(1));
Y([alpha(:);x(:)]) = M*x;
yn = Y(pi/8,4,5);
t(alpha,x(1),x(2),y(1),y(2));
simplify(t(pi/8,1,0,yn(1),yn(2)))
simplify(t(pi/8,4,0,yn(1),yn(2)))

t([alpha(:);x(:);y(:)]) = finverse(s,alpha);
Y([alpha(:);x(:)]) = M*x;
yn = Y(pi/8,4,5);
l=1
ra = vpa(simplify(t(1,4,5,yn(1),yn(2))))
rb = vpa(simplify(t(ra^-1,4,5,yn(1),yn(2))))
rc = vpa(simplify(t(pi/8,4,5,yn(1),yn(2))))
vpa(s(pi/8,4,5,yn(1),yn(2)))

solve(s(alpha,4,5,yn(1),yn(2)) == 41)
%%
syms x y
test1(x,y) = y/x;
test2 = finverse(test1,x);
test1(test2(x,y),y)


%%
syms alpha real
ROT(alpha) = [[cos(alpha);sin(alpha)],[-sin(alpha);cos(alpha)]];
R = symbolicTensor(ROT,[2 2]);
R.apply(pi/8)
%%
syms x y real
Xpos(x,y) = [x;y];
X = symbolicTensor(Xpos,[2 2]);

%%



syms l1 l2 real
LAM(l1,l2) = [[l1;0],[0;l2]];
iLAM = inv(LAM);
L = symbolicTensor(LAM,[2 2],iLAM);


syms M iM rho theta x y real
CPmap(rho,theta) = [rho*cos(theta);rho*sin(theta)];
iCPmap(x,y) = [(x.^2 + y.^2).^-.5;atan2(y,x)];
CP = symbolicTensor(CPmap,[1 2],iCPmap');



[POL(:,:,1),POL(:,:,2)] = ndgrid(linspace(0,6,10),linspace(-pi,pi,100));
polarDomain = geometry2D('','',POL);
Y.plot();


[CAR(:,:,1),CAR(:,:,2)] = ndgrid(linspace(-20,20,40),linspace(-20,20,40));
cartDomain = geometry2D('','',CAR);



alpha = tensor(pi/4);
lV1 = tensor([2]);
lV2 = tensor([.5]);

close all
figure
TS1 = [R X];
Y = TS1.applyForward(cartDomain,alpha);
iY = TS1.applyInverse(Y,alpha);
Y.plot({'r' 'r'});
cartDomain.plot({'k' 'k'});
waitforbuttonpress
iY.plot({'b' 'b'});


figure;
TS1 = [R L CP];
Y = TS1.applyForward(polarDomain,alpha,lV1,lV2);
iY = TS1.applyInverse(Y,alpha,lV1,lV2);
Y.plot({'r' 'r'});
polarDomain.plot({'k' 'k'});
waitforbuttonpress
iY.plot({'b' 'b'});
%%
figure;
test1 = CP.applyForward(polarDomain);
test1.plot()
test2 = CP.applyInverse(test1);
test2.plot()

%%
LX.plot({'k' 'k'});
%LRX.plot({'m' 'm'});
r = TS(1);
for e = 2:numel(TS)
    r = r.dot(TS(e),[1 1]);
end
rp = r.applyForward(cartDomain,lV1,lV2,alpha);
rp.plot({'m' 'm'});
%%
[CAR(:,:,1),CAR(:,:,2)] = ndgrid(linspace(-20,20,40),linspace(-20,20,40));
cartDomain = geometry2D('','',CAR);
Y = R(cartDomain);
Y.plot();
%%
N = symbolicTensor.dotProduct(R,T,[2 1]);

%%
% make the inner part of the grid
[POL(:,:,1),POL(:,:,2)] = ndgrid(linspace(0,6,10),linspace(-pi,pi,100));
domain1X = tensorContainer(tensor(X1(:,:,1)));
domain1Y = tensorContainer(tensor(X1(:,:,2)));
D1 = [domain1X domain1Y];

% make the inner part of the grid
[X2(:,:,1),X2(:,:,2)] = ndgrid(linspace(6,12,10),linspace(-pi,pi,100));
domain1X = tensorContainer(tensor(X2(:,:,1)));
domain1Y = tensorContainer(tensor(X2(:,:,2)));
D2 = [domain1X domain1Y];

% make the inner part of the grid
[X3(:,:,1),X3(:,:,2)] = ndgrid(linspace(12,24,10),linspace(-pi,pi,100));
domain1X = tensorContainer(tensor(X3(:,:,1)));
domain1Y = tensorContainer(tensor(X3(:,:,2)));
D3 = [domain1X domain1Y];


para = tensorContainer(tensor([pi/8 1.5 1 100 100]));


R1 = cT.evaluate(D1,para);
R2 = cT.evaluate(D2,para);
R3 = cT.evaluate(D3,para);

figure;
hold on

mI.plot();
subI = mI.evaluate(R1);

R1.plot({'r' ,'r'})
R2.plot({'g' ,'g'})
R3.plot({'b' ,'b'})



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit testing for tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % constructor tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % create a basic tensor with co-low signature
    t1 = tensor(rand(3,4,5));
    % create a basic tensor with high signature
    t2 = tensor(rand(3,1,2));
    % create a basic tensor with high signature
    t3 = tensor(rand(3,1));
    t4 = tensor(rand(1,3));
    % create a tensor with defined signature
    t5 = tensor(rand(3,4,5));

    %%%%%%%%%%%%%%%%%%%%%%%%
    % dot product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % test the dot product along the dims
    r = tensor.dotProduct(t1,t3,[1 1]);
    r = tensor.dotProduct(t1,t2,[1 1]);
    r = tensor.dotProduct(t1,t5,[3 3]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % tensor product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    s = tensor.tensorProduct(t1,t3);
    q = tensor.tensorProduct(t1,t2);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % permute tests  note that the signature needs to change too
    %%%%%%%%%%%%%%%%%%%%%%%%
    q.permute([1 2 4 3 5]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % folding tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    fT = tensor(rand(3,4,5,6,7,8));
    fT.fold({3 [4 5] [2 1 6]});
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % plus test
    %%%%%%%%%%%%%%%%%%%%%%%%
    a = tensor(rand(3,4,5,6,7,8));
    b = tensor(rand(3,4,5,6,7,8));
    c = a + b;
    d = a + 1;
    e = 1 + b;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q0 = tensor(rand(3,4,5));
    q0.fold({[1:q0.order()]});
    q0 = tensor(rand(3,4,5));
    1<q0
    q1 = tensor(rand(3,4,5));
    q1>1

    q0 = tensor(rand(3,4,5));
    q1 = tensor(rand(3,4,5));
    r = 1<q0|q1>1;

    q0 = tensor(rand(60,1));
    q1 = tensor(rand(60,1));
    Q = tensor(rand(60,60));
    r = 1<q0|Q|q1>1;



    q0 = zeros(11,1);
    q0(3) = 1;
    q0 = tensor(q0);
    q1 = zeros(11,1);
    q1(3) = 1;
    q1 = tensor(q1);
    Q = tensor(diag(-5:5));
    r = 1<q0|Q|q1>1;
    

    %%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit testing for tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % constructor tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % create a basic tensor with co-low signature
    t1 = tensor(rand(3,4,5));
    % create a basic tensor with high signature
    t2 = tensor(rand(3,1,2),1);
    % create a basic tensor with high signature
    t3 = tensor(rand(3,1),1);
    t4 = tensor(rand(1,3),1);
    % create a tensor with defined signature
    t5 = tensor(rand(3,4,5),[0 1 1]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % dot product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % test the dot product along the dims
    r = tensor.dotProduct(t1,t3,[1 1]);
    r = tensor.dotProduct(t1,t2,[1 1]);
    r = tensor.dotProduct(t1,t5,[3 3]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % tensor product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    s = tensor.tensorProduct(t1,t3);
    q = tensor.tensorProduct(t1,t2);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % permute tests  note that the signature needs to change too
    %%%%%%%%%%%%%%%%%%%%%%%%
    q.permute([1 2 4 3 5]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % folding tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    fT = tensor(rand(3,4,5,6,7,8),[1 1 1 0 0 1]);
    fT.fold({3 [4 5] [2 1 6]});
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % plus test
    %%%%%%%%%%%%%%%%%%%%%%%%
    a = tensor(rand(3,4,5,6,7,8),[1 1 1 0 0 1]);
    b = tensor(rand(3,4,5,6,7,8),[1 1 1 0 0 1]);
    c = a + b;
    d = a + 1;
    e = 1 + b;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q0 = tensor(rand(3,4,5));
    q0.fold({[1:q0.order()]});
    q0 = tensor(rand(3,4,5));
    1<q0
    q1 = tensor(rand(3,4,5));
    q1>1

    q0 = tensor(rand(3,4,5),[1 1 1]);
    q1 = tensor(rand(3,4,5));
    r = 1<q0|q1>1;

    q0 = tensor(rand(60,1),[1]);
    q1 = tensor(rand(60,1),[0]);
    Q = tensor(rand(60,60),[0 1]);
    r = 1<q0|Q|q1>1;



    q0 = zeros(11,1);
    q0(3) = 1;
    q0 = tensor(q0,[1]);
    q1 = zeros(11,1);
    q1(3) = 1;
    q1 = tensor(q1,[0]);
    Q = tensor(diag(-5:5),[0 1]);
    r = 1<q0|Q|q1>1;
    %% symbolic
    
    
    
    syms d x0 x g
    d(x,x0) = dirac(x-x0);
    d = subs(d,x0,10)
    d = subs(d,x,10)
    d(x,x0) = x + x0;
    d = int(d,x0);
   
    
    d(x) = dirac(x);
    d = int(d,x,-inf,inf);
    
    d(x,x0) = x*dirac(x-x0);
    d = subs(d,x0,10);
    d = int(d,x,-inf,inf);
    
    
    d(x) = dirac(x);
    r(g) = int(g,x,-inf,inf);
    a = subs(r,g,d)
    int(d,x,-inf,inf);
    
    
    
    syms d x0 x X x1 x2 c F
    X(x1,x2) = (x1^.5)*(x2^.5);
    psi1(x1,c) = dirac(x1-c);
    psi2(x2,c) = dirac(x2-c);
    F(x1,x2,c) = (x1^.5)*(x2^.5)*dirac(x1-c)*dirac(x2-c);
    p = int(int(F,x1,-inf,inf),x2,-inf,inf);
    
    
    %%
    clear all
    syms f rho theta alpha x y xpos ypos pos npos p v w
    
    
    rot(alpha,x,y) = [[cos(alpha) -sin(alpha)];[sin(alpha) cos(alpha)]]*[x;y];
    
    
    
    xpos(rho,theta) = rho*cos(theta);
    ypos(rho,theta) = rho*sin(theta);
    
    
    
    
    pos(rho,theta) = [rho*cos(theta);rho*sin(theta)];
    new = subs(rot,[x y],[xpos ypos]);
    newN(alpha,rho,theta) = subs(rot,[x y],[xpos ypos]);
    what = compose(rot,xpos,x,[rho theta]);
    what = compose(what,ypos,y,[rho theta]);
    rot(alpha,npos) = [[cos(alpha) -sin(alpha)];[sin(alpha) cos(alpha)]]*npos;
    
    %%
    clear all
    syms y
    M = sym('M', [2 2]);
    v = sym('v', [2 1]);
    w = sym('w', [1 2]);
    y
    y(w)= w*M*v;
    
    clear all 
    syms w1 w2 v1 v2 a1 a2 b1 b2
    A = sym('A', [2 2]);
    B = sym('B', [2 2]);
    v = sym('v', [2 1]);
    w = sym('w', [1 2]);
    a = sym('a', [2 1]);
    b = sym('b', [1 2]);
    
    
    
    f(b)= b*A*a;
    g(w)= w*B*v;
    
    
     T = sym('T', [2 2 3]);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    