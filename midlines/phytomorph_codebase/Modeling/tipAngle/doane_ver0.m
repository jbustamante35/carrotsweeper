% load data
f{1} = '/home/nate/Downloads/151029 WT am copy 2.csv';
f{2} = '/home/nate/Downloads/151114 3.3-9am.csv';
f{3} = '/home/nate/Downloads/151001 082825 am.csv';
f{4} = '/home/nate/Downloads/151004 085606 am copy 2.csv';
for e = 1:numel(f)
    D(e).data = csvread(f{e});
    D(e).data(:,end) = [];
    if e == 1
        D(e).data(end,:) = [];
        D(e).data(:,3) = [];
    end
end
D(1).name = 'WT';
D(2).name = '3.3';
D(3).name = 'pot1';
D(4).name ='pot2';
%% get means and std
G = [];
M = [];
for e = 1:numel(D)
    U(:,e) = mean(D(e).data,2);
    S(:,e) = std(D(e).data,1,2)*size(D,2)^-.5;
    LEG{e} = D(e).name;
    G = [G e*ones(1,size(D(e).data,2))];
    M = [M D(e).data];
end

close all
figure;
plot(U);
legend(LEG);
figure;
errorbar(U,S);
legend(LEG);
figure;
plot(M);
figure;
plot(U)
%% pull on group X, Y
[mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL_T(M,size(M,1));
close all
figure
plot(mL);
[mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL_T(M,6);
plot(mS);
cp1 = find(G==1);
cp2 = find(G==2);
cp = [cp1 cp2];
subD = mC(:,cp);
subG = G(cp);
[lambda] = myLDA(subD',subG);
l1 = lambda'*mC(:,cp1);
l2 = lambda'*mC(:,cp2);
[h p] = ttest2(l1,l2)
full_lambda = PCA_BKPROJ_T(lambda,mE,mU) - mU;
figure;
plot(full_lambda);
figure;
plot(mE);
lambda/norm(lambda);
scale = 7;
X = scale*linspace(-max(l1),max(l1),100);
mu = mean(l1);
sigma = std(l1);
Y(1,:) = normpdf(X,mu,sigma);
mu = mean(l2);
sigma = std(l2);
Y(2,:) = normpdf(X,mu,sigma);
%% split others 1 - 3
cp1 = find(G==1);
cp2 = find(G==3);
l1 = lambda'*mC(:,cp1);
l2 = lambda'*mC(:,cp2);
[h p] = ttest2(l1,l2)
mu = mean(l2);
sigma = std(l2);
Y(3,:) = normpdf(X,mu,sigma);
%% split other 1 - 4
cp1 = find(G==1);
cp2 = find(G==4);
l1 = lambda'*mC(:,cp1);
l2 = lambda'*mC(:,cp2);
[h p] = ttest2(l1,l2)
mu = mean(l2);
sigma = std(l2);
Y(4,:) = normpdf(X,mu,sigma);
%% make fun plot
close all
plot(Y')
legend(LEG)
%% split wt vs all mut - carefull - reassign grouping vector
cp1 = find(G==1);
cp2 = find(G~=1);
G(cp1) = 1;
G(cp2) = 2;
%% is real - hold out
lambda = [];
R = 20;
N = 1;
h = [];
p = [];
for N = 1:6
    for rep = 1:100
        L1 = [];
        L2 = [];
        h1 = [];
        h2 = [];
        for r = 1:R
            cp1 = find(G==1);
            cp2 = find(G==2);
            %r1 = randi(numel(cp1),1,N);
            %r2 = randi(numel(cp2),1,N);
            r1 = randperm(numel(cp1));
            r2 = randperm(numel(cp2));
            r1 = r1(1:N);
            r2 = r2(1:N);
            h1(r,:) = cp1(r1);
            h2(r,:) = cp2(r2);
            cp1(r1) = [];
            cp2(r2) = [];
            cp = [cp1 cp2];
            subD = mC(:,cp);
            subG = G(cp);
            [lambda(:,r)] = myLDA(subD',subG);
        end

        for e = 1:size(lambda,2)
            if lambda(:,1)'*lambda(:,e) < 0
                lambda(:,e) = -lambda(:,e);
                flip = 1;
            end
        end
        %{
        close all
        figure;
        plot(mC(1,cp1),mC(2,cp1),'k.');
        hold on
        plot(mC(1,cp2),mC(2,cp2),'r.');
        for e = 1:R
            vec = [lambda(1,e),lambda(2,e)];
            vec = vec / norm(vec);
            quiver(vec(1),vec(2),200)
        end
        waitforbuttonpress
        %}
        for e = 1:size(lambda,2)
            l1 = lambda(:,e)'*mC(:,h1(e,:));
            l2 = lambda(:,e)'*mC(:,h2(e,:));
            %[mini_h(N,rep,e) mini_p(N,rep,e)] = ttest2(l1,l2);
            L1 = [L1 l1];
            L2 = [L2 l2];
        end
        [h(N,rep) p(N,rep)] = ttest2(L1,L2);
    end
end
%%




