function   [log_lh,mahalaD]=wdensity_phyto(X, mu, Sigma, p, sharedCov, CovType,L,logDetSigma)
%WDENSITY Weighted conditional density and mahalanobis distance.
%   LOG_LH = WDENSITY(...) returns log of component conditional density
%   (weighted by the component probability) of X. LOG_LH is a N-by-K matrix
%   LOG_LH, where K is the number of Gaussian components. LOG_LH(I,J) is
%   log (Pr(point I|component J) * Prob( component J))
%
%   [LOG_LH, MAHALAD]=WDENSITY(...) returns the Mahalanobis distance in
%   the N-by-K matrix MAHALAD. MAHALAD(I,J) is the Mahalanobis distance of
%   point I from the mean of component J.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2011/05/09 01:28:14 $

    log_prior = log(p);
    [n,d]=size(X);
    k=size(mu,1);
    log_lh = zeros(n,k);
    mahalaD = zeros(n,k);
    %logDetSigma = -Inf;
    for j = 1:k
        

        Xcentered = bsxfun(@minus, X, mu(j,:));
       
        if CovType == 2
            xRinv = Xcentered /L{j} ;
        else
            xRinv = bsxfun(@times,Xcentered , (1./ L{j}));
        end
        mahalaD(:,j) = sum(xRinv.^2, 2);

        log_lh(:,j) = -0.5 * mahalaD(:,j) +...
            (-0.5 *logDetSigma(j) + log_prior(j)) - d*log(2*pi)/2;
        %get the loglikelihood for each point with each component
        %log_lh is a N by K matrix, log_lh(i,j) is log \alpha_j(x_i|\theta_j)
    end

   