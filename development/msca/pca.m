
function [T,P,ssq,a1,a2,a3,a4] = pca(X,center,~,ncomp)
if center
    X = X-mean(X);
end

[u,s,v] = svd(X, "econ");
T = u*s;
P = v;
ssq = diag(s).^2;
ssq = [ssq ssq./sum(ssq).*100 cumsum(ssq./sum(ssq)).*100];

T = T(:,1:ncomp);
P = P(:,1:ncomp);
ssq = ssq(1:ncomp,:);

a1 = 0; a2 = 0; a3 = 0; a4 = 0;