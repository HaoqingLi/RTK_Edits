function [samp,dif,samp_mean,samp_var]=cubature_linear(mu,var,y,H,real_B_t)
S=chol(var);
n=length(mu);
univec1=eye(n);
univec2=-eye(n);
univec=sqrt(n).*[univec1 univec2];
mumatrix=repmat(mu,1,2*n);
samp=mumatrix+S*univec;

samp_mean=mean(samp,2);
samp_var=0;
for i=1:2*n
    samp_var=samp_var+(samp(:,i)-samp_mean)*(samp(:,i)-samp_mean)';
end
samp_var=samp_var./(2*n);

B_t=0;
for i=1:2*n
    B_t=B_t+(y-H*samp(:,i))*(y-H*samp(:,i))';
%      B_t=B_t+(samp(:,i))*(samp(:,i))';
end
B_t=B_t./(2*n);
dif=real_B_t-B_t;

end