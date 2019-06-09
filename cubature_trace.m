function [B_t]=cubature_trace(mu,var,y,obj,satPos, satRefPos,waveLengVec2,basePosition)
S=chol(var,'lower');
n=length(mu);
univec1=eye(n);
univec2=-eye(n);
univec=sqrt(n).*[univec1 univec2];
mumatrix=repmat(mu,1,2*n);
samp=mumatrix+S*univec;
B_t=0;
for i=1:2*n
    x=samp(:,i);
    h_x_new_new = observation_est( x, satPos, satRefPos, x(obj.sizeState_+1:end), waveLengVec2, obj.D_, basePosition); % Observation model
    B_t=B_t+(y-h_x_new_new)'*inv(obj.R_)*(y-h_x_new_new);
end
B_t=B_t./(2*n);
samp_mean=mean(samp,2);
samp_var=0;
for i=1:2*n
    samp_var=samp_var+(samp(:,i)-samp_mean)*(samp(:,i)-samp_mean)';
end
samp_var=samp_var./(2*n);

end