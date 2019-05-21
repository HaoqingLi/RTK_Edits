function [obj]=KF_OD_loop(obj,z,h_z,H,basePosition,iGNSS)
%% Prediction
       basePosition = [4041839.1018   537121.6018  4888452.5105]; 
       obj.state_(1:3)=obj.state_(1:3)-basePosition';
       x_old_old=obj.state_;
       tau=1;
       z_i=1;
       e0=0.9;
       f0=0.1;
       ti=1;
       shreshold=1e-7;
  %% Update
  
     while(tau>shreshold)

         if(z_i<1e-10)
           x_new_new=obj.state_;
           P_new_new=obj.P_;
         else
            obj.S_ = H*obj.P_*H' + obj.R_; 
            obj.K_ = obj.P_*H'/(obj.S_+obj.R_./z_i);   
%             x_new_new=obj.state_+ obj.K_*(z-h_z);
            x_new_new=obj.state_+ obj.K_*(z-H*obj.state_);
            dif=z-H*obj.state_;
            P_new_new=obj.P_- obj.K_*(obj.S_+obj.R_./z_i)* obj.K_';
         end
            e=e0+z_i;
            f=f0+1-z_i;
            b_t_hat=z*z'+H*(P_new_new+x_new_new*x_new_new')*H'-z*(H*x_new_new)'-(H*x_new_new)*z';
%             [samp,dif,samp_mean,samp_var]=cubature(x_new_new,P_new_new,z,H,b_t_hat);
            pi_t=psi(e)-psi(e+f);
            pi_1_t=psi(f)-psi(e+f);
            P_1=exp(pi_t-0.5*trace(b_t_hat/(obj.R_)));
            P_0=exp(pi_1_t);
            z_i=(P_1)/(P_1+P_0);
            if(ti>1)
            tau=norm(x_new_new-x_old_old)/norm(x_old_old);
            end
            x_old_old=x_new_new;
            P_old_old=P_new_new;
            ti=ti+1;
            if(iGNSS<10)
                z_i=1;
            end
     end
        x_old_old(1:3)=x_old_old(1:3)+basePosition';
        obj.state_= x_old_old;
        obj.P_ = P_old_old;
          


end