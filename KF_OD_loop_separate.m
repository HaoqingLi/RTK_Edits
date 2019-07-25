function [obj,z_i]=KF_OD_loop_new(obj,z,h_z,H,basePosition,iGNSS,z_DDDelta,satPos, satRefPos,waveLengVec2)
%% Prediction
       basePosition = [4041839.1018   537121.6018  4888452.5105]; 
%        obj.state_(1:3)=obj.state_(1:3)-basePosition';
       x_old_old=obj.state_;
       tau=1;
       z_i=ones(1,2);

       e0_C=0.995;
       f0_C=0.005;
       e0_L=0.92;
       f0_L=0.08;

       ti=1;
       shreshold=1e-6;

  %% Update
  
     while(tau>shreshold)
         

         if(any(z_i<1-1e-4))
           x_new_new=obj.state_;
           P_new_new=obj.P_;
         else
            R=obj.R_;
            R(1:length(R)/2,1:length(R)/2)=R(1:length(R)/2,1:length(R)/2)./z_i(1);
            R(length(R)/2+1:end,length(R)/2+1:end)=R(length(R)/2+1:end,length(R)/2+1:end)./z_i(2);
            obj.S_ = H*obj.P_*H' ;
            obj.K_ = obj.P_*H'/(obj.S_+R);   
            x_new_new=obj.state_+ obj.K_*(z-h_z);
            dif=z-H*obj.state_;
            P_new_new=obj.P_- obj.K_*(obj.S_+R)* obj.K_';
         end
         

            e_C=e0_C+z_i(1);
            f_C=f0_C+1-z_i(1);

            e_L=e0_L+z_i(2);
            f_L=f0_L+1-z_i(2);

         
            z=z+z_DDDelta;
            x_new_new(1:3)=x_new_new(1:3)-basePosition';
%             b_t_hat=z*z'+H*(P_new_new+x_new_new*x_new_new')*H'-z*(H*x_new_new)'-(H*x_new_new)*z';
            x_new_new(1:3)=x_new_new(1:3)+basePosition';
            [samp,samp_mean,samp_var,b_t_hat]=cubature(x_new_new,P_new_new,z,obj,satPos, satRefPos,waveLengVec2,basePosition);
            
%            tr=cubature_trace(x_new_new,P_new_new,z,obj,satPos, satRefPos,waveLengVec2,basePosition);

            pi_t_C=psi(e_C)-psi(e_C+f_C);
            pi_1_t_C=psi(f_C)-psi(e_C+f_C);
            pi_t_L=psi(e_L)-psi(e_L+f_L);
            pi_1_t_L=psi(f_L)-psi(e_L+f_L);
            R_KF=obj.R_;
            R_KF_C=1e0.*R_KF(1:length(R_KF)/2,1:length(R_KF)/2);
            R_KF_L=1e0.*R_KF(length(R_KF)/2+1:end,length(R_KF)/2+1:end);
            b_t_C=b_t_hat(1:length(R_KF)/2,1:length(R_KF)/2);
            b_t_L=b_t_hat(length(R_KF)/2+1:end,length(R_KF)/2+1:end);
            

            P_1_C=exp(pi_t_C-0.5*trace(b_t_C/(R_KF_C)));
            P_0_C=exp(pi_1_t_C);
            z_i(1)=(P_1_C)/(P_1_C+P_0_C);
            
            
            P_1_L=exp(pi_t_L-0.5*trace(b_t_L/(R_KF_L)));
            P_0_L=exp(pi_1_t_L);
            z_i(2)=(P_1_L)/(P_1_L+P_0_L);
            
            if(ti>1)
            tau=norm(x_new_new-x_old_old)/norm(x_old_old);
            end
            z=z-z_DDDelta;
            x_old_old=x_new_new;
            P_old_old=P_new_new;
            ti=ti+1;
             if(iGNSS==6)
                1;
             end
     end
        R=obj.R_;
        z_i_=z_i;
        obj.state_= x_old_old;
        obj.P_ = P_old_old;
          


end