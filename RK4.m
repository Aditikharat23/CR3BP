function [t_track,state_track] = RK4(state,h,tend,u) 
%State input and output is a row vector
t_track(1)=0;
state_track(1,:)= state;
t = 0;
iterations = floor(tend/h) +1;
for i  = 1:iterations
    
k1 = h*fn.cr3bp(t,state,u);
k2 = h*fn.cr3bp(t,state + 0.5*k1,u);
k3 = h*fn.cr3bp(t,state + 0.5*k2,u);
k4 = h*fn.cr3bp(t,state+ k3,u);


state = state + (k1+2*k2+2*k3+k4)/6 ;

state_track(i+1,:) = state;
t = t+h;
t_track(end+1) = t;
end
end