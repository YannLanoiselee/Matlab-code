function [List_D_Local] = Local_D(Trajectory,T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Calculates the local exponent \alpha by calculating over T steps the power law
%exponent of the TAMSD over windows of size [1,WindowMax]
%The first T/2 and last T/2 values exponent are not calculated
% 
    N=size(Trajectory,1);
    
    D=zeros(N,1);
    for j=T/2+1:N-T/2
     D(j,1)=sum(var(diff(Trajectory(j-T/2:j+T/2,:),1,1),1)./2);%max([Resultatx(1,1),Resultaty(1,1),Resultatz(1,1)]);
    end
    List_D_Local=D;
    
end

