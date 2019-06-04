function [List_Drift_local] = Local_drift(Trajectory,T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Calculates the local exponent \alpha by calculating over T steps the power law
%exponent of the TAMSD over windows of size [1,WindowMax]
%The first T/2 and last T/2 values exponent are not calculated
% 
    N=size(Trajectory,1);
    
    Drift=zeros(N,1);
    for j=T/2+1:N-T/2
     Drift(j,1)=(sum((Trajectory(j+T/2,:)-Trajectory(j-T/2,:)).^2)).^0.5;%max([Resultatx(1,1),Resultaty(1,1),Resultatz(1,1)]);
    end
    List_Drift_local=Drift;
end

