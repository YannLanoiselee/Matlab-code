function [Local_Measures,S,Z] = Local_measures_calc(X,Y,T)
%UNTITLED Summary of this function goes here
%% calculation of the local measures
N=size(X,1);
%% volume and anisotropy of ocnvex hull
method={'diameter','volume','anisotropy'}; 
Local_Measures=zeros(N,4);
cp_test=0;
for n_method=1:2
cp_test=cp_test+1;
    [~,~,Q_n] = Test_Convex_Hull([X(1:N,1),Y(1:N,1)] ,T/2,method{1,n_method});
Local_Measures(:,cp_test)=Q_n;
end
%% local D
List_D_Local=Local_D([X,Y],T);

Local_Measures(:,3)=List_D_Local.^0.5;
%% local drift

Local_Measures(:,4)=Local_drift([X,Y],T);

%% local exponent
[List_Exponent_Local,~]=TAMSD_Trajectory_Local([X,Y],T,T/2);
Local_Measures(:,5)=List_Exponent_Local;

for i=T+1:N-T
    Local_Measures(i-T:i+T,3:5)=Local_Measures(i-T:i+T,3:5)+repmat(Local_Measures(i,3:5),2*T+1,1)./(2*T+1);%Time(1,pos_T);%mean(Historique_Hull(i,:));
end
%% calculation of the distance matrix
Z=Local_Measures(T:end-T,:);
% Z=Z.*repmat(std(Z),size(Z,1),1)./repmat(mean(Z),size(Z,1),1)./repmat(std(diff(Z)),size(Z,1),1);
Z=Z./repmat(mean(Z),size(Z,1),1);
S=zeros(size(Z,1),size(Z,1));
for i=1:size(Z,1)
    S(:,i)=(sum((Z-repmat(Z(i,:),size(Z,1),1)).^2,2)).^0.5;
%     S(:,i)=abs(repmat(S_n(i,1),N,1)-S_n);
end


end

