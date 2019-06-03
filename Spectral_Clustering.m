function [IDX] = Spectral_Clustering(W)


% A tutorial on spectral clustering
% ulrike von Luxburg
% normalized spectral clustering 
% fom Jordan and Weiss

%% Compute the degree matrix
    D = diag(sum(W,1));
%% Matrix L
L=eye(size(W,1))-D^(-1/2)*W*D^(-1/2);
%%eigen value decomposition
[eigVectors,eigValues] = eig(L);
%% Largest variation in the cluster numbers
e_val=sum(abs(eigValues),2);
[~,idx_cl]=max(abs(e_val(2:end,1)-e_val(1:end-1,1)));
cluster_number=idx_cl;%
neig=cluster_number;
%% select cluster_number largest eigen vectors
    nEigVec = eigVectors(:,1:neig);
%% construct the normalized matrix U from the obtained eigen vectors
for i=1:size(nEigVec,1)
    upo = sqrt(sum(nEigVec(i,:).^2));    
    Z(i,:) = nEigVec(i,:)./ upo; 
end
%% clustering
[IDX,~] = kmeans(Z,neig);%eva.OptimalK); 

end

