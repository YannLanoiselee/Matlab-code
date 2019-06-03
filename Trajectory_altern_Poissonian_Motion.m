function [ X,Y,list_model ] = Trajectory_altern_Poissonian_Motion( N,M,lambda )
% Trajectory : (N,M) vector containing M columns of trajectories
% List_Sigma : List of sigma value at a given time
% half_period : duration of a half period
% half_period_number : number of half periods

% Altern between
% Bm n_model=1
% OUP n_model=2
% fBm n_model=3
number_models=3;



%% fBm 1
sig_fBm1=1;H1=0.7;
%% fBm 2
sig_fBm2=1;H2=0.35;
%% Ornstein Uhlenbeck
sig_OUP=1;
k=1;
%% BM + drift
sig_BM1=1; drift1=0.2;

%% BM only
sig_BM2=0.5; drift2=0;
%% Sin
L=lambda;
Amp=2;
%%
X=zeros(N,M);
Y=zeros(N,M);
list_model=zeros(N,M);

for m=1:M
    list_time= cumsum(poissrnd(lambda,N,1));
    pos=list_time<=N;
    list_time=list_time(pos);
    if list_time(end,1)<N
        list_time=[list_time; N];
    end
    list_time=diff([0;list_time]);
    Position=0;% repére
    n_model=mod((1:size(list_time,1))',number_models)+1;%5*ones(size(list_time,1));%%randi([1 3],size(list_time,1),1);%3*ones(size(list_time,1));%
    
    for Ti=1:size(list_time,1)
        v=n_model(Ti,1);
        Timephase=list_time(Ti,1);
        if Position>0
            X_last=X(Position,m);
            Y_last=Y(Position,m);
        else
            X_last=0;
            Y_last=0;
        end
        
        
        if v==1 % fBm
            X(Position+1:Position+Timephase,m)=X_last+Trajectory_fBm( Timephase,H1,sig_fBm1^2/2);
            Y(Position+1:Position+Timephase,m)=Y_last+Trajectory_fBm( Timephase,H1,sig_fBm1^2/2 );
        elseif v==2 % fBm
            X(Position+1:Position+Timephase,m)=X_last+Trajectory_fBm( Timephase,H2,sig_fBm2^2/2);
            Y(Position+1:Position+Timephase,m)=Y_last+Trajectory_fBm( Timephase,H2,sig_fBm2^2/2 );
        elseif v==3 % OUP
            X(Position+1:Position+Timephase,m)= Trajectory_OUP(Timephase,1,k,sig_OUP,X_last,X_last );
            Y(Position+1:Position+Timephase,m)= Trajectory_OUP( Timephase,1,k,sig_OUP,Y_last,Y_last );
        elseif v==4 % Brownian
            angle=2*pi*rand();%
            
            traj_driftx= cos(angle)*cumsum(drift1*ones(Timephase,1));
            traj_drifty=sin(angle)*cumsum(drift1*ones(Timephase,1));
            
            X(Position+1:Position+Timephase,m)= traj_driftx+Trajectory_BM( Timephase,1,sig_BM1,X_last );
            Y(Position+1:Position+Timephase,m)= traj_drifty+Trajectory_BM( Timephase,1,sig_BM1,Y_last );
        elseif v==5 % fBm
            angle=2*pi*rand();%
            
            traj_driftx= cos(angle)*cumsum(drift2*ones(Timephase,1));
            traj_drifty=sin(angle)*cumsum(drift2*ones(Timephase,1));
            
            X(Position+1:Position+Timephase,m)= traj_driftx+Trajectory_BM( Timephase,1,sig_BM2,X_last );
            Y(Position+1:Position+Timephase,m)= traj_drifty+Trajectory_BM( Timephase,1,sig_BM2,Y_last );
            %             sens=(rand()-1/2);
            %             sens=sens./abs(sens);
            %                 X(Position+1:Position+Timephase,m)=X_last-Amp+Amp*cos(sens*(0:Timephase-1)*2*pi/L)';
            %                 Y(Position+1:Position+Timephase,m)=Y_last+Amp*sin(sens*(0:Timephase-1)*2*pi/L)';
            
        end
        
        list_model(Position+1:Position+Timephase,m)=v*ones( Timephase,1);
        Position=Position+Timephase;
    end
end
%%
    function [ Trajectory ] = Trajectory_OUP( N,M,k,sig,mu,x0 )
        %Generatesd an OUP Process
        dt = 1e-2;
        t = 1:dt:N; t=repmat(t',1,M);% Time vector
        % Set random seed
        x=x0;
        for i = 1:length(t)-1
            x(i+1,:) = x(i,:)+k*(mu-x(i,:))*dt+sig*sqrt(dt)*randn(1,M);
        end
        Trajectory=x(1:1/dt:end,:);
    end
%%
    function[Trajectory]=Trajectory_BM( N,M,sig,x0 )
        Trajectory=x0+cumsum(sig*randn(N,M),1);
    end
%%
    function [ Trajectory ] = Trajectory_fBm( N,H,D )
        % Gives a trajectory of fractionnal brownian motion
        % N : number of points
        % H : Hölder exponent
        % D : Diffusion coefficient
        % Trajectory : Column vector containing a brownian trajectory
        %              with given h and D
        
        % Davies and Harte method (fft)
        if H==0.5
            Trajectory=cumsum((2*D)^0.5*randn(N,1));
        else
            C=VecteurCovariance(N,H,D);
            
            L2=fft(C(:));
            
            %Calculus of Z=QL^0.5W
            Vk=randn(N,2);
            wk=zeros(2*N,1);
            wk(1,1)=(L2(1)/(2*N))^0.5*Vk(1,1);
            wk(N+1,1)=(L2(N+1)/(2*N))^0.5*Vk(1,2);
            for i=2:N
                wk(i,1)=(L2(i)/(4*N))^0.5*(Vk(i,1)+1i*Vk(i,2));
                wk(N+i,1)=(L2(N+i)/(4*N))^0.5*(Vk(2*N+2-(i+N),1)-1i*Vk(2*N+2-(i+N),2));
            end
            Z=real(fft(wk));% Vecteur de bruit gaussien fractionnaire
            Trajectory = cumsum(Z(1:N));%Mouvement brownien fractionnaire
        end
        function [ Vector ] = VecteurCovariance( N,H,D )
            %Function which gives the first column of the Covariance matrix (2Nx2N) d'un
            %Mouvement Brownien fractionnaire
            Vector=zeros(2*N,1); % Premiére ligne de la matrice d'autocovariance
            for ju=1:N
                Vector(ju,1)=AutoCovariance(ju-1,H,D);%i-1;% %de 0 a N-1
            end
            Vector(N+1,1)=0;
            for gu=N+2:2*N
                Vector(gu,1)=AutoCovariance(2*N-gu+1,H,D);%2*N-i+1;% %de N-1 à 0
            end
        end
        function [ Result ] = AutoCovariance(k,H,D)
            %Function that gives autocovariance of a fractionnal gaussian
            %noise
            % with k l'indice du saut et H l'exposant de Holder
            Result=D*((abs(k-1))^(2*H)-2*(abs(k))^(2*H)+(abs(k+1))^(2*H));
        end
        
    end
end