    
clc;
clear all;
close all;
% prompt = 'Enter the number of passengers = ';
% N = input(prompt); 
carmat = [2 3];
cars=[];
pcombpos=[];
dcombpos=[];
carnum=[];
combtimespace=[];
for loops = 1:3
N = carmat(randi(1));
pos=[];
eta=[];
tau=[];
val=[];
timespace=[];
tour=[];
car=100*rand(1,2);
pos=100*rand(N*2,2);
D=[];


% pos=[5.62618912427925,72.2143728021138;44.3679286866576,82.9293071878729;29.6470403866946,24.0455898905690;0.0425542411259583,49.3188166862448;98.8521370163346,2.86878517907319;87.4655301269997,82.5943040673200;14.2617146436496,80.3818800669476;85.0566509436744,12.6219189828099;44.3783933010078,93.1199385332450;90.4875444571824,99.2084006947154;92.7397168816461,88.2598909045976;46.9534457195288,21.9782432431265;71.3354928154916,87.7083107253795;9.80398116732768,99.5052736198553;51.4412744586874,17.2349553305288;69.4834012315226,82.5662827388037;84.3986878655096,30.1815834398412;81.5686052867703,92.6136003928314;14.3009964850661,27.5207701089181;75.8341289433974,77.8374071244509];
% pos = [72.7349580929689,90.6673424927724;45.9411338748033,88.2143323876479;23.5417502459506,94.7244907813225;70.3620501148626,41.0331547317111;90.3354665098209,56.2650043097417;37.8988673199119,46.3829004866527;34.7316578852631,17.0068281107392;40.6624667905148,84.3897472719139;93.3109383399309,48.1242487160148;18.2816974646092,16.4297491328674;47.0264675147589,77.4953244442090;0.0149817115251194,26.8802609576352;12.9875852605737,38.6217779666874;80.2092466507340,34.7706298031206;66.5129707486672,39.1254594175180;22.5124129029250,73.6064329605265;22.7360302058110,78.8034489268214;95.1170741827675,18.5474076468766;48.8997389329240,0.540358304318545;66.9987348780218,14.4508856229825];
% xp=pos(1:N,1);
% yp=pos(1:N,2);
% xd=pos(N+1:2*N,1);
% yd=pos(N+1:2*N,2);

x = [car(1);  pos(:,1)];
y = [car(2);  pos(:,2)];
n = numel(x);
    for i=1:n-1
        for j=i+1:n  
            D(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            D(j,i)=D(i,j);
        end
    end
    val.n=n;
    val.x=x;
    val.y=y;
    val.D=D;

Cost=@(tour) TourLength(tour,val);

nVar=val.n;
% conditions
it = 20;      
Ant  = val.n;      % check for number of `ants
Q     = 100;
tau0=10^-6;	
alpha=1;        
beta=5;         
rho=0.7;       
e = 5;


% Initialization
eta=1./val.D;          
tau=tau0*ones(nVar,nVar); 
Bestpath=zeros(it,1);    
empty_ant.Tour=[];
empty_ant.Cost=[];
ant=repmat(empty_ant,Ant,1);
BestSol.Cost=inf;


% ACO Main Loop

for it=1:it
    
    % Move Ants
    for k=1:Ant        
        % pool for requests 
        nonpool = N+2:(2*N)+1;
        ant(k).Tour = 1;
%         ant(k).Tour=randi([1 N]); 
        pool = 1:N+1;
%         nonpool(nonpool==ant(k).Tour(end)+N)=[];
        for l=2:(2*N)+1           
            i=ant(k).Tour(end);            
            P=tau(i,:).^alpha.*eta(i,:).^beta;            
            P(ant(k).Tour)=0; 
            P(nonpool)=0;
            P=P/sum(P);           
            j=RouletteWheel(P);          
            ant(k).Tour=[ant(k).Tour j];
            if j<=N
            pool = [pool j+N]; 
            end
            nonpool(nonpool==j+N)=[];
        end       
        ant(k).Cost=Cost(ant(k).Tour);       
        if ant(k).Cost<BestSol.Cost
            BestSol=ant(k);
        end     
    end
    %disp(BestSol)
    tau=(1-rho)*tau;
    
    %  Pheromone Update
    for k=1:Ant
        tour=ant(k).Tour;
        bst = BestSol.Tour;
        tour=[tour tour(1)]; 
        bst=[bst bst(1)];
        for l=1:nVar
            i=tour(l);
            j=tour(l+1);
            ind = find(bst==i,1);
            if (j==bst(ind+1))
                tau(i,j)=tau(i,j)+ Q/ant(k).Cost + (e*Q/BestSol.Cost);
            else
                tau(i,j)=tau(i,j)+Q/ant(k).Cost;
            end
        end
    end
    

    
    Bestpath(it)=BestSol.Cost;
%     disp(['Iteration ' num2str(it) ': Best path length = ' num2str(Bestpath(it))]);
%     figure(1);
     tour=[BestSol.Tour BestSol.Tour(1)];
%      disp(tour)
%     plot(val.x(tour),val.y(tour),'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6])
%     
%     xlabel('x');
%     ylabel('y');
%     grid on;
%     drawnow;
    
end
% Results
% figure;
% plot(Bestpath)
% xlabel('Iteration number');
% ylabel('length of Best Tour');
% grid on;

check=interval(tour,D);


for i =1:(length(tour)-1)
   timespace(tour(i)) = check((i));
end

cars=[car; cars];
pcombpos=[pos(1:N,:); pcombpos];
dcombpos=[pos(N+1:end,:); dcombpos];
combtimespace=[timespace(2:end)'; combtimespace];
carnum = [N; carnum];
end
combpos=[pcombpos; dcombpos];
function L=TourLength(tour,model)

    n=numel(tour);

    tour=[tour tour(1)];
    
    L=0;
    for i=1:n
        L=L+model.D(tour(i),tour(i+1));
    end

end
function j=RouletteWheel(P)

    r=rand;
    
    C=cumsum(P);
    
    j=find(r<=C,1,'first');

end

function out = interval(tour,D)
out=zeros(1,length(tour(1:end-1)));
for i =1 :length(tour(1:end-1))
    for k = 1:i-1
        out(i)=D(tour(k),tour(k+1))+out(i);
    end
end
end