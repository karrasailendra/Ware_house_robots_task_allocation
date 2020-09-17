function out = ACOTW(car,pos,twcomb)
% prompt = 'Enter the number of passengers = ';
% N = input(prompt); 
% pos=100*rand(N*2,2);
% clear all
% clc
% pos=[5.62618912427925,72.2143728021138;44.3679286866576,82.9293071878729;29.6470403866946,24.0455898905690;0.0425542411259583,49.3188166862448;98.8521370163346,2.86878517907319;87.4655301269997,82.5943040673200;14.2617146436496,80.3818800669476;85.0566509436744,12.6219189828099;44.3783933010078,93.1199385332450;90.4875444571824,99.2084006947154;92.7397168816461,88.2598909045976;46.9534457195288,21.9782432431265;71.3354928154916,87.7083107253795;9.80398116732768,99.5052736198553;51.4412744586874,17.2349553305288;69.4834012315226,82.5662827388037;84.3986878655096,30.1815834398412;81.5686052867703,92.6136003928314;14.3009964850661,27.5207701089181;75.8341289433974,77.8374071244509];
% pos=[94.4365048138129,60.9792045709351;53.3989167664674,40.4214791574910;94.5952652761107,23.2318927975681;23.1650581282914,21.3948728531963;74.9694224876883,63.9074574528537;37.5944867961994,19.7053397090811;66.2295777417929,61.2069769012359;98.2708463182729,37.5963945444074;92.6782767130573,59.0935070157827;21.2231583804602,83.4996926162185];
% pos = [40,100,100,100,0,80,70,30,0,0; 0,0,30,50,100,0,100,100,80,30]';
% time=0;
% car=[89.8240902727620,41.9334745300729];
% car=[43.3783933010078 92.1199385332450];
% xp=pos(1:N,1);
% yp=pos(1:N,2);
% xd=pos(N+1:2*N,1);
% yd=pos(N+1:2*N,2);
N = length(pos(:,1))/2;
x = [car(1); pos(:,1)];
y = [car(2); pos(:,2)];
n = numel(x);


timespace= twcomb; % change this for function



% N=(n-1)/2;
    for i=1:n-1
        for j=i+1:n  
            D(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            D(j,i)=D(i,j);
        end
    end
%     for k=1:2*N-1
%     time = time+D(k+1,k);
%     end
%     time = time+ D(k+1,1);
    
%     for i = 1:2*N
%         timespace(i)=rand(1)*0.7*time;
%     end
% %     timespace=sort(timespace);
%     timespace=[52.2905797700914,10.1906367183525,117.638991189021,75.8571550062872,207.868653750089,371.994487617707,40.4044333837806,190.973712052215,0,353.076197825428,364.253950160035,135.068437184377,330.635849561065,454.515751978195,257.408301844122,392.478967799234,173.401724934767,341.983909200116,101.904381534482,384.561013073019];
%     timepool = 1:length(timespace);
%    timespace= [62.8639308188104;191.659578080091;144.480832889154;66.3807021973016;26.9756154588469;190.741651969206;80.2801103910090;225.261029255319;113.453927804677;161.408193666687];   
  
% 
    twcomb=[0,20];
    for i = 1:2*N
        tol=20;
        twcomb=[twcomb; timespace(i)-tol, timespace(i)+tol];
%         twcomb = [twcomb; 0 300];
    end
    val.n=n;
    val.x=x;
    val.y=y;
    val.D=D;

Cost=@(tour) TourLength(tour,val);

nVar=val.n;
% conditions
it = 80;      
Ant  =5*val.n;      % check for number of `ants
Q     = 100;
tau0=10^0;	
alpha=1;        
beta=0;               
e = 5;
gamma=1;
rho=0.3; % evapouration rate


% Initialization
eta=1./val.D;          
tau=tau0*ones((nVar),(nVar)); 
Bestpath=zeros(it,1);    
empty_ant.Tour=[];
empty_ant.Cost=[];
ant=repmat(empty_ant,Ant,1);
BestSol.Cost=inf;
out.Cost = 0;

% ACO Main Loop

for it=1:it
    
    % Move Ants
    shit = 0;
    for k=1:Ant  
%         disp(k)
        % pool for requests
        nonpool = N+2:(2*N)+1;
        ant(k).Tour=1; 
        pool = 1:N+1;
%         nonpool(nonpool==ant(k).Tour(end)+N)=[];
        omega=[];
        troll = 1;
        
        for l=2:2*(N)+1           
            i=ant(k).Tour(end);
            for a1 = 1: 2*N+1
                if twcomb(a1,2)>= timecost(ant(k).Tour,twcomb,D)+ D(i,a1)
                  omega(i,a1)= exp((((timecost(ant(k).Tour,twcomb,D)+ D(i,a1)-twcomb(a1,1)))/(twcomb(a1,2)-(timecost(ant(k).Tour,twcomb,D) + D(i,a1))+1)));
                else 
                    if a1 ~= 1
                    end
                  omega(i,a1) = 0;  
                end
            end
%             disp(omega)
%             P1=tau(i,:).^alpha.*eta(i,:).^beta;
%             max(omega(i,:))
            P=(tau(i,:).^alpha.*eta(i,:).^beta.*omega(i,:).^gamma);
%             P=tau(i,:).^alpha.*omega(i,:).^gamma;
%             P=[1 P];
%             disp(P)
            P(ant(k).Tour)=0; 
            P(nonpool)=0;
%             if max(P)==inf
%            disp(P)     
%             end
            if sum(P) == 0
                troll = 0;
                shit=shit+1;
%                 disp(l)
                break
            end
            ing = P;
            P=(P/sum(P));           
            j=RouletteWheel(P);
            if isempty(j)
                j = find(ing==inf);
            end
            ant(k).Tour=[ant(k).Tour j];
%             disp(j)
            if j<=N+1
            pool = [pool j+N]; 
            end
            if isempty(j)
            end
            nonpool(nonpool==j+N)=[];
        end       
        if troll == 0
            troll = 1;
                
        else
        ant(k).Cost=timecost(ant(k).Tour,twcomb,D);    
%         disp(k)
        if ant(k).Cost<BestSol.Cost
            BestSol=ant(k);
        end 
        end
         
       
        tour=ant(k).Tour;
        for l=1:length(tour)-1
            i=tour(l);
            j=tour(l+1);
%             ind = find(bst==i,1);
%             if (j==bst(ind+1))
%                 tau(i,j)=(1-rho)*tau(i,j)+ Q/ant(k).Cost + (e*Q/BestSol.Cost);
%             else
                tau(i,j)=(1-rho)*tau(i,j);
%             end
        end
        if length(tour) == 11
        end
        
        
    end
    
    %disp(BestSol)
    tau=(1-rho)*tau;
    
    %  Pheromone Update
    if (BestSol.Cost) > 1000000
        out.Cost = 500000;
        out.Tour = [];
        break
    else
    for k=1:Ant
%         if length(ant(k).Tour)== (2*N)+1
         tour=ant(k).Tour;
        bst = BestSol.Tour;
%         tour=[tour tour(1)]; 
   if isempty(ant(k).Cost)
     ant(k).Cost = inf;
    end
%         bst=[bst bst(1)];
        for l=1:length(tour)-1
            i=tour(l);
            j=tour(l+1);
%             ind = find(bst==i,1);
%             if (j==bst(ind+1))
%                 tau(i,j)=(1-rho)*tau(i,j)+ Q/ant(k).Cost + (e*Q/BestSol.Cost);
%             else
                tau(i,j)=(1-rho)*tau(i,j)+Q/ant(k).Cost;
%             end
        end
%         end
    end
    end
    Bestpath(it)=BestSol.Cost;
    out.Tour=[];
%     disp(['Iteration ' num2str(it) ': Best path length = ' num2str(Bestpath(it))]);
%     figure(1);
%      tour=[BestSol.Tour BestSol.Tour(1)];
% %      disp(BestSol.Tour)
%     plot(val.x(tour),val.y(tour),'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6])
%     
%     xlabel('x');
%     ylabel('y');
%     grid on;
%     drawnow;
%     Ant  = 2*val.n;
end
% Results
% figure;
% plot(Bestpath)
% xlabel('Iteration number');
% ylabel('length of Best Tour');
% grid on;


 
% check=interval(BestSol.Tour,D);
% if isempty(BestSol.Tour)
%   out.cool = 0;
% else
%       checkmat = timecostmat(BestSol.Tour,twcomb,D);
%     cool=0;
% 
%     for i = 2:length(checkmat) + 1
%          checklist(BestSol.Tour(i)-1)=checkmat(i-1);
%     end
% 
%     for i = 1:length(checkmat)
%        if  checklist(i) > twcomb(i+1,1) && checklist(i) < twcomb(i+1,2)
%            cool = cool +1;
%        end
%     end
% 
%     out.cool = cool;
% end
% timecostbest = timecost(BestSol.Tour,twcomb,D);
% for i=1:length(tour)-1
% if check(i)<=twcomb(tour(i),2) && check(i)>=twcomb(tour(i),1)
% 
%     dire(i)= 1;
% else
%     dire(i)= 0;
% end
% end

if out.Cost == 500000
else
out.Cost = BestSol.Cost;

end

end
    
function L=TourLength(tour,model)

    n=numel(tour);

%     tour=[tour tour(1)];
    
    L=0;    
    for i=1:n-1
        L=L+model.D(tour(i),tour(i+1));
    end

end
function j=RouletteWheel(P)

    r=rand;
    
    C=cumsum(P);
    
    j=find(r<=C,1,'first');
    
%     k=max(P);
%     j=find(P==k);

end

function out = interval(tour,D)
out=zeros(1,length(tour(1:end-1)));
for i =1 :length(tour(1:end-1))
    for k = 1:i-1
        out(i)=D(tour(k),tour(k+1))+out(i);
    end
end
end

function out = timecost(tour,twcomb,D)
out=0;
  for i =1:length(tour)-1
   if out >= twcomb(tour(i),1) && out <= twcomb(tour(i),2)
       out = out+D(tour(i),tour(i+1));
   elseif out < twcomb(i,1)
       out = twcomb(i,1);
   else
       out =inf;
   end
  end
end

function outmat = timecostmat(tour,twcomb,D)
out=0;
outmat=[];
  for i =1:length(tour)-1
   if out >= twcomb(tour(i),1) && out <= twcomb(tour(i),2)
       out = out+D(tour(i),tour(i+1));
   elseif out < twcomb(i,1)
       out = twcomb(i,1);
   else
       out =inf;
   end
   outmat=[outmat out];
  end
end