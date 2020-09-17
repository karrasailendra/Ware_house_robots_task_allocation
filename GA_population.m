clc;
clearvars -except combpos combtimespace cars
close all;  
% prompt = 'Enter the number of passengers = ';
% P = 30; % input(prompt); 
P = length(combpos)/2;
% prompt = 'Enter the number of cars = ';
% C = 7; % input(prompt);  
C = length(cars);
% pos=100*rand(P*2,2);
% pos = [54.9280929584129,0.0957510931683614;3.98483910996506,18.8467220299546;55.1829849040330,75.2805161379527;68.1007116477246,62.1258881083963;28.5047622473183,84.4815252412195;72.1894936924880,20.9143168102400;2.11975006329034,63.8867548113356;66.7926777384821,2.38784058081196;52.3254759287707,1.87883440705329;31.5789811719396,82.8968685630402;28.0043767027660,69.4163292304162;39.4877227805107,74.3424520742846;41.4888897697545,8.09922825028926;85.8798202571836,35.4655936724098;93.1446387911041,26.5020793445935;33.9055258100080,34.9547665867157;35.4919517584656,40.2224857414105;47.3672096618537,97.1126961418939;98.4212496528695,39.6430537331353;46.2936151075491,5.27669598215339;38.2356743391658,56.2589751301206;28.0695400557028,5.28728399013004;17.7401623029976,62.1615069337770;18.1269466004408,89.7197856613693;32.8861021452982,65.7074097875966;51.1623840165549,13.5706746706972;0.651649719757252,49.9123403399490;45.4019106804204,50.0721251822814;40.5984607815625,71.6922139280067;44.5738259237303,1.79183266917720;90.2727582006921,69.0735210616644;63.1369307952129,72.6784782385606;40.3043985284266,69.0779753437661;10.7234165649967,18.6857409298853;34.0982504614510,92.0225178152608;72.4490596460086,33.9597230680280;50.1701094170752,61.5780579620943;86.2292272722789,87.2495071821673;81.4105654128070,91.0579295537315;4.51064480659030,79.3088417577659;16.1103108799612,68.2280918609845;90.1057326220262,96.9945495716043;97.5567006121722,65.2075898119630;37.0554193383152,98.3530263348776;26.3078547920594,30.0908386674635;37.8959942095892,3.26897370593465;39.2265291855669,96.9307638996925;81.8689624971681,48.0715044742785;48.5664657457337,88.2694222463659;33.4862112316332,32.4170546363623;11.0266555338439,31.6334338794788;58.5473090012815,95.9469907836223;24.4449102907890,5.39994493028778;8.51456660873603,92.3507118946229;48.9239325438764,32.9786978237172;5.48846259415408,16.1982552034487;4.08470976182047,70.7562570540990;87.7484928570313,0.359956066255684;11.6445630896545,69.3641410711331;6.69570205766585,88.6558954653879];

pos = combpos;
population=100;
new_generation=0.3*population;
fitness_selection = population;
itern = 400;
mut=0;

for i=1:population
parent(i).mat=zeros(C,P);
end
l=0;

 for i=1:P
     for l = 1:population
     k=randi(C);
     parent(l).mat(k,i)=1;
     end
 end
for i= 1:(population) % check for second iteration for child 3 and 4

    max_sum_parent = max(sum(parent(i).mat));
    min_sum_parent = min(sum(parent(i).mat));
    while max_sum_parent ~= 1 || min_sum_parent == 0
         parent(i).mat =  horicorrection(parent(i).mat);
%          l=l+1;
%          disp(max(sum(child(i).mat')))
         max_sum_parent = max(sum(parent(i).mat));
         min_sum_parent = min(sum(parent(i).mat));
    end
    
    max_sum_parent = max(sum(parent(i).mat'));
    while max_sum_parent> 3
         parent(i).mat =  correction(parent(i).mat);
%          l=l+1;
%          disp(max(sum(child(i).mat')))
         max_sum_parent = max(sum(parent(i).mat'));
    end
    
end 
for i = 1 : population
    parent(i).tour = tournum(parent(i).mat,C,P);
    mat11 = [];
    temptw=[];
    temp=[];
    for k = 1: C
    mat11 = cell2mat(struct2cell(parent(i).tour(k)));
    for g1 = 1:length(mat11)
        mat11 = [mat11 mat11(g1)+length(combpos)/2];
    end
    temp = routetoxy(mat11,pos,P);
    for g = 1:(length(mat11))
%         if g <= (length(mat11)/2)
            temptw(g,1) = combtimespace(mat11(g));
%         else
%             temptw(g,1) = combtimespace(mat11(g)+length(mat11)) ; 
%         end
    end
%     child(i).xy(:,1,k)=temp(:,1);
%     child(i).xy(:,2,k)=temp(:,2);
if isempty(temp)
   opttourlength(k) = 0;
else
   out = ACOTW(cars(k,:),temp,temptw);
   opttourlength(k)=out.Cost;
   parent(i).tour(k).xy=out.Tour;
end
    end
    parent(i).tourlength=sum(opttourlength);
end
T = struct2table(parent);
sortedT = sortrows(T, 'tourlength'); 
parent = table2struct(sortedT);
for iter = 1 : itern 
    tourlengthmat=[];
    opttourlength=[];
    child=[];
    inctour=[];
    nextchild=[];
    temp=[];
    tempchild=[];
    out=[];
%  s1=sum(parent(1).mat');
%  s2=sum(parent(2).mat');
% for i= 1:(population) % check for second iteration for child 3 and 4
%     max_sum_parent = max(sum(parent(i).mat'));
%     while max_sum_parent> 5
%          parent(i).mat =  correction(parent(i).mat);
% %          l=l+1;
% %          disp(max(sum(child(i).mat')))
%          max_sum_parent = max(sum(parent(i).mat'));
%     end
%     max_sum_parent = max(sum(parent(i).mat));
%     min_sum_parent = min(sum(parent(i).mat));
%     while max_sum_parent ~= 1 || min_sum_parent == 0
%          parent(i).mat =  horicorrection(parent(i).mat);
% %          l=l+1;
% %          disp(max(sum(child(i).mat')))
%          max_sum_parent = max(sum(parent(i).mat));
%          min_sum_parent = min(sum(parent(i).mat));
%     end
%     
% end 

% for i = 1 : population
%     parent(i).tour = tournum(parent(i).mat,C,P);
%     mat11 = [];
%     temptw=[];
%     temp=[];
%     for k = 1: C
%     mat11 = cell2mat(struct2cell(parent(i).tour(k)));
%     for g1 = 1:length(mat11)
%         mat11 = [mat11 mat11(g1)+length(combpos)/2];
%     end
%     temp = routetoxy(mat11,pos,P);
%     for g = 1:(length(mat11))
% %         if g <= (length(mat11)/2)
%             temptw(g,1) = combtimespace(mat11(g));
% %         else
% %             temptw(g,1) = combtimespace(mat11(g)+length(mat11)) ; 
% %         end
%     end
% %     child(i).xy(:,1,k)=temp(:,1);
% %     child(i).xy(:,2,k)=temp(:,2);
% if isempty(temp)
%    opttourlength(k) = 0;
% else
%    out = ACOTW(cars(k,:),temp,temptw);
%    opttourlength(k)=out.Cost;
%    parent(i).tour(k).xy=out.Tour;
% end
%     end
%     parent(i).tourlength=sum(opttourlength);
% end
T = struct2table(parent);
sortedT = sortrows(T, 'tourlength'); 
parent = table2struct(sortedT);

for j = 1:2: new_generation
    a1= randi(fitness_selection);
    a2= randi(fitness_selection);
    if rand(1)>0
    out = horisplit(parent(a1).mat,parent(a2).mat);
    else
        out = vertisplit1(parent(a1).mat,parent(a2).mat);
    end
    if rand(1)<0.7
       child(j).mat=mutation(out.child1);
       child(j+1).mat=mutation(out.child2);
    else
    child(j).mat=out.child1;
    child(j+1).mat=out.child2;
    end
%     child(population+j).mat=parent(j).mat;
%     child(population+j+1).mat=parent(j+1).mat;
end
childi = child;
for i= 1:(new_generation) % check for second iteration for child 3 and 4

    max_sum_child = max(sum(child(i).mat));
    min_sum_child = min(sum(child(i).mat));
    while max_sum_child ~= 1 || min_sum_child == 0
         child(i).mat =  horicorrection(child(i).mat);
%          l=l+1;
%          disp(max(sum(child(i).mat')))
         max_sum_child = max(sum(child(i).mat));
         min_sum_child = min(sum(child(i).mat));
    end
        max_sum_child = max(sum(child(i).mat'));
    while max_sum_child > 3
         child(i).mat =  correction(child(i).mat);
%          l=l+1;
%          disp(max(sum(child(i).mat')))
          max_sum_child = max(sum(child(i).mat'));
    end
    
end

thin = [];
thick = [];

for i = 1:new_generation
    thin(i) = max(sum(child(i).mat));
    thinner(i) = min(sum(child(i).mat));
    thick(i) = max(sum(child(i).mat'));
    if thinner(i) == 0 || thin(i) ~= 1 || thick(i) > 5
    end
end

for i = 1 : new_generation
    child(i).tour = tournum(child(i).mat,C,P);
    mat11 = [];
    temptw=[];
    temp=[];
    for k = 1: C
    mat11 = cell2mat(struct2cell(child(i).tour(k)));
    for g1 = 1:length(mat11)
        mat11 = [mat11 mat11(g1)+length(combpos)/2];
    end
    temp = routetoxy(mat11,pos,P);
    for g = 1:(length(mat11))
%         if g <= (length(mat11)/2)
            temptw(g,1) = combtimespace(mat11(g));
%         else
%             temptw(g,1) = combtimespace(mat11(g)+length(mat11)) ; 
%         end
    end
%     child(i).xy(:,1,k)=temp(:,1);
%     child(i).xy(:,2,k)=temp(:,2);
if isempty(temp)
   opttourlength(k) = 0;
else
   out = ACOTW(cars(k,:),temp,temptw);
   opttourlength(k)=out.Cost;
   child(i).tour(k).xy=out.Tour;
end
    end
    child(i).tourlength=sum(opttourlength);
end
new_population = table2struct(sortrows([struct2table(parent); struct2table(child)],'tourlength'));
parent = new_population(1:population);

thin = [];
thick = [];

for i = 1:new_generation
    thin(i) = max(sum(child(i).mat));
    thinner(i) = min(sum(child(i).mat));
    thick(i) = max(sum(child(i).mat'));
    if thinner(i) == 0 || thin(i) ~= 1 || thick(i) > 5
    end
end
besttourlength(iter)=parent(1).tourlength;
disp(iter)
ccheck = 0;
for h =1
    for j = 1: population
        if parent(h).mat == parent(j).mat
            ccheck = ccheck+1;
        end
    end
end

if ccheck == population
    
end
    
end
%%
% for i = 1: length(parent(1).tour)
%     plot(parent(1).tour(i).xy(1,:),parent(1).tour(i).xy(2,:),'o-');
%     hold on
% end
% drawnow
% hold off
% 
% x=1:iter;
% figure()
% plot(x,besttourlength)
%%


%%  calculating route for ACO
function col = tournum(child,C,P)
    idx = find(child == 1);
    idx = idx';

    for i=1:C
    car(i).tour = [];
    end
    for l=1:P
    r = rem(idx(l),C);
    if r==0
        r=C;
    end
    car(r).tour= [car(r).tour idx(l)];
%     car(r,length(car(r,:))+1)=idx(l);
    end
    for k = 1:C
        [row,col(k).tour]=ind2sub(size(child),car(k).tour);
    end
  
 end
%%
function out = correction(child)
    s1=sum(child');
    maxk=max(s1);
    maxpass=find(s1 == maxk);
    idn=find(child(maxpass(1),:) == 1); % put randomperm for multiple maximum 
    pk=randperm(length(idn));
    child(maxpass(1),idn(pk(1)))=0;
    mink=min(s1);
    minpass=find(s1 == mink);
    minp=randperm(length(minpass));
    child(minpass(minp(1)), idn(pk(1)))=1;
    out = child;
end
function out = horicorrection(child)
    s1=sum(child);
    maxk=max(s1);
    while maxk >= 2
    maxpass=find(s1 == maxk);
    for i=1:length(maxpass)
         idn=find(child(:,maxpass(i)) == 1); % put randomperm for multiple maximum 
         pk=randperm(length(idn));
         child(idn(pk(1)),maxpass(i))=0;
         s1=sum(child);
         maxk=max(s1);
    end
    end
    mink=min(s1);
    while mink == 0
    zeropass=find(s1 == 0);
    for i=1:length(zeropass)
         idn=find(child(:,zeropass(i)) == 0); % put randomperm for multiple maximum 
         pk=randperm(length(idn));
         child(idn(pk(1)),zeropass(i))=1;
         s1=sum(child);
         mink=min(s1);
    end
%     minpass=find(s1 == mink);
%     minp=randperm(length(minpass));
%     child(minpass(minp(1)), idn(pk(1)))=1;
    end
    out = child;
end

function out = vertisplit1(parent1,parent2)
S=size(parent1);
k1 = randi(length(parent1));
k2 = randi(length(parent1));
k = sort([k1 k2]);
child1=zeros(S);
child2=zeros(S);
child1(:,1:k(1))=parent1(:,1:k(1));
child1(:,k(1)+1:k(2))= parent2(:,k(1)+1:k(2));
child1(:,k(2)+1:S(2))= parent1(:,k(2)+1:S(2));

child2(:,1:k(1))=parent2(:,1:k(1));
child2(:,k(1)+1:k(2))= parent1(:,k(1)+1:k(2));
child2(:,k(2)+1:S(2))= parent2(:,k(2)+1:S(2));

out.child1=child1;
out.child2=child2;
end

function out = vertisplit2(parent1,parent2)

S=size(parent1);
k = randi(length(parent1));
child1=zeros(S);
child2=zeros(S);
child1(:,1:k)=parent1(:,1:k);
child1(:,k+1:S(2))= parent2(:,k+1:S(2));
child2(:,1:k)=parent2(:,1:k);
child2(:,k+1:S(2))= parent1(:,k+1:S(2));
out.child1=child1;
out.child2=child2;

end
function out = horisplit(parent1,parent2)
S=size(parent1);
k = randi((S(1)));
if k == S(1)
    k=k-1;
end
child1=zeros(S);
child2=zeros(S);
child1(1:k,:)=parent1(1:k,:);
child1(k+1:S(1),:)= parent2(k+1:S(1),:);
child2(1:k,:)=parent2(1:k,:);
child2(k+1:S(1),:)= parent1(k+1:S(1),:);
out.child1=child1;
out.child2=child2;
end
function out = mutation(parent1)
% j=size(parent1);
% parent1=zeros(j);
% for i=1:j(2)
%      k=randi(j(1));
%      parent1(k,i)=1;
% end
%  out=parent1;
j= size(parent1);
for m=1:3
    k=randi(j(2));
    maxpass=find(parent1(:,k) == 1);
    parent1(maxpass,k)=0;
    l=randi(j(1));
    parent1(l,k)=1;
end
out=parent1;

end
%% ACO function

function out = Antcolony(pos)
x = pos(:,1);
y = pos(:,2);
n = numel(x);
N=n/2;
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
it = 30;      
Ant  = val.n;      % check for number of `ants
Q     = 100;
tau0=10^-6;	
alpha=1;        
beta=5;         
rho=0.5;       
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
        nonpool = N+1:2*N;
        ant(k).Tour=randi([1 N]); 
        pool = [1:N ant(k).Tour(end)+N];
        nonpool(nonpool==ant(k).Tour(end)+N)=[];
        for l=2:2*N           
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
    

    out.cost=BestSol.Cost; 
     Bestpath(it)=BestSol.Cost;
%      disp(['Iteration ' num2str(it) ': Best path length = ' num2str(Bestpath(it))]);
%      figure(1);
      tour=[BestSol.Tour];
%      plot(val.x(tour),val.y(tour),'-s','MarkerSize',10,...
%      'MarkerEdgeColor','red',...
%      'MarkerFaceColor',[1 .6 .6])
    out.tour(1,:)= val.x(tour);
    out.tour(2,:)=val.y(tour);
%     
%     xlabel('x');
%     ylabel('y');
%     grid on;
%     drawnow;
    
end
end

%%
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
function out=routetoxy(route,pos,P)
k =route(1:length(route)/2);
if isempty(k)
out=[];
else
for i= 1:length(k)
    xp(i)=pos(k(i),1);
    yp(i)=pos(k(i),2);
    xd(i)=pos(k(i)+P,1);
    yd(i)=pos(k(i)+P,2);
end

x = [xp xd];
y = [yp yd];
out(:,1)=x';
out(:,2)=y';
end
end