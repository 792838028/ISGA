function [gBestScore,gBest,cg] = ISGA (Np, Max_iter, Xmin, Xmax, dim, fobj)
t =1;
Xmax = Xmax.*ones(1,dim);
Xmin = Xmin.*ones(1,dim); 
Pos = initialization(Np,dim,Xmax,Xmin);
L_low=20*ones(Np,dim);
lw=65*ones(Np,dim);
Vel = zeros(Np,dim);
gBest = zeros(1,dim);
gBestScore = inf;
cg = zeros(1,Max_iter);
Elite_pool=[];
for i = 1:Np
    fitness(i) = fobj(Pos(i,:));
end
[~,index] = sort(fitness);
gBest = Pos(index(1),:);
gBestScore=fitness(index(1));
cg(1) = gBestScore;
Elite_pool(1,:)=gBest;
Elite_pool(2,:)=Pos(index(2),:);
Elite_pool(3,:)=Pos(index(3),:);
beta =3/2; 
while t < Max_iter
    [index_FDB]=fitnessDistanceBalance(Pos, fitness);
    coe = (4*(t/Max_iter))/exp(4*(t/Max_iter));
         if (t/Max_iter)<=0.85
             fi = (rand()+((t/Max_iter)))*pi;
         else
             fi = (1+(t/Max_iter))*pi;    
         end

    for i = 1:Np
        for j=1:dim
        acc = ((gBest(j) - Pos(i,j)) - 1.29*Vel(i,j).^2*sin(fi))*10^-2;
        Vel(i,j) = coe*Vel(i,j) + acc;
        end
    end

   [~,index] = sort(fitness);
   for i = 1: Np
       New_Pos(i,:) = Pos(index(i),:);
       New_Vel(i,:) = Vel(index(i),:);
   end
   Pos = New_Pos;
   Vel = New_Vel;
    a = 4*rand() - 2;
    b = 3*rand() -1.5;
    c = 2*rand() - 1;

    for i =1:Np
        aa(i,:) = Pos(i,:).*fitness(i);
        bb(i) = Np*fitness(i);
    end
    Xc = sum(aa)/sum(bb);
    Pos = Pos+Vel;

   
%% Exploration phase
    if fi < pi
            for i = 1:Np
                for j=1:dim
                if i<=1/5*Np
                    Pos(i,j) =Pos(i,j)+a*(Pos(index_FDB,j)-Pos(i,j))+Vel(i,j);
                elseif  1/5*Np< i && i < 4/5*Np
                    Pos(i,j) =Pos(i,j)+a*(Pos(index_FDB,j)-Pos(i,j))+b*(Xc(j)-Pos(i,j))-c*(Pos(Np,j)+Pos(i,j))+Vel(i,j);  
                else
                    Pos(i,j) =Pos(i,j)+a*(Pos(index_FDB,j)-Pos(i,j))+b*(Xc(j)-Pos(i,j))+Vel(i,j);
                end
                Flag4ub(i,j)=Pos(i,j)>Xmax(j);
                Flag4lb(i,j)=Pos(i,j)<Xmin(j);
                Pos(i,j)=(Pos(i,j).*(~(Flag4ub(i,j)+Flag4lb(i,j))))+Xmax(j).*Flag4ub(i,j)+Xmin(j).*Flag4lb(i,j);
                end
                fitness_1(i) = fobj(Pos(i,:));
                   t=t+1;
                if fitness_1(i)<gBestScore
                   gBest = Pos(i,:);
                   gBestScore = fitness_1(i);

                end
                   cg(t)=gBestScore;
            end
%% Lead goose rotation mechanism
            [~,index_1] = sort(fitness_1);
              for i = 1:round(1/5*Np)
                for j=1:dim
                k1=randperm(3,1);
                k2=randperm(3,1);
                k3=randperm(5,1);
                k4=randperm(round(1/5*Np),1);
                k5=rand;
                Pos0(i,j) = (1-k5)*Pos(index_1(k1),j)+k5*(Elite_pool(k2,j)+Pos(index_1(1),j))/2+(2*rand-0.5)*(Pos(index_1(1),j)-(Pos(index_1(k3),j)+Pos(index_1(i),j))/2);
                Flag4ub(i,j)=Pos0(i,j)>Xmax(j);
                Flag4lb(i,j)=Pos0(i,j)<Xmin(j);
                Pos0(i,j)=(Pos0(i,j).*(~(Flag4ub(i,j)+Flag4lb(i,j))))+Xmax(j).*Flag4ub(i,j)+Xmin(j).*Flag4lb(i,j);
                end
                fitness_2(i) = fobj(Pos0(i,:));
                if fitness_2(i)<fitness_1(index_1(i))
                   Pos(index_1(i),:)=Pos0(i,:);
                   fitness_1(index_1(i))=fitness_2(i);
                end
                t=t+1;
                if fitness_1(index_1(i))<gBestScore
                   gBest = Pos(index_1(i),:);
                   gBestScore = fitness_1(index_1(i));

               end
                cg(t)=gBestScore;
              end
              fitness=fitness_1;
%% Exploitation phase
    else
        if rand>0.5
            for i = 1:Np
                for j=1:dim
                Pos(i,j) = Pos(i,j) + (Pos(i,j) - gBest(j))*rand;
                end
            end
        else
            for i = 1:Np
                for j=1:dim
                Pos(i,j) = gBest(j) + (Pos(i,j) - gBest(j))*rand*Brownian(1);
                end
            end
        end
%% Honk-guiding mechanism
        sum1=zeros(1,dim);
        for i=1:Np
            sum1=sum1+Pos(i,:);
        end
        X_centroid=sum1/Np;
        for i = 1:Np
                for j=1:dim
                r(i,j)=50*(gBest(j)-Pos(i,j))/(Xmax(j)-Xmin(j));
                if r(i,j) < 0.28
                   r(i,j) = 0.28;
                end
                Lp(i,j)=lw(i,j)-20*log10(r(i,j))-11;
                L(i,j)=(Lp(i,j)-L_low(i,j))/(lw(i,j)-L_low(i,j));
                end
                Pos(i,:) = Pos(i,:) + (1-L(i,:)).*(gBest - Pos(i,:)).*rand(1,dim).*Brownian(dim)+(1.5*rand(1,dim)-0.5).*(X_centroid-Pos(i,:));%跟随强壮个体试图跳出局部解
                for j=1:dim
                Flag4ub(i,j)=Pos(i,j)>Xmax(j);
                Flag4lb(i,j)=Pos(i,j)<Xmin(j);
                Pos(i,j)=(Pos(i,j).*(~(Flag4ub(i,j)+Flag4lb(i,j))))+Xmax(j).*Flag4ub(i,j)+Xmin(j).*Flag4lb(i,j);
                end
        end   
    end
%% Outlier boundary

    [~,index]=sort(fitness);
    [fmax,B]=max( fitness );
    [~,A]=min( fitness );
    fitness_avg=sum(fitness,'all')/Np;
    low_limit = ceil(numel(index)*0.2);
    up_limit = floor(numel(index)*1);
    mid_index = low_limit:up_limit;
    rand_mid_index = mid_index(randperm(length(mid_index)));
    num_selected = max(1, round(length(Np)*0.5));
    b2 = index(rand_mid_index(1:num_selected));
    for q =  1  : length(b2)
        if( fitness( index( b2(q) ) )>(fitness_avg) )
            Pos(index( b2(q)),:)=Pos(A,:)+(randn(1,dim)).*(Pos(A,:)-Pos(index(b2(q)),:));
        else
            Pos(index(b2(q)),:) =Pos(index(b2(q)),:)-(2*rand(1,dim)-1).*(Pos(index(b2(q)),:)-Pos(B,:))/(fmax-fitness(index(b2(q))))+levy(dim,beta);
        end
                for j=1:dim
                Flag4ub(index( b2(q)),j)=Pos(index( b2(q)),j)>Xmax(j);
                Flag4lb(index( b2(q)),j)=Pos(index( b2(q)),j)<Xmin(j);
                Pos(index( b2(q)),j)=(Pos(index( b2(q)),j).*(~(Flag4ub(index( b2(q)),j)+Flag4lb(index( b2(q)),j))))+Xmax(j).*Flag4ub(index( b2(q)),j)+Xmin(j).*Flag4lb(index( b2(q)),j);
                end
        if fi<pi
        fitness(index( b2(q))) = fobj(Pos(index( b2(q)),:));
        t=t+1;
        if fitness(index( b2(q)))<gBestScore
                   gBest = Pos(index( b2(q)),:);
                   gBestScore = fitness(index( b2(q)));
        end
                cg(t)=gBestScore;
        end
    end
%% 
    if fi>=pi
        for i=1:Np
        fitness(i) = fobj(Pos(i,:));
        t=t+1;
        if fitness(i)<gBestScore
                   gBest = Pos(i,:);
                   gBestScore = fitness(i);
        end
                cg(t)=gBestScore;
        end
    end
    [~,index]=sort(fitness);
    Elite_pool(1,:)=gBest;
    Elite_pool(2,:)=Pos(index(2),:);
    Elite_pool(3,:)=Pos(index(3),:);
end
if t>=Max_iter
cg=cg(end-(Max_iter-1):end);
end
end
%%
function o = Brownian(dim)
    T = 1;
    r = T/dim;
    dw = sqrt(r)*randn(1,dim);
    o = cumsum(dw);
end
function Levy_step = levy(dimension,beta)
alpha_u = (gamma(1+beta)*sin(pi*beta/2)/(gamma(((1+beta)/2)*beta*2^((beta-1)/2))))^(1/beta);
alpha_v = 1;
u=normrnd(0,alpha_u^2,[1 dimension]);
v=normrnd(0,alpha_v^2,[1 dimension]);
Levy_step = 0.01.*u./(abs(v).^(1/beta));
end