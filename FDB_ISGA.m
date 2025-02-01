function [gBestScore,gBest,cg] = FDB_ISGA(Np, Max_iter, Xmin, Xmax, dim, fobj)
Np=50;
t =1;
Xmax = Xmax.*ones(1,dim);
Xmin = Xmin.*ones(1,dim); 
Pos = initialization(Np,dim,Xmax,Xmin);

L_low=20*ones(Np,dim);
lw=65*ones(Np,dim);
cs=round(Max_iter*0.05);
Vel = zeros(Np,dim);

gBest = zeros(1,dim);
gBestScore = inf;
cg = zeros(1,Max_iter);
Elite_pool=[];
for i = 1:Np
    fitness(i) = fobj(Pos(i,:));
end
[~,index] = sort(fitness);%根据适应度值排序
gBest = Pos(index(1),:);
gBestScore=fitness(index(1));
cg(1) = gBestScore;
%精英池
Elite_pool(1,:)=gBest;
Elite_pool(2,:)=Pos(index(2),:);
Elite_pool(3,:)=Pos(index(3),:);

%莱维飞行参数
beta =3/2; 
while t < Max_iter
    [index_FDB]=fitnessDistanceBalance(Pos, fitness);
    coe = (4*(t/Max_iter))/exp(4*(t/Max_iter));%速度的权重因子
         % if (t/Max_iter)<=0.15
         %     fi = (t/Max_iter)*pi;
         % elseif (t/Max_iter)>0.15 && (t/Max_iter)<0.85
         %     fi = (rand()+((t/Max_iter)))*pi;%雪雁人字形夹角θ
         % else
         %     fi = 2*(t/Max_iter)*pi;    
         % end 
         if (t/Max_iter)<=0.85
             fi = (rand()+((t/Max_iter)))*pi;%雪雁人字形夹角θ
         else
             fi = (1+(t/Max_iter))*pi;    
         end

     % fi = (rand()+(0.05+0.9*(t/Max_iter)))*pi;%雪雁人字形夹角θ
    for i = 1:Np
        for j=1:dim
        % if rand<0.5
        % acc = ((Pos(index_FDB,j) - Pos(i,j)) - 1.29*Vel(i,j).^2*sin(fi))*10^-2;%加速度：最佳位置和个体位置差再减去空气阻力
        % else
        acc = ((gBest(j) - Pos(i,j)) - 1.29*Vel(i,j).^2*sin(fi))*10^-2;%加速度：最佳位置和个体位置差再减去空气阻力
        % end
        Vel(i,j) = coe*Vel(i,j) + acc;%速度更新
        end
    end

   [~,index] = sort(fitness);%根据适应度值排序
   for i = 1: Np
       New_Pos(i,:) = Pos(index(i),:);%按照顺序对种群和对应的速度进行排序
       New_Vel(i,:) = Vel(index(i),:);
   end
   Pos = New_Pos;%排完放回去
   Vel = New_Vel;
    a = 4*rand() - 2;
    b = 3*rand() -1.5;
    c = 2*rand() - 1;

    for i =1:Np
        aa(i,:) = Pos(i,:).*fitness(i);
        bb(i) = Np*fitness(i);
    end
    Xc = sum(aa)/sum(bb);
    Pos = Pos + Vel;

   
%% 全局搜索
    if fi < pi%夹角小于π
            for i = 1:Np
                for j=1:dim
                if i<=1/5*Np
                % if rand<0.5
                    Pos(i,j) = Pos(i,j)  + a*(Pos(index_FDB,j) - Pos(i,j)) + Vel(i,j);
                % else
                %     Pos(i,j) = Pos(i,j)  + a*(gBest(j) - Pos(i,j)) + Vel(i,j);%最佳个体位置更新:只负责寻找最佳位置    
                % end
                elseif  1/5*Np< i && i < 4/5*Np
                    % if rand<0.5
                    Pos(i,j) = Pos(i,j)  + a*(Pos(index_FDB,j) - Pos(i,j)) + b*(Xc(j) -Pos(i,j)) - c*(Pos(Np,j) + Pos(i,j))  + Vel(i,j);  
                    % else
                    % Pos(i,j) = Pos(i,j)  + a*(gBest(j) - Pos(i,j)) + b*(Xc(j) -Pos(i,j)) - c*(Pos(Np,j) + Pos(i,j))  + Vel(i,j);%其他个体：受到最佳位置、群体平均位置、以及最差位置的影响  
                    % end
                else
                    % if rand<0.5
                    Pos(i,j) = Pos(i,j)  + a*(Pos(index_FDB,j) - Pos(i,j)) + b*(Xc(j) -Pos(i,j)) + Vel(i,j);
                    % else
                    % Pos(i,j) = Pos(i,j)  + a*(gBest(j) - Pos(i,j)) + b*(Xc(j) -Pos(i,j)) + Vel(i,j);%最差个体位置：受最佳位置和群体平均位置的影响  
                    % end
                end
                Flag4ub(i,j)=Pos(i,j)>Xmax(j);%边界检查
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
%% 头雁竞争机制（头雁池）
            [~,index_1] = sort(fitness_1);%根据适应度值排序
            % pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',1);
            % rr = random(pd,length(1:round(1/5*Np)),1);    %Generate 1 Cauchy random number
            % if fitness_1(index_1(1)) >= gBestScore
              for i = 1:round(1/5*Np)
                for j=1:dim
                k1=randperm(3,1);
                k2=randperm(3,1);
                k3=randperm(5,1);
                k4=randperm(round(1/5*Np),1);
                k5=rand;
                Pos0(i,j) = (1-k5)*Pos(index_1(k1),j)+k5*(Elite_pool(k2,j)+Pos(index_1(1),j))/2+(2*rand-0.5)*(Pos(index_1(1),j)-(Pos(index_1(k3),j)+Pos(index_1(i),j))/2);
                Flag4ub(i,j)=Pos0(i,j)>Xmax(j);%边界检查
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
            % end
%% 局部搜索
    else%夹角大于π即切换一字型队列局部搜索   
        if rand>0.5
            for i = 1:Np
                for j=1:dim
                % if rand<0.5
                % Pos(i,j) = Pos(i,j) + (Pos(i,j) - Pos(index_FDB,j))*rand;%跟随强壮个体试图跳出局部解
                % else
                Pos(i,j) = Pos(i,j) + (Pos(i,j) - gBest(j))*rand;%跟随强壮个体试图跳出局部解
                % end
                end
            end
        else
            for i = 1:Np
                for j=1:dim
                % if rand<0.5
                % Pos(i,j) = Pos(index_FDB,j) + (Pos(i,j) - Pos(index_FDB,j))*rand*Brownian(1);%布朗运动跳出局部解
                % else
                Pos(i,j) = gBest(j) + (Pos(i,j) - gBest(j))*rand*Brownian(1);%布朗运动跳出局部解
                % end
                end
            end
        end
%% 叫声引导机制
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
                Flag4ub(i,j)=Pos(i,j)>Xmax(j);%边界检查
                Flag4lb(i,j)=Pos(i,j)<Xmin(j);
                Pos(i,j)=(Pos(i,j).*(~(Flag4ub(i,j)+Flag4lb(i,j))))+Xmax(j).*Flag4ub(i,j)+Xmin(j).*Flag4lb(i,j);
                end

        end   
    end
%% 警戒机制(换成防止变成离群之雁)警惕值讲成和群体之间的距离边界，超过就是离群

    [~,index]=sort(fitness);
    [fmax,B]=max( fitness );
    [fmin,A]=min( fitness );
    fitness_avg=sum(fitness,'all')/Np;
    % 避免小数索引
    low_limit = ceil(numel(index)*0.2); % 20%位置的索引
    up_limit = floor(numel(index)*1); % 80%位置的索引
    % 选取20%到80%之间的索引
    mid_index = low_limit:up_limit;
    % 随机排列这个范围内的索引
    rand_mid_index = mid_index(randperm(length(mid_index)));
    % 选择这些索引中的20%个
    num_selected = max(1, round(length(Np)*0.5)); % 至少选择一个元素
    b2 = index(rand_mid_index(1:num_selected));
    for q =  1  : length(b2)      % Equation (5)
        if( fitness( index( b2(q) ) )>(fitness_avg) )
            Pos(index( b2(q)),:)=Pos(A,:)+(randn(1,dim)).*(Pos(A,:)-Pos(index(b2(q)),:));
        else
            Pos(index(b2(q)),:) =Pos(index(b2(q)),:)-(2*rand(1,dim)-1).*(Pos(index(b2(q)),:)-Pos(B,:))/(fmax-fitness(index(b2(q))))+levy(dim,beta);
        end
                for j=1:dim
                Flag4ub(index( b2(q)),j)=Pos(index( b2(q)),j)>Xmax(j);%边界检查
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

%% 检查迭代
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
    %% Update the elite pool
    [~,index]=sort(fitness);
    Elite_pool(1,:)=gBest;
    Elite_pool(2,:)=Pos(index(2),:);
    Elite_pool(3,:)=Pos(index(3),:);
end
if t>=Max_iter
cg=cg(end-(Max_iter-1):end);
end
end



function o = Brownian(dim)
    T = 1;
    r = T/dim;
    dw = sqrt(r)*randn(1,dim);
    o = cumsum(dw);
end
function Levy_step = levy(dimension,beta)
% nvar : 求解变量的个数
% beta = 1.5  常数
% 常数X
% beta = 3/2;
% 方差alpha_u
alpha_u = (gamma(1+beta)*sin(pi*beta/2)/(gamma(((1+beta)/2)*beta*2^((beta-1)/2))))^(1/beta);
% 方差alpha_v
alpha_v = 1;
% u ,v服从正态分布
u=normrnd(0,alpha_u^2,[1 dimension]);
v=normrnd(0,alpha_v^2,[1 dimension]);
Levy_step = 0.01.*u./(abs(v).^(1/beta));

end