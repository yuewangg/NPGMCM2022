function [gbest,record,records]=pso(obj,D)
%% 参数初始化
c1=1.4;
c2=1.4;
maxgen=25;             % 迭代次数
sizepop=100;             % 种群规模
Vmax=0.05;
Vmin=-0.05;
popmax=0.125;
popmin=-0.125;
w=0.8; 

VisDis=0.25;
rand('seed',1); 

%% 产生初试粒子速度
for i=1:sizepop
    pop(i,:)=(popmax-popmin)*rand(1,D)+popmin;
    V(i,:)=(Vmax-Vmin)*rand(1,D)+Vmin;
    fitness(i)=obj(pop(i,:));
end

%% 个体极值和群体极值
[bestfitness, bestindex]=min(fitness);
gbest=pop(bestindex,:);   %全局最佳,将取得最优适应度值的粒子定义为全局最佳粒子
pbest=pop;                %个体最佳 
fitnesszbest=bestfitness; %全局最佳适应度值

%% 迭代
for i=1:maxgen 
    i
    for j=1:sizepop
        % 速度更新
        V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));
        % 速度边界处理
        V(j,:)=min(max(V(j,:),Vmin),Vmax); 
        
        % 种群更新
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,:)=min(max(pop(j,:),popmin),popmax);
        
        %适应度值更新
        fitness(j)=obj(pop(j,:));
    end
    
    for j=1:sizepop
        fitnesstmp=fitness(j);
        for j2=1:sizepop
            % 个体最优更新
            if fitness(j2)<fitnesstmp && norm(pop(j,:)-pop(j2,:))<VisDis
                pbest(j,:)=pop(j2,:);
                fitnesstmp=fitness(j2);
            end
        end
        % 群体最优更新
        if fitness(j)<fitnesszbest
            gbest=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    record(1,i)=fitnesszbest; 
    records(i,:)=gbest;
 end

% %% 绘图
% figure;
% plot(record);
% xlabel('迭代次数');
% ylabel('GDOP'); 

 