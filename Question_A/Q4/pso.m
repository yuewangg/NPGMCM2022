function [gbest,record,records]=pso(obj,D)
%% ������ʼ��
c1=1.4;
c2=1.4;
maxgen=25;             % ��������
sizepop=100;             % ��Ⱥ��ģ
Vmax=0.05;
Vmin=-0.05;
popmax=0.125;
popmin=-0.125;
w=0.8; 

VisDis=0.25;
rand('seed',1); 

%% �������������ٶ�
for i=1:sizepop
    pop(i,:)=(popmax-popmin)*rand(1,D)+popmin;
    V(i,:)=(Vmax-Vmin)*rand(1,D)+Vmin;
    fitness(i)=obj(pop(i,:));
end

%% ���弫ֵ��Ⱥ�弫ֵ
[bestfitness, bestindex]=min(fitness);
gbest=pop(bestindex,:);   %ȫ�����,��ȡ��������Ӧ��ֵ�����Ӷ���Ϊȫ���������
pbest=pop;                %������� 
fitnesszbest=bestfitness; %ȫ�������Ӧ��ֵ

%% ����
for i=1:maxgen 
    i
    for j=1:sizepop
        % �ٶȸ���
        V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));
        % �ٶȱ߽紦��
        V(j,:)=min(max(V(j,:),Vmin),Vmax); 
        
        % ��Ⱥ����
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,:)=min(max(pop(j,:),popmin),popmax);
        
        %��Ӧ��ֵ����
        fitness(j)=obj(pop(j,:));
    end
    
    for j=1:sizepop
        fitnesstmp=fitness(j);
        for j2=1:sizepop
            % �������Ÿ���
            if fitness(j2)<fitnesstmp && norm(pop(j,:)-pop(j2,:))<VisDis
                pbest(j,:)=pop(j2,:);
                fitnesstmp=fitness(j2);
            end
        end
        % Ⱥ�����Ÿ���
        if fitness(j)<fitnesszbest
            gbest=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    record(1,i)=fitnesszbest; 
    records(i,:)=gbest;
 end

% %% ��ͼ
% figure;
% plot(record);
% xlabel('��������');
% ylabel('GDOP'); 

 