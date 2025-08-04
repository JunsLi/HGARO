    %%% Artificial Rabbits Optimization (ARO) for 23 functions %%%
% function [BestX,BestF,HisBestF,Ave,Std]=HARO(F_index,MaxIt,nPop,fobj)
function [BestX,BestF,HisBestF,DNASet]=HARO(nPop,MaxIt,Low,Up,Dim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FunIndex: Index of function.                       %
    % MaxIt: Maximum number of iterations.               %
    % PopSize: Size of population.                       %
    % PopPos: Position of rabbit population.             %
    % PopFit: Fitness of population.                     %
    % Dim: Dimensionality of prloblem.                   %
    % BestX: Best solution found so far.                 %
    % BestF: Best fitness corresponding to BestX.        %
    % HisBestF: History best fitness over iterations.    %
    % Low: Low bound of search space.                    %
    % Up: Up bound of search space.                      %
    % R: Running operator.                               %
    % L:Running length.                                  %
    % A: Energy factor.                                  %
    % H: Hiding parameter.                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Low,Up,Dim]=FunRange(F_index);
% PopPos=zeros(nPop,Dim);
% PopFit=zeros(nPop,1);
 Up = Up.*ones(1,Dim);
   Low = Low.*ones(1,Dim); 

PopPos=initialization(nPop,Dim,Up,Low);
PopFit=zeros(nPop,1);
% X1=PopPos(:,1:Dim);
% X2=PopPos(:,2*Dim+1:3*Dim);
% PopFit=zeros(size(PopPos,1),1);
% Run_no=30;
% mm=100;
% for mn= 1:mm
%  for irun=1:Run_no
% Best_XSF = zeros(1,Dim);  % Determine the vale of SF Best Fitness
NewDNA=[];

% for i=1:nPop
%     PopPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
% %     PopFit(i)=BenFunctions(PopPos(i,:),F_index,Dim);
% end
% %================！！！！！！！！！！！！！！！！！饥饿游戏初始参数
Destination_fitness=inf;%change this to -inf for maximization problems
Worstest_fitness=-inf;
% AllFitness = inf*ones(nPop,1);%record the fitness of all positions
% VC1=ones(nPop,1);
hungry0 = zeros(1,size(PopPos,1));
tempPosition0=zeros(nPop,Dim);

weight3 = ones(nPop,Dim);%hungry weight of each position 
weight4=ones(nPop,Dim);%hungry weight of each position 

hungry = zeros(1,size(PopPos,1));
tempPosition=zeros(nPop,Dim); 
count=0;
% 
% %===============================================！！！！！！！！！
BestF=inf;
BestX=[];
% for i=1:nPop
%     if PopFit(i)<=BestF
%         BestF=PopFit(i);
%         BestX=PopPos(i,:);
%     end
% end

HisBestF=zeros(MaxIt,1);

standard=(-42.9)*(Dim-1);
tic
NewDNA=[];
DNASet=[];
  tic
for It=1:MaxIt
  
     VC2 = 0.03; %The variable of variation control 
    sumHungry = 0;%record the sum of each hungry
    
    Direct1=zeros(nPop,Dim);
    Direct2=zeros(nPop,Dim);
    theta=2*(1-It/MaxIt);
    %%===================s==============================
    %%%%%%%%%%%%%%%%原算法更新前符合要去的序列个数
PopPos001=PopPos;
     PopPos001=unique(PopPos001,'row','stable');
    NewDNA001=[];
    % %%%%%%%%%%%%%%%%%%%%%%%GC含量
    px= size(PopPos001,1);
index = 1;
temp=[];
for q = 1:px
    if   CTcontent(PopPos001(q,:))~=1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
     PopPos001(temp,:)=[];
end
%=====================
%四C碱基不连续约束
% %==================
index2=1;
temp=[];
pxx = size(PopPos001,1);
for q = 1:pxx
    if   Discontinuous(PopPos001(q,:))==1
        temp(index2) =q;
        index2 = index2 + 1;
    end
end
if length(temp)~=0
     PopPos001(temp,:)=[];
end
   %%%%% %添加末端约束（左边不能是TC,右边不能是CT(末端是T，此末端是C)）
% if size(NewDNA1,1)==0
%     NewDNA1=[];
% else
% x7 = size(PopPos001,1);
% index = 1;
% temp=[];
% % XX3=NewDNA1(:,2*Dim+1:3*Dim);
% for q = 1:x7
%     if  Endlimit(PopPos001(q,:))==1
%         temp(index) =q;
%         index = index + 1;
%     end
% end
% if length(temp)~=0
%    PopPos001(temp,:)=[];
% end   
%%%%%%%%%%%======================
% % end
%   %=================%Add 扭应力约束===================
% if size(NewDNA001,1)==0
%     NewDNA001=[];
% else
x5 = size(PopPos001,1);
index = 1;
temp=[];
% XX3=NewDNA1(:,2*Dim+1:3*Dim);
for q = 1:x5
    if  CTC_TCTcontent(PopPos001(q,:))==1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
   PopPos001(temp,:)=[];
end   
    %=====================================================
     %   %I-motif free constraint 
%  if size(NewDNA1,1)==0
%     NewDNA1=[];
% else
x6 = size(PopPos001,1);
index = 1;
temp=[];
% Xxn=NewDNA1(:,1:Dim);
for q = 1:x6
    if  Imotif(PopPos001(q,:))==0
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
   PopPos001(temp,:)=[];
end  
    %=======================================================================
%            PopPos001=unique(PopPos001,'row','stable');
              PopFit001=[];
         for i=1:size(PopPos001,1)
        PopFit001(i)=Eint(PopPos001(i,:));
         end
     for i=1:size(PopPos001,1)
         
%          if PopFit001(i)>standard
        if PopFit001(i)<=standard
            NewDNA001=[PopPos001(i,:);NewDNA001];
        end
     end  
  %%%%%%%%%%%%%%%  更新前，约束筛选，最终输出        %%%%%%%%%%%%
    %==========================================================
    newpop=size(PopPos,1);
    for i=1: newpop
        L=(exp(1)-exp(((It-1)/MaxIt)^2))*(sin(2*pi*rand)); %Eq.(3)
        rd=ceil(rand*(Dim));
        Direct1(i,randperm(Dim,rd))=1;
        c=Direct1(i,:); %Eq.(4)
        R=L.*c; %Eq.(2)
        A=2*log(1/rand)*theta;%Eq.(15)
        if A>1
            K=[1:i-1 i+1:  newpop];
            UV=size(PopPos,1);
            RandInd=K(randi([1 UV-1]));
%             newPopPos=PopPos(RandInd,:)+R.*( PopPos(i,:)-PopPos(RandInd,:))+round(0.5*(0.05+rand))*randn; %Eq.(1) 迂回觅食策略
%             原论文中就这两个策略， 可以加1-2个策略提高寻优性能
              newPopPos(i,:)=PopPos(RandInd,:)+R.*( PopPos(i,:)-PopPos(RandInd,:))...
                +round(0.5*(0.05+rand))*randn; %Eq.(1)迂回觅食
        else
            Direct2(i,ceil(rand*Dim))=1;
            gr=Direct2(i,:); %Eq.(12)
            H=((MaxIt-It+1)/MaxIt)*randn; %Eq.(8)
            b=PopPos(i,:)+H*gr.*PopPos(i,:); %Eq.(13)
            newPopPos(i,:)=PopPos(i,:)+ R.*(rand*b-PopPos(i,:)); %Eq.(11) 随机隐藏策略
        end  
%         newPopPos=SpaceBound(PopPos,Up,Low);
%         newPopFit=BenFunctions(newPopPos,F_index,Dim);

 %=========新位置，原来的2，3，一部分变成了小数， 把这些数变成2或3
         a=[2 3 ];
         for j=1:Dim
         if newPopPos(i,j)<2.5
            newPopPos(i,j) = 2;
         elseif newPopPos(i,j)>2.5
             newPopPos(i,j) = 3;
             else
            newPopPos(i,j) = a(randperm(length(a)));
        end
         end
         
%%%%%%%%%     这些对所以位置进行适应度函数测试， 这样是不对的， 
%%%%%%%%%%%%%     需要对满足一定约束分析适应度大小      
% %             for i=1:size(PopPos,1)
%                 newPopFit(i)=Eint(newPopPos(i,:));
% %             end 
% 
%         if newPopFit(i)<PopFit(i)
%             PopFit(i)=newPopFit(i);
%             PopPos(i,:)=newPopPos(i,:);
%         end
%%%%%错误的计算适应度和位置替换， 导致了部分序列全部替换成适应度最大的序列，
    end
        %%===================s==============================
    %%%%%%%%%%%%%%%%兔子更新后符合要去的序列个数
PopPos002= newPopPos;
     PopPos002=unique(PopPos002,'row','stable');
    NewDNA002=[];
    % %%%%%%%%%%%%%%%%%%%%%%%GC含量
    px= size(PopPos002,1);
index = 1;
temp=[];
for q = 1:px
    if   CTcontent(PopPos002(q,:))~=1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
     PopPos002(temp,:)=[];
end
%=====================
%四C碱基不连续约束
% %==================
index2=1;
temp=[];
pxx = size(PopPos002,1);
for q = 1:pxx
    if   Discontinuous(PopPos002(q,:))==1
        temp(index2) =q;
        index2 = index2 + 1;
    end
end
if length(temp)~=0
     PopPos002(temp,:)=[];
end

   %%%%% %添加末端约束（左边不能是TC,右边不能是CT(末端是T，次末端是C)）
x8 = size(PopPos002,1);
index = 1;
temp=[];
% XX3=NewDNA1(:,2*Dim+1:3*Dim);
% for q = 1:x8
%     if  Endlimit(PopPos002(q,:))==1
%         temp(index) =q;
%         index = index + 1;
%     end
% end
% if length(temp)~=0
%    PopPos002(temp,:)=[];
% end   
%%%%%%%%%%%==================================
%    %%%%%%%%%%  %Add 扭应力约束
x5 = size(PopPos002,1);
index = 1;
temp=[];
% XX3=NewDNA1(:,2*Dim+1:3*Dim);
for q = 1:x5
    if  CTC_TCTcontent(PopPos002(q,:))==1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
   PopPos002(temp,:)=[];
end   
    %========================================================================
     %   %I-motif free constraint 
%  if size(NewDNA1,1)==0
%     NewDNA1=[];
% else
x6 = size(PopPos002,1);
index = 1;
temp=[];
% Xxn=NewDNA1(:,1:Dim);
for q = 1:x6
    if  Imotif(PopPos002(q,:))==0
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
   PopPos002(temp,:)=[];
end   
%  end
    
    
    %=======================================================================
%            PopPos001=unique(PopPos001,'row','stable');
              PopFit002=[];
         for i=1:size(PopPos002,1)
        PopFit002(i)=Eint(PopPos002(i,:));
         end
      
     for i=1:size(PopPos002,1)
%          if PopFit002(i)>standard
        if PopFit002(i)<=standard
            NewDNA002=[PopPos002(i,:);NewDNA002];
        end
     end  
  %%%%%%%%%%%%%%%  更新前，约束筛选，最终输出        %%%%%%%%%%%%
    %==========================================================
     NewDNA003=[NewDNA001;NewDNA002];
     NewDNA003=unique(NewDNA003,'row','stable');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        PopFit003=[];
         for i=1:size(NewDNA003,1)
        PopFit003(i)=Eint(NewDNA003(i,:));
         end


    for i=1:size(NewDNA003,1)
        if PopFit003(i)<BestF
            BestF=PopFit003(i);
            BestX=newPopPos(i,:);%% 关注一下， 这个BestX只是一行的数据
        end
    end
%   %==========================
%，最佳位置不符合要求 ，只是最好的适应度， 但是最好的不符合，，，，，

%=====================
%     a=[2 3 ];
%      for i=1:size(PopPos,1)
%          for j=1:size(PopPos,2)
%          if PopPos(i,j)<2.5
%             PopPos(i,j) = 2;
%          elseif PopPos(i,j)>2.5
%              PopPos(i,j) = 3;
%              else
%             PopPos(i,j) = a(randperm(length(a)));
%         end
%          end
%     end
%     
% %    %%%%%%%%%%%%%饥饿开始
   for i=1:size(newPopPos,1)
        % Check if solutions go outside the search space and bring them back
%         Flag4Up=PopPos(i,:)>Up;
%         Flag4Low=PopPos(i,:)<Low;
%         PopPos(i,:)=(PopPos(i,:).*(~()))+Up.*Flag4Up+Low.*Flag4Low;
        
        
      
        AllFitness(i) = Eint(newPopPos(i,:));
%         AllFitness(i) = F_index(PopPos(i,:));
%         newPopFit=BenFunctions(newPopPos,F_index,Dim);
        
   end

     [PopFit003Sorted,IndexSorted] = sort(PopFit003);
    bestFitness = PopFit003Sorted(1);
    worstFitness = PopFit003Sorted(size(NewDNA003,1));
    
    %Update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=newPopPos(IndexSorted(1),:);
        Destination_fitness = bestFitness;
        count=0;
    end
    
    if worstFitness > Worstest_fitness
        Worstest_fitness = worstFitness;
    end
    
    for i = 1:size(newPopPos,1)
         %calculate the variation control of all positions
         VC1(i) = sech(abs(AllFitness(i)-Destination_fitness));    
         %calculate the hungry of each position
        if Destination_fitness == AllFitness(i)
            hungry(1,i) = 0;
            count = count+1;
            tempPosition(count,:)=newPopPos(i,:);
        else
            temprand = rand();
            c = (AllFitness(i)-Destination_fitness)/(Worstest_fitness-Destination_fitness)*temprand*2*(Up-Low);
            if c<100
                b=100*(1+temprand);
            else
                b=c;
            end   
            hungry(1,i) = hungry(1,i)+ max(b); 
            sumHungry = sumHungry + hungry(1,i);
        end
    end 
    
    %calculate the hungry weight of each position
    for i=1:size(newPopPos,1)
        for j=2:size(newPopPos,2)
                weight3(i,j) = (1-exp(-abs(hungry(1,i)-sumHungry)))*rand()*2;
                if rand()<VC2
                    weight4(i,j) = hungry(1,i)*size(newPopPos,1)/sumHungry*rand();
                else
                    weight4(i,j) = 1;
                end
        end
        
    end
    
%     
%     % Update the Position of search agents
    shrink=2*(1-It/MaxIt); % a decreases linearly fron 2 to 0 
    for i=1:size(newPopPos,1)
        if rand<VC2
            newPopPos(i,:) = newPopPos(i,j)*(1+randn(1));
        else
            A = randi([1,count]);
            for j=1:size(newPopPos,2)
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    newPopPos(i,j) = weight4(i,j)*tempPosition(A,j)+vb*weight3(i,j)*abs(tempPosition(A,j)-PopPos(i,j));
                else
                    newPopPos(i,j) = weight4(i,j)*tempPosition(A,j)-vb*weight3(i,j)*abs(tempPosition(A,j)-PopPos(i,j));
                end
            end
        end
    end
% %  %%%%   %-===============================饥饿结束
%     
      A=[2 3 ];
     for i=1:size(newPopPos,1)
         for j=1:size(newPopPos,2)

         if newPopPos(i,j)<2.5
            newPopPos(i,j) = 2;
         elseif newPopPos(i,j)>2.5
             newPopPos(i,j) = 3;
             else
            newPopPos(i,j) = A(randi(numel(A),1,1));
        end
         end
    end
    %==================================================================================================================================
   %%%%%%%%%%%%%===========================   
    PopPos004=newPopPos;
    NewDNA004=[];
    %%%%%%%%%%%%%%%%%%%去重%%%%%%%%%%%%
      PopPos004=unique(PopPos004,'row','stable');
    % %%%%%%%%%%%%%%%%%%%%%%%GC含量%%%%%%%%%%%
    px= size(PopPos004,1);
index = 1;
temp=[];
for q = 1:px
    if   CTcontent(PopPos004(q,:))~=1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
     PopPos004(temp,:)=[];
end
%%%%%%%%%%%%%%%四C碱基不连续约束
% %==================
index2=1;
temp=[];
pxx = size(PopPos004,1);
for q = 1:pxx
    if   Discontinuous(PopPos004(q,:))==1
        temp(index2) =q;
        index2 = index2 + 1;
    end
end
if length(temp)~=0
     PopPos004(temp,:)=[];
end  
% %%%%% %添加末端约束（左边不能是TC,右边不能是CT(末端是T，次末端是C)）
% x11 = size(PopPos004,1);
% index = 1;
% temp=[];
% % XX3=NewDNA1(:,2*Dim+1:3*Dim);
% for q = 1:x11
%     if  Endlimit(PopPos004(q,:))==1
%         temp(index) =q;
%         index = index + 1;
%     end
% end
% if length(temp)~=0
%    PopPos004(temp,:)=[];
% end   
% %%%%%%%%%%%======================
    %=================扭应力约束=======================================================

x5 = size(PopPos004,1);
index = 1;
temp=[];
% XX3=NewDNA1(:,2*Dim+1:3*Dim);
for q = 1:x5
    if  CTC_TCTcontent(PopPos004(q,:))==1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
   PopPos004(temp,:)=[];
end   
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%            PopPos2=unique(PopPos2,'row','stable');
%               PopFit2=[];
         for i=1:size(PopPos004,1)
             PopFit004(i)=Eint(PopPos004(i,:));
         end

%            NewDNA1=[];
     for i=1:size(PopPos004,1)
        if PopFit004(i)<=standard
%           if SFfitness(i)>standard
            NewDNA004=[PopPos004(i,:);NewDNA004];
        end
     end  
%      %精英挑选
%          [B,min_index] = min(PopFit004);
%            BestF=B;
%            BestPos= PopPos004(min_index,:);%Pick out the best particles and their fitness functions
  %==========================================================   
  
  
%     Add torsional stress constraint
NewDNA1=NewDNA004;
if size(NewDNA1,1)==0
    NewDNA1=[];
else
x5 = size(NewDNA1,1);
index = 1;
temp=[];
% XX3=NewDNA1(:,2*Dim+1:3*Dim);
for q = 1:x5
    if  CTC_TCTcontent(NewDNA1(q,:))==1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
   NewDNA1(temp,:)=[];
end   
end
% ========================================
  
% ===========================================
      %I-motif free constraint 
 if size(NewDNA1,1)==0
    NewDNA1=[];
else
x = size(NewDNA1,1);
index = 1;
temp=[];
% Xxn=NewDNA1(:,1:Dim);
for q = 1:x
    if  Imotif(NewDNA1(q,:))==0
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
   NewDNA1(temp,:)=[];
end   
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %饥饿over
  %%%%%%%%%%把更新后的位置传到下一代的位置========
        PopPos=newPopPos;
  %%%%%%%%=================================
%        NewDNA=[NewDNA003;NewDNA004];
       NewDNA=NewDNA003;
       DNASet=[DNASet;NewDNA];
     DNASet = unique( DNASet,'rows','stable');
     for i=1:size(DNASet,1)
     DNAFitness(i) = Eint(DNASet(i,:));
     end
%      for w=1:1
%          if NewDNA==[]
% %              break
%     [DNAFitnessSorted,IndexSorted] = sort(DNAFitness);
%      DNAbestFitness = DNAFitnessSorted(1);
     
      [Firstfitness,xxx]=min( DNAFitness);
      [Endfieness,xx]=max( DNAFitness);

     Destination_fitness= Firstfitness;
      Worstest_fitness=Endfieness;
      
     
      DNAelite=DNASet(xxx,:);
      DNArubbish=DNASet(xx,:);

    
      fprintf('In the %d generation of evolution, the size of DNASet is %d\n',It,size( DNASet,1));
    toc
end
%     HisBestF(It)=BestF;
% Arry_Fitness(It)=BestF;
% toc
end
% Arry_Fitness(irun)=BestF;
% 
% display(['Run num : ', num2str(irun)]);
% % display(['The best solution obtained by HARO is : ', num2str(Best_pos),10]);
% display(['The best optimal value of the objective funciton found by HARO is : ', num2str(BestF),10]);
% display(sprintf('==============================='));

%  end
 
 
%   Ave=mean(Arry_Fitness);
% Std=std(Arry_Fitness);
% 
% % display(['Run mn : ', num2str(mn)]);
% % display(['The best solution obtained by HARO is : ', num2str(Best_pos),10]);
% % display(['The best optimal value of the objective funciton found by HARO is : ', num2str(BestF),10]);
%     display(['The Ave by HARO is : ', num2str(Ave,10)]);
%     display(['The Std by HARO is : ', num2str(Std,10)]);
% end
 
% end
 

