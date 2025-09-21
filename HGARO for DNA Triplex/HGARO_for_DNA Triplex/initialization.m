
% This function initialize the first population of search agents
%
function PopPos=initialization(nPop,Dim,Up,Low)
       
Boundary_no= size(Up,2); % numnber of boundaries


% t=[2 3 ];
% t(randint(1,1,[1 length(t)]))


if Boundary_no>1
        for i=1:Dim   
%              T=[2,3];
%              t=T(randi(1,[1 length(T)]));
%        DNASet(:,i)=ones(nPop,1).*T(randi(1,[1 length(T)])); 
       
       DNASet(:,i)=floor(unifrnd(2,4,nPop,1));
%         DNASet(:,i)=round(rand(nPop,1).*((Up(i)-Low(i))./2))+Low(i); 
         end
end
%unique(A,'row','stable')
% DNASet=unique(DNASet,'row','stable');
% XSF = unique(DNASet,'rows','stable');
%
% for i=1:size(DNASet,1)
%         % Check if solutions go outside the search space and bring them back
%         Flag4ub=DNASet(i,:)>ub;
%         Flag4lb=DNASet(i,:)<lb;
%         DNASet(i,:)=(DNASet(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;       
% end
%
%  for p = 1:size(DNASet,1)
%      for q = 1:size(DNASet,2)-1
%          if DNASet(p,q) == DNASet(p,q+1)
%              del_idx(p,:) = p;
%          end
%      end
%  end
%  del_idx(del_idx(:,end)==0) = [];
%  DNASet(del_idx,:) = [];
% %%%%%%%%%%%%%%%%%%%%
% px = size(DNASet,1);
% index = 1;
% temp=[];
% for q = 1:px
%     if   CTcontent(DNASet(q,:))~=1
%         temp(index) =q;
%         index = index + 1;
%     end
% end
% if length(temp)~=0
%     DNASet(temp,:)=[];
% end
% %=====================
% %
% % %==================
% index2=1;
% temp=[];
% pxx = size(DNASet,1);
% for q = 1:pxx
%     if   Discontinuous(DNASet(q,:))==1
%         temp(index2) =q;
%         index2 = index2 + 1;
%     end
% end
% if length(temp)~=0
%     DNASet(temp,:)=[];
% end
% %======================
% pxxx = size(DNASet,1);
% index3=1;
% temp=[];
% for q = 1:pxxx
%     if   Stack(DNASet(q,:))~=1
%         temp(index3) =q;
%         index3 = index3 + 1;
%     end
% end
% if length(temp)~=0
%     DNASet(temp,:)=[];
% end

% pxxxx = size(DNASet,1);
% Set2=[];
% for i=1:pxxxx
%     for j=1:Dim
%         if DNASet(i,j)==3
%         Set2(i,j)=3;
%         else
%          Set2(i,j)=2;
%         end
%     end
% end
%
% Set2(:,:)=Set2(:,end:-1:1);
% % %loop settings
% % TSet=zeros(pxxxx,5);
% % TSet(TSet==0)=3;
% 
% %
% % DNASet=[DNASet,TSet];
% DNASet0=[DNASet,Set2];

%-------------------------------------------
% cc=size(DNASet0,1);
%  DNASet1=zeros(cc,Dim);
%  DNASet2=[];
% %  TA=[0,3];
% %  CG=[1,2];
% for j=1:size(DNASet0,1)
%         for i=1:Dim   
% %             for t=1:cc
%                 if DNASet0(j,i)==3 
%                 
%                 DNASet1(j,i)=0;
% %                 DNASet1(j,i)=TA(randi(numel(TA),1));
%             else
% %                 DNASet1(j,i)=CG(randi(numel(CG),1));
%                 DNASet1(j,i)=1;
%                 end
% %             end
%         end
% %       D=  repmat(DNASet0(j,:),cc,1);
%      DD=[DNASet0(j,:),DNASet1(j,:)];
%      DNASet2=[DNASet2;DD];
% %      DNASet2=unique(DNASet2,'row','stable');
% %       DNASet2=DD;
% end
%-----------------------------------


% pxx = size(DNASet2,1);
% index4 = 1;
% temp=[];
% DNASet3=[];
% DNASet3=DNASet2(:,2*dim+1:3*dim);
% for q = 1:pxx
%     if   AContent(DNASet3(q,:))==1
%         temp(index4) =q;
%         index4 = index4 + 1;
%     end
% end
% if length(temp)~=0
%     DNASet2(temp,:)=[];
% end
%=========================================
% temp=[];
% index5 = 1;
% pxxx=size(DNASet2,1);
% DNASet3=DNASet2(:,2*dim+1:3*dim);
% for q = 1:pxxx
%     if   CTC_TCTcontent(DNASet3(q,:))==1
%         temp(index5) =q;
%         index5 = index5 + 1;
%     end
% end
% if length(temp)~=0
%     DNASet2(temp,:)=[];
% end
%=======================================================
%  if size(DNASet2,1)==0
%     DNASet2=[];
% else
%  x7 = size(DNASet2,1);
% index = 1;
% temp=[];
% XX6=DNASet2(:,2*dim+1:3*dim);
% for q = 1:x7
%     if  Discontinuous2(XX6(q,:))==1
%         temp(index) =q;
%         index = index + 1;
%     end
% end
% if length(temp)~=0
%    DNASet2(temp,:)=[];
% end 
%           end

%=================================================
PopPos =DNASet;
% PopPos= unique(DNASet,'rows','stable');

