
 
function [Low,Up,Dim,fobj]=FunRange(FunIndex)
 
    Dim=30;

    switch FunIndex
        
        case 'F1'
            fobj = @F1;
            Low=-100;Up=100;
            
            
        case 'F2'
            fobj = @F2;
            Low=-10;Up=10;
            
            
        case 'F3'
            fobj = @F3;
            Low=-100;Up=100;
            
            
        case 'F4'
            fobj = @F4;
            Low=-100;Up=100;
            
            
        case 'F5'
            fobj = @F5;
            Low=-30;Up=30;
            
            
        case 'F6'
            fobj = @F6;
            Low=-100;Up=100;
            
            
        case 'F7'
            fobj = @F7;
            Low=-1.28;Up=1.28;
            
            
        case 'F8'
            fobj = @F8;
            Low=-500;Up=500;
            
            
        case 'F9'
            fobj = @F9;
            Low=-5.12;Up=5.12;
            
            
        case 'F10'
            fobj = @F10;
            Low=-32;Up=32;
            
            
        case 'F11'
            fobj = @F11;
            Low=-600;Up=600;
            
            
        case 'F12'
            fobj = @F12;
            Low=-50;Up=50;
            
            
        case 'F13'
            fobj = @F13;
            Low=-50;Up=50;
            
            
        case 'F14'
            fobj = @F14;
            Low=-65.536;Up=65.536;Dim=2;
            
            
        case 'F15'
            fobj = @F15;
            Low=-5;Up=5;Dim=4;
            
            
        case 'F16'
            fobj = @F16;
            Low=-5;Up=5;Dim=2;
            
            
        case 'F17'
            fobj = @F17;
            Low=[-5 0];Up=[10 15];Dim=2;
            
            
        case 'F18'
            fobj = @F18;
            Low=-2;Up=2;Dim=2;
            
            
        case 'F19'
            fobj = @F19;
            Low=0;Up=1;Dim=3;
            
            
        case 'F20'
            fobj = @F20;
            Low=0;Up=1;Dim=6;
            
            
        case 'F21'
            fobj = @F21;
            Low=0;Up=10;Dim=4;
            
            
        case 'F22'
            fobj = @F22;
            Low=0;Up=10;Dim=4;
            
            
        otherwise
            fobj = @F23;
            Low=0;Up=10;Dim=4;
    
    end
% end

% function Fit=BenFunctions(X,FunIndex,Dim)
% 
% 
% 
% 
%       switch FunIndex
%      
%           %%%%%%%%%%%%%%%%%%%%%%%%%%unimodal function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           
%           %Sphere
%           case 'F1'
%               Fit=sum(X.^2);
%               
%               %Schwefel 2.22
%           case 'F2'
%               Fit=sum(abs(X))+prod(abs(X));
%               
%               %Schwefel 1.2
%           case 'F3'
%               Fit=0;
%               for i=1:Dim
%                   Fit=Fit+sum(X(1:i))^2;
%               end
%               
%               %Schwefel 2.21
%           case 'F4'
%               Fit=max(abs(X));
%               
%               %Rosenbrock
%           case 'F5'
%               Fit=sum(100*(X(2:Dim)-(X(1:Dim-1).^2)).^2+(X(1:Dim-1)-1).^2);
%               
%               %Step
%           case 'F6'
%               Fit=sum(floor((X+.5)).^2);
%               
%               %Quartic
%           case 'F7'
%               Fit=sum([1:Dim].*(X.^4))+rand;
%               
%               
%               %%%%%%%%%%%%%%%%%%%%%%%%%%multimodal function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               
%               %Schwefel
%           case 'F8'
%               Fit=sum(-X.*sin(sqrt(abs(X))));
%               
%               %Rastrigin
%           case 'F9'
%               Fit=sum(X.^2-10*cos(2*pi.*X))+10*Dim;
%               
%               %Ackley
%           case 'F10'
%               Fit=-20*exp(-.2*sqrt(sum(X.^2)/Dim))-exp(sum(cos(2*pi.*X))/Dim)+20+exp(1);
%               
%               %Griewank
%           case 'F11'
%               Fit=sum(X.^2)/4000-prod(cos(X./sqrt([1:Dim])))+1;
%               
%               %Penalized
%           case 'F12'
%               a=10;k=100;m=4;
%               Dim=length(X);
%               Fit=(pi/Dim)*(10*((sin(pi*(1+(X(1)+1)/4)))^2)+sum((((X(1:Dim-1)+1)./4).^2).*...
%                   (1+10.*((sin(pi.*(1+(X(2:Dim)+1)./4)))).^2))+((X(Dim)+1)/4)^2)+sum(k.*...
%                   ((X-a).^m).*(X>a)+k.*((-X-a).^m).*(X<(-a)));
%               
%               %Penalized2
%           case 'F13'
%               a=10;k=100;m=4;
%               Dim=length(X);
%               Fit=.1*((sin(3*pi*X(1)))^2+sum((X(1:Dim-1)-1).^2.*(1+(sin(3.*pi.*X(2:Dim))).^2))+...
%                   ((X(Dim)-1)^2)*(1+(sin(2*pi*X(Dim)))^2))+sum(k.*...
%                   ((X-a).^m).*(X>a)+k.*((-X-a).^m).*(X<(-a)));
%               
%               %%%%%%%%%%%%%%%%%%%%%%%%%%fixed-dimensionalmultimodalfunction%%%%%%%%%%%%%%
%               
%               %Foxholes
%           case 'F14'
%               a=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
%                   -32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
%               for j=1:25
%                   b(j)=sum((X'-a(:,j)).^6);
%               end
%               Fit=(1/500+sum(1./([1:25]+b))).^(-1);
%               
%               %Kowalik
%           case 'F15'
%               a=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
%               b=[.25 .5 1 2 4 6 8 10 12 14 16];b=1./b;
%               Fit=sum((a-((X(1).*(b.^2+X(2).*b))./(b.^2+X(3).*b+X(4)))).^2);
%               
%               %Six Hump Camel
%           case 'F16'
%               Fit=4*(X(1)^2)-2.1*(X(1)^4)+(X(1)^6)/3+X(1)*X(2)-4*(X(2)^2)+4*(X(2)^4);
%               
%               %Branin
%           case 'F17'
%               Fit=(X(2)-(X(1)^2)*5.1/(4*(pi^2))+5/pi*X(1)-6)^2+10*(1-1/(8*pi))*cos(X(1))+10;
%               
%               %GoldStein-Price
%           case 'F18'
%               Fit=(1+(X(1)+X(2)+1)^2*(19-14*X(1)+3*(X(1)^2)-14*X(2)+6*X(1)*X(2)+3*X(2)^2))*...
%                   (30+(2*X(1)-3*X(2))^2*(18-32*X(1)+12*(X(1)^2)+48*X(2)-36*X(1)*X(2)+27*(X(2)^2)));
%               
%               %Hartman 3
%           case 'F19'
%               a=[3 10 30;.1 10 35;3 10 30;.1 10 35];c=[1 1.2 3 3.2];
%               p=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
%               Fit=0;
%               for i=1:4
%                   Fit=Fit-c(i)*exp(-(sum(a(i,:).*((X-p(i,:)).^2))));
%               end
%               
%               %Hartman 6
%           case 'F20'
%               af=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
%               cf=[1 1.2 3 3.2];
%               pf=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
%                   .2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
%               Fit=0;
%               for i=1:4
%                   Fit=Fit-cf(i)*exp(-(sum(af(i,:).*((X-pf(i,:)).^2))));
%               end
%               
%               %Shekel 5
%           case 'F21'
%               a=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
%               c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
%               Fit=0;
%               for i=1:5
%                   Fit=Fit-1/((X-a(i,:))*(X-a(i,:))'+c(i));
%               end
%               
%               %Shekel 7
%           case 'F22'
%               a=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
%               c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
%               Fit=0;
%               for i=1:7
%                   Fit=Fit-1/((X-a(i,:))*(X-a(i,:))'+c(i));
%               end
%               
%               %Shekel 10
%           otherwise
%               a=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
%               c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
%               Fit=0;
%               for i=1:10
%                   Fit=Fit-1/((X-a(i,:))*(X-a(i,:))'+c(i));
%               end
% 
%       end
