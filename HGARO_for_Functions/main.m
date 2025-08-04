%--------------------------------------------------------------------------
%%% Artificial Rabbits Optimization (ARO) for 23 functions %%%
% ARO code v1.0.
% Developed in MATLAB R2011b
% --------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BestX:The best solution                  %
% BestF:The best fitness                   %
% HisBestF:History of the best fitness     %
% FunIndex：Index of functions             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
dimSize = 30;
MaxIteration=1000;
PopSize=500;
% FunIndex=20;
FunIndex='F5';
    [Low,Up,Dim,fobj]=FunRange(FunIndex);
    [BestX,BestF,HisBestF,Ave,Std]=HARO(FunIndex,MaxIteration,PopSize,fobj);
    % display(['FunIndex=', num2str(FunIndex)]);
%     display(['The best fitness of F',num2str(FunIndex),' is: ', num2str(BestF)]);
    display(['The Ave by HARO is : ', num2str(Ave,10)]);
    display(['The Std by HARO is : ', num2str(Std,10)]);
    
    %display(['The best solution is: ', num2str(BestX)]);

    if BestF>0
        semilogy(HisBestF,'r','LineWidth',2);
    else
        plot(HisBestF,'r','LineWidth',2);
    end

    xlabel('Iterations');
    ylabel('Fitness');
    title(['F',num2str(FunIndex)]);


