function ST = Eint (candidate,Dim)
Dim=12;
E=0;

for i=1:Dim-1
    if candidate(1,i) == 3 && candidate(1,i+1) == 3
        E=E-34;
    elseif(candidate(1,i) == 3 && candidate(1,i+1) == 2)
         E=E-44.4;
    elseif candidate(1,i) == 2 && candidate(1,i+1) == 3
         E=E-44.4;
        else
        E=E-48.8;
    end
end
ST = E;
% if stack <= standard
%     ST = 1;
% else
%     ST =0;
% end

  %A=0
  %G=1
  %C=2
  %T=3