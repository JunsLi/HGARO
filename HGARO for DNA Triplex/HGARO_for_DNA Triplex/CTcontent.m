
function  flag = CTcontent(candidate,Dim)
[px,py] = size(candidate);
Dim=12;
count = 0;
for i = 1:py
    if candidate(i) == 2 
           count = count + 1;
    end
end
% if  count>=4 && count<= 8
%  if  count==(py/2)
% if  count==6
if count/Dim >=0.4 &&   count/Dim <=0.6

% if  count==10
    flag = 1;
else
    flag =0;

end
  %A=0
  %G=1
  %C=2
  %T=3
