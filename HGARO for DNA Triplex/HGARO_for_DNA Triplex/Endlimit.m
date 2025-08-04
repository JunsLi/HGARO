%ĞòÁĞÄ©¶ËÏŞÖÆ
function  flag = Endlimit(candidate,Dim)
Dim=10; 
py = size(candidate,2);
count = 0;

%     if candidate(1) ==1 && candidate(Dim) ==1
    if candidate(1,1) ==3 && candidate(1,2) ==2  
%     if candidate(1,1) ==3 && candidate(1,2) ==2  &&  candidate(1,Dim) ==3 && candidate(1,Dim-1) ==2   
        
%     if  candidate(1) ==1 
           count = count + 1;
    end
    if candidate(1,Dim) ==3 && candidate(1,Dim-1) ==2  
%     if  candidate(1) ==1 
           count = count + 1;
    end
    

if count >0
    flag = 1;
else
    flag =0;
end
  %A=0
  %G=1
  %C=2
  %T=3
