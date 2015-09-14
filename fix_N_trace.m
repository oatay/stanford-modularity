%function that fixes outlyers in traces
function [F1N]=fix_N_trace(F1N)

LEN=length(F1N);

to_fix=zeros(1,LEN);
for i=2:LEN-1
    B=F1N(1,i-1)-F1N(1,i);
    A=F1N(1,i)  -F1N(1,i+1);
    if sign(B.*A)==-1 && (abs(A)>4000) && abs(B)>4000
    to_fix(1,i)=1;
    end   
end

X=find(to_fix==1);
for i=1:length(X)
    j=X(i);
    F1N(1,j)= 0.5.*(F1N(1,j-1)+ F1N(1,j+1));
end