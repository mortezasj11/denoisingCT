% horizontal vector is shifted towards right
function [M]=fshift(M,shiftnum);

n=size(M);
row=n(1);
col=n(2);

if (shiftnum<1)|(shiftnum>=col)
    display('shiftnum is wrong!');
    return;
end

Mtemp=zeros(n);
for i=1:col-shiftnum
    Mtemp(i+shiftnum)=M(i);
end
for i=(col-shiftnum+1):col
    Mtemp(i+shiftnum-col)=M(i);
end

M=Mtemp;
end
