n0=5;dim=14;
a=strings(n0,300);num=zeros(n0,300);
for i=1:1:n0
    copy1=WT_NoStim(i).IDs;%%%%%%%%%%%%%%%%%%%%%%
    siz=length(copy1);
    k=0;
    for j=1:1:siz
        copy2=copy1{j};
        siz2=length(copy2);
        if siz2>0
            for l=1:1:1
                copy3=copy2(l);
                copy4=copy3{1};
                if ischar(copy4)
                    k=k+1;
                    a{i,k}=copy4;
                    num(i,k)=j;
                end
            end
            
        end
    end
end
c=a(1:1,1:100);

for i=2:1:n0
    c=intersect(c,a(i:i,1:100));
end
c1=c;
d1=zeros(n0,200);
d=zeros(1,200);
for k=1:1:n0
    for i=2:1:dim+1
        siz=length(copy1);
        for j=1:1:300
            if c1(i)==a(k,j)
                d(i-1)=num(k,j);
                break;
            end
        end
    end
d1(k:k,:)=d;
end
state1=WT_NoStim(1).States;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state2=WT_NoStim(2).States;
state3=WT_NoStim(3).States;
state4=WT_NoStim(4).States;
state5=WT_NoStim(5).States;