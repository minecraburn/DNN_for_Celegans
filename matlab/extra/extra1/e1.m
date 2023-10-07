sttr2_merge=zeros(6,6);
sttr2_merge(1:4,1:4)=sttr2(1:4,1:4);
sttr2_merge(6,6)=sttr2(7,7);
sttr2_merge(1:4,6)=sttr2(1:4,7);
sttr2_merge(6,1:4)=sttr2(7,1:4);
sttr2_merge(5,1:4)=sttr2(5,1:4)+sttr2(6,1:4);
sttr2_merge(5,6)=sttr2(5,7)+sttr2(6,7);
sttr2_merge(1:4,5)=sttr2(1:4,5)+sttr2(1:4,6);
sttr2_merge(6,5)=sttr2(7,5)+sttr2(7,6);
sttr2_merge=sttr2_merge/sum(sum(sttr2_merge));
sttr2_merge=sttr2_merge+1e-10;
sttr2_merge=sttr2_merge/sum(sum(sttr2_merge));
stpdata_merge=zeros(6,1);
stpdata_merge(1:4)=stpdata(1:4);
stpdata_merge(5)=stpdata(5)+stpdata(6);
stpdata_merge(6)=stpdata(7);
stpdata_merge=stpdata_merge/sum(sum(stpdata_merge));
kl_array=zeros(30,1);
klp_array=zeros(30,1);
% x_array=[1e-6,3e-6,6e-6,1e-5,3e-5,6e-5,1e-4,3e-4,6e-4,1e-3,3e-3,6e-3,1e-2,3e-2,6e-2,1e-1,3e-1,6e-1,1e-0];
x_array=zeros(30,1);
for i=1:1:30
    x_array(i)=1e-6*(1.5^(i-1));
end
for i=1:1:30
    sttr3=squeeze(tranp(i,:,:));
    stp3=squeeze(tranlinep(i,:));
    sttr3_merge=zeros(6,6);
    sttr3_merge(1:4,1:4)=sttr3(1:4,1:4);
    sttr3_merge(6,6)=sttr3(7,7);
    sttr3_merge(1:4,6)=sttr3(1:4,7);
    sttr3_merge(6,1:4)=sttr3(7,1:4);
    sttr3_merge(5,1:4)=sttr3(5,1:4)+sttr3(6,1:4);
    sttr3_merge(5,6)=sttr3(5,7)+sttr3(6,7);
    sttr3_merge(1:4,5)=sttr3(1:4,5)+sttr3(1:4,6);
    sttr3_merge(6,5)=sttr3(7,5)+sttr3(7,6);
    stp3_merge(1:4)=stp3(1:4);
    stp3_merge(5)=stp3(5)+stp3(6);
    stp3_merge(6)=stp3(7);
    stp3_merge=stp3_merge/sum(sum(stp3_merge));
    sttr3_merge=sttr3_merge/sum(sum(sttr3_merge));
    sttr3_merge=sttr3_merge+1e-10;
    sttr3_merge=sttr3_merge/sum(sum(sttr3_merge));
    kl=0;
    klp=0;
    for j=1:6
        for k=1:6
            kl=kl-sttr2_merge(j,k)*log(sttr3_merge(j,k));
        end
        klp=klp-stpdata_merge(j)*log(stp3_merge(j));
    end
    kl_array(i)=kl;
    klp_array(i)=klp;
end
semilogx(x_array,kl_array);
hold on
semilogx(x_array,klp_array);
plot([6e-4,6e-4],[1.5,5.5]);
legend("KL of state transition probability","KL of state probability")