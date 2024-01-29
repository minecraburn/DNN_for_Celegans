linep_normal=zeros(15,6);
linep_normal(:,1:4)=linep(:,1:4);
linep_normal(:,5)=linep(:,5)+linep(:,6);
linep_normal(:,6)=linep(:,7);
for i=1:1:15
    linep_normal(i,:)=linep_normal(i,:)/sum(sum(linep_normal(i,:)));
end
for i=2:1:15
    klp(i)=0;
    klp2(i)=0;
    for j=1:1:6
     klp(i)=klp(i)+linep_normal(1,j)*log(linep_normal(1,j)/linep_normal(i,j));
    end
    klp2(i)=klp2(i)+(linep_normal(1,1)+linep_normal(1,2))*log((linep_normal(1,1)+linep_normal(1,2))/(linep_normal(i,1)+linep_normal(i,2)));
    for j=3:1:6
     klp2(i)=klp2(i)+linep_normal(1,j)*log(linep_normal(1,j)/linep_normal(i,j));
    end
end
henzhou={'Forward','Slow','Dorsal turn','Ventral turn','Reverse','Sustain Reverse'};
zongzhou={'None',"AIBL","AIBR","ALA","AVAL","AVAR","AVBL","AVBR","AVER","RIBL","RID","RIML","RIMR","RIVL","RIVR"};
m=16;n=7;
column_name=strcat(henzhou);
row_name=strcat(zongzhou);

set(figure(10000),'position',[20 20 600 400]);
format shortG
uitable(gcf,'Data',linep_normal,'Position',[100 100 441 280],'Columnname', column_name,'Rowname', row_name);
format short
%% 
connect=[119+119,126+126,11+11,493+493,478+478,198+185,203+203,167+167,92+92,39+42,86+90,115+119,22+26,31+36];

bar([klp(2:15)/sum(klp(2:15));connect/sum(connect)]');
name=["AIBL","AIBR","ALA","AVAL","AVAR","AVBL","AVBR","AVER","RIBL","RID","RIML","RIMR","RIVL","RIVR"];
set(gca,'XTick',1:14,'xticklabel',name);
set(gca,'Yticklabel',[]);
legend("KL divergence", "Degree")


%%

xlabel("Log KL divergence",'FontSize',13);
ylabel("Log degree",'FontSize',13);
lgklp=log(klp);
lgcon=log(connect);
xva=0;
xva(1,1:14)=lgklp;
xva(2,1:14)=1;
a=lgcon/xva;
b=a(2);
a=a(1);
hold on
scatter(log(klp),log(connect));
plot([-2.8,0.8],[-2.8*a+b,0.8*a+b]);
xlim([-3,1]);
ylim([2,8]);
set(gca,'FontSize',13);
% legend("r=0.6101")