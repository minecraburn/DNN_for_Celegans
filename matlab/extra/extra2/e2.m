linep_normal=zeros(15,6);
linep_normal(:,1:4)=linep(:,1:4);
linep_normal(:,5)=linep(:,5)+linep(:,6);
linep_normal(:,6)=linep(:,7);
for i=1:1:15
    linep_normal(i,:)=linep_normal(i,:)/sum(sum(linep_normal(i,:)));
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