
sttr=zeros(9,9);
sttr2=zeros(7,7);
stpdata=zeros(7,1);
temrec=zeros(10000,4);
temid=zeros(4,1);
temlsid=0;
temco=0;
temok=0;
for k=1:1:n0
    state1=WT_NoStim(k).States;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp2=size(state1);
    for i=1:1:temp2(2)
        if state1(i)==6
            state1(i)=5;
        end
    end
    for i=1:1:temp2(2)-1
        if state1(i)~=8
            stpdata(state1(i))=stpdata(state1(i))+1;
        end
        if state1(i)~= state1(i+1)
            sttr(state1(i),state1(i+1))=sttr(state1(i),state1(i+1))+1;
            if state1(i)== 8 || state1(i+1)==8
                break
            end
            sttr2(state1(i),state1(i+1))=sttr2(state1(i),state1(i+1))+1;
            if state1(i)~=2
                temco=0;
                temok=0;
            end
            if temok==1
                temco=temco+1;
            end
            if (state1(i)==3 || state1(i)==4) && state1(i+1)==2
                temok=1;temlsid=state1(i)-2;
                sttr2(state1(i),state1(i+1))=sttr2(state1(i),state1(i+1))-1;
                sttr2(state1(i),5)=sttr2(state1(i),5)+1;
            end
            if state1(i+1)==5  && state1(i)==2 && temok==1
                temid(temlsid)=temid(temlsid)+1;
                temrec(temid(temlsid),temlsid)=temco;
                sttr2(state1(i),state1(i+1))=sttr2(state1(i),state1(i+1))-1;
%                 sttr(temlsid+2,2)=sttr(temlsid+2,2)-1;
%                 sttr(temlsid+2,9)=sttr(temlsid+2,9)+1;
%                 sttr(2,state1(i+1))=sttr(2,state1(i+1))-1;
%                 sttr(9,state1(i+1))=sttr(9,state1(i+1))+1;
                temco=0;
                temok=0;
            end
%             if state1(i)==8
%                 "to"
%                 i
%                 state1(i+1)
%             end
%              if state1(i+1)==8
%                 "from"
%                 k
%                 i
%                 state1(i)
%             end           
        end
    end
end

stnum=zeros(1,8);
strec=zeros(3000,8);
for i=1:1:2000
    stnum(state1(i))=stnum(state1(i))+1;
    strec(stnum(state1(i)),state1(i))=i;
end

condis=[0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08];
% sco1=score(1:3137,:);

%聚类预备工作：查看初次聚类停止的合适距离
% avdis=0;
% lenum1=0;
% see=8;
% for i=1:1:stnum(see)
%     if norm(sco1(strec(2,see):strec(2,see),:)-sco1(strec(i,see):strec(i,see),:))<0.15
%         lenum1=lenum1+1;
%         avdis=avdis+norm(sco1(strec(2,see):strec(2,see),:)-sco1(strec(i,see):strec(i,see),:));
%     end
% end
% avdis=avdis/lenum1;

%初次聚类 理论最大复杂度为（n^3m）约10^10计算量，实际上经过优化后得到一定量级的下降
%目标：数组leishu,为k聚类奠定基础。


% leishu=zeros(1,8);
% starr=zeros(1500,8);
% stconu=zeros(1500,8);
% stcogo=zeros(14,1500,8);
% disco2=zeros(1500,1500,8);
% for i=1:1:8
%     labigdis=0;
%     leishu(i)=stnum(i);
%     for j=1:1:stnum(i)
%         starr(j,i)=j;
%         stcogo(:,j,i)=sco1(strec(j,i),:);
%         stconu(j,i)=1;
%         for k=1:1:stnum(i)
%             disco2(j,k,i)=norm(sco1(strec(j,i),:)-sco1(strec(k,i),:));
%         end
%     end
%     while labigdis<condis(i)
%         leishu(i)=leishu(i)-1;
%         temp3=10;tempj=0;tempk=0;
%         for j=1:1:stnum(i)
%             while(starr(starr(j,i),i)~=starr(j,i))
%                 starr(j,i)=starr(starr(j,i),i);
%             end
%         end
%         for j=1:1:stnum(i)
%             if starr(j,i)~=j
%                 continue;
%             end
%             for k=j+1:1:stnum(i)
%                 if starr(k,i)~=k
%                     continue;
%                 end
%                 if disco2(j,k)<temp3
%                     temp3=disco2(j,k);
%                     tempj=j;
%                     tempk=k;    
%                 end
%             end
%         end
%         stcogo(:,tempj,i)=stcogo(:,tempj,i)*(stconu(tempj,i)/(stconu(tempj,i)+stconu(tempk,i)))+stcogo(:,tempk,i)*(stconu(tempk,i)/(stconu(tempj,i)+stconu(tempk,i)));
%         stconu(tempj,i)=stconu(tempj,i)+stconu(tempk,i);
%         starr(tempk,i)=tempj;
%         for j=1:1:stnum(i)
%             if starr(j,i)~=j
%                 continue;
%             end
%             test=stcogo(:,tempj,i)-stcogo(:,j,i);
%             disco2(j,tempj)=norm(stcogo(:,tempj,i)-stcogo(:,j,i));
%             disco2(tempj,j)=norm(stcogo(:,tempj,i)-stcogo(:,j,i));
%         end
%         labigdis=temp3;
%     end
% end

% stcoar=zeros(150,8);
%  for i=1:1:8
%      for j=1:1:stnum(i)
%          while(starr(starr(j,i),i)~=starr(j,i))
%                  starr(j,i)=starr(starr(j,i),i);
%          end
%      end
%      k=0;
%      for j=1:1:stnum(i)
%          if starr(j,i)==j
%              k=k+1;
%              stcoar(k,i)=j;
%          end
%      end
%  end



leishu=[55,61,34,124,47,61,80,3];



%开始k聚类

% starr2=zeros(1500,8);
% stconu2=zeros(150,8);
% stcogo2=zeros(14,150,8);
% 
% for i=1:1:8
%     test=0;
%     eder1=0;
%     for j=1:1:leishu(i)
%         stcogo2(:,j,i)=sco1(strec(j*round(stnum(i)/leishu(i)-1),i),:);
%         stconu2(j,i)=0;
%     end
%     while eder1==0
%         test=test+1;
%         eder1=1;
%         for j=1:1:stnum(i)
%             temp4=10;templei=0;
%             for k=1:1:leishu(i)
%                 if norm(sco1(strec(j,i),:)-stcogo2(:,k,i))<temp4
%                     temp4=norm(sco1(strec(j,i),:)-stcogo2(:,k,i));
%                     templei=k;
%                 end
%             end
%             if templei~=starr2(j,i)
%                 eder1=0;
%                 starr2(j,i)=templei;
%             end
%         end
%         for j=1:1:leishu(i)
%             stcogo2(:,j,i)=zeros(1,14);
%             stconu2(j,i)=0;
%         end
%         for j=1:1:stnum(i)
%             stcogo2(:,starr2(j,i),i)=stcogo2(:,starr2(j,i),i)+sco1(strec(j,i),:)';
%             stconu2(starr2(j,i),i)=stconu2(starr2(j,i),i)+1;
%         end
%         for j=1:1:leishu(i)
%             if stconu2(j,i)~=0
%                 stcogo2(:,j,i)=stcogo2(:,j,i)/stconu2(j,i);
%             end
%         end
%         if test>100
%             break;
%         end
%     end
% end
