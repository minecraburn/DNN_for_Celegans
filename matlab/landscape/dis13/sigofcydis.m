%i=1000 1002
D=0.000413509842474;
sig=D*eye(14);
winl=200;
cy1=cycle1;
cy2=cycle2;

% for i=winl+1:1:4000
%     dis1=pdist2(cy1(i,:),cy1(winl,:));
%     dis0=pdist2(cy1(i-1,:),cy1(winl,:));
%     dis2=pdist2(cy1(i+1,:),cy1(winl,:));
%     if dis1<dis0 && dis1<dis2 && dis1<0.1
%         i
%         dis1
%     end
% end

% winl2=1;
% for i=winl2+1:1:8000
%     dis1=pdist2(cy2(i,:),cy2(winl2,:));
%     dis0=pdist2(cy2(i-1,:),cy2(winl2,:));
%     dis2=pdist2(cy2(i+1,:),cy2(winl2,:));
%     if dis1<dis0 && dis1<dis2 && dis1<0.01
%         i
%         dis1
%     end
% end

% rec=0
% for i=winl2+1:1:8000
%     dis1=pdist2(cy2(i,:),cy2(i+100,:));
%     if dis1<1e-2
%         i
%         dis1
%         rec=rec+1;
%         if rec>=5
%            break
%         end
%     end
% end
%%
% plot(cycle1(1:200,1),cycle1(1:200,2));
% hold on
% plot(cycle2(1:200,1),cycle2(1:200,2));

%% 
temp1=squeeze(Fi(1,:));
temp1=temp1/sqrt(temp1*temp1');
temp2=rand(13,14);
temp2=temp2-temp2*temp1'*temp1;
temp2=orth(temp2')';
temp3=[temp1;temp2];
temp4=temp3*temp3';
%%
% setting field
have_cy2=1;
retr_cy1=1;


%coe0.02
% winl=1000;
% winl2=200;
% winr=1834;
% winr2=929;

%coe0.04
% winl=200;
% winl2=1;
% winr=1033;
% winr2=2300;

%coe0.08
% winl=200;
% winl2=1;
% winr=968;
% winr2=400;

%coe0.1
% winl=200;
% winr=1046;
% have_cy2=0;

%coe0.2
winl=1;
winl2=1;
winr=200;
winr2=200;
retr_cy1=0;

%special
% winr2=7495;
h=figure;
movex=0; movey=0;coe3=0.5;coe2=1.8;dicoe=2;dicoe2=dicoe;test_sigu2=10;

% % coe0.02 
% xsidel=-0.8;xsider=0.9;ysidel=-0.6;ysider=0.3;%%ysider=0.3;
% step=0.0025;
% pc=1;
% nea=1;neb=2;test_sig=4e-8;test_sig2=4e-8;
% coe3=0.1;coe2=1;

% % coe0.04 
% xsidel=-0.8;xsider=0.9;ysidel=-0.6;ysider=0.3;%%ysider=0.3;
% step=0.0025;
% pc=1;
% nea=1;neb=2;test_sig=2e-8;test_sig2=5e-9;test_sigu2=10;
% coe3=0.01;coe2=0.6;

% % coe0.08 
% xsidel=-0.8;xsider=0.9;ysidel=-0.6;ysider=0.3;%%ysider=0.3;
% step=0.0025;
% pc=1;
% nea=1;neb=2;test_sig=3e-8;test_sig2=4e-9;test_sigu2=10;
% coe3=0.01;coe2=0.6;

% % coe0.1
% xsidel=-0.8;xsider=0.9;ysidel=-0.6;ysider=0.3;%%ysider=0.3;
% step=0.0025;
% pc=1;
% nea=1;neb=2;test_sig=3e-8;test_sig2=3e-8;
% coe2=0.0005; coe3=0;
 
% % coe0.2 
xsidel=-0.8;xsider=0.9;ysidel=-0.6;ysider=0.3;%%ysider=0.3;
step=0.0025;
pc=1;
nea=1;neb=2;test_sig=4e-8;test_sig2=4e-8;
coe2=5e-15; coe3=0;

% coe0.3 
% xsidel=-0.8;xsider=0.9;ysidel=-0.6;ysider=0.3;%%ysider=0.3;
% step=0.0025;
% pc=1;
% nea=1;neb=2;test_sig=4e-8;test_sig2=4e-8;
% coe2=0; coe3=0;

%fig8a
% xsidel=0;xsider=2;ysidel=0;ysider=1.2;%%ysider=0.3;
% step=0.0025;
% pc=0;
% nea=14;neb=2;test_sig=4e-8;test_sig2=4e-8;
% coe2=0; coe3=0;

cy_n="cycle1";
re_la=1;re_tr=0;have_tr=0;have_legend=0;
re_cy1=1;
% qianzhui="dis/1/";
top_g=50;

D_sigma=0.0055483;
% D_sigma=sqrt(0.0055483);
% D_sigma=sqrt(mean(simsigma(:,1,1))^2+mean(simsigma(:,2,2))^2);

%%
if re_cy1==1
    all_Sig=zeros(12001,14,14);
    all_Sig1=zeros(12001,14,14);
    all_Sig2=zeros(12001,14,14);
    for i=1:1:12000
        temp1=squeeze(Fi(i,:));
        temp1=temp1/sqrt(temp1*temp1');
        ok=0;
        while ok==0
            temp2=rand(13,14);
            temp2=temp2-temp2*temp1'*temp1;
            temp2=orth(temp2')';
            tempn=size(temp2);
            if tempn(1)==13
                ok=1;
            end
        end
        parF=squeeze(parFij(i,:,:));
        parF=temp2*parF*temp2';
        Inmins1=eye(13);
        tempA=kron(Inmins1,parF')+kron(parF',Inmins1);
        [tempL,tempU]=lu(tempA);
        vecb=reshape(-2*D*Inmins1,[],1);
        vecSig=tempU\(tempL\vecb);
        sigless=reshape(vecSig,13,13);
        %     res=parF'*sigless+sigless*parF+2*D*Inmins1;
        sig=temp2'*sigless*temp2;
        all_Sig1(i,:,:)=sig+D*temp1'*temp1;
    end
    all_Sig1(12001,:,:)=all_Sig1(2,:,:);
    for s=1:1:4
        if retr_cy1<0.5
            break
        end
        for i=1:1:12000
            temp1=squeeze(Fi(i,:));
            temp1=temp1/sqrt(temp1*temp1');
            ok=0;
            while ok==0
                temp2=rand(13,14);
                temp2=temp2-temp2*temp1'*temp1;
                temp2=orth(temp2')';
                tempn=size(temp2);
                if tempn(1)==13
                    ok=1;
                end
            end
            parF=squeeze(parFij(i,:,:));
            parF=temp2*parF*temp2';
            Inmins1=eye(13);
            tempA=kron(Inmins1,parF')+kron(parF',Inmins1);
            [tempL,tempU]=lu(tempA);
            tempB=temp2*squeeze(all_Sig1(i+1,:,:)-all_Sig1(i,:,:))*temp2';
            vecb=reshape(-2*D*Inmins1+tempB/2,[],1);
            vecSig=tempU\(tempL\vecb);
            sigless=reshape(vecSig,13,13);
            %     res=parF'*sigless+sigless*parF+2*D*Inmins1;
            sig=temp2'*sigless*temp2;
            all_Sig1(i,:,:)=sig+D*temp1'*temp1;
        end
        all_Sig1(12001,:,:)=all_Sig1(2,:,:);
    end
    all_Sig1=all_Sig1/10;
    if have_cy2>0
        for i=1:1:12000
            temp1=squeeze(Fi2(i,:));
            temp1=temp1/sqrt(temp1*temp1');
            ok=0;
            while ok==0
                temp2=rand(13,14);
                temp2=temp2-temp2*temp1'*temp1;
                temp2=orth(temp2')';
                tempn=size(temp2);
                if tempn(1)==13
                    ok=1;
                end
            end
            parF=squeeze(parFij2(i,:,:));
            parF=temp2*parF*temp2';
            Inmins1=eye(13);
            tempA=kron(Inmins1,parF')+kron(parF',Inmins1);
            [tempL,tempU]=lu(tempA);
            vecb=reshape(-2*D*Inmins1,[],1);
            vecSig=tempU\(tempL\vecb);
            sigless=reshape(vecSig,13,13);
            %     res=parF'*sigless+sigless*parF+2*D*Inmins1;
            sig=temp2'*sigless*temp2;
            all_Sig2(i,:,:)=sig+D*temp1'*temp1;
        end
        all_Sig2=all_Sig2/10;
    end
%     for i=1:1:12000
%         sig=(sig*squeeze(parFij(i,:,:))+squeeze(parFij(i,:,:))'*sig+D)/20;
%         all_Sig(i,:,:)=sig+Fi(i,:)'*Fi(i,:)*D;
%     end
end


all_N=size(simsigma);
all_N=all_N(1);

if pc==1
    read_simpoi=simpoi;
    read_simsigma=simsigma;
    readall_Sig1=all_Sig1;
%     read_score=score;
    readall_Sig2=all_Sig2;
    read_cy1=cy1;
    read_cy2=cy2;
end

if pc==0
    read_simpoi=simpoi*invcoe.*cc3+cc0;
    read_simsigma=simsigma*0;
    for i= 1:1:all_N
        read_simsigma(i,:,:)=diag(cc3)'*invcoe'*squeeze(simsigma(i,:,:))*invcoe*diag(cc3);
    end
    read_score=score*invcoe.*cc3+cc0;
    read_cy1=cy1*invcoe.*cc3+cc0;
    read_cy2=cy2*invcoe.*cc3+cc0;
    readall_Sig1=all_Sig1;
    readall_Sig2=all_Sig2;
    for i= 1:1:winr
        readall_Sig1(i,:,:)=diag(cc3)'*invcoe'*squeeze(all_Sig1(i,:,:))*invcoe*diag(cc3);
        readall_Sig2(i,:,:)=diag(cc3)'*invcoe'*squeeze(all_Sig2(i,:,:))*invcoe*diag(cc3);
    end
end

readall_Sig1=readall_Sig1/dicoe;
readall_Sig2=readall_Sig2/dicoe2;
read_simsigma=read_simsigma/dicoe;


all_N=size(simpoi);
all_N=all_N(1);
xla=[xsidel:step:xsider];
yla=[ysidel:step:ysider];
if re_la==1
    zla1=0*meshgrid(xla,yla);
    zla2=0*meshgrid(xla,yla);
    zla3=0*meshgrid(xla,yla);
end
[xla,yla]=meshgrid(xla,yla);
if re_la==1
    j=0;
    temp_use1=0;
    for k=winl:1:winr
        my_Sigma=squeeze(readall_Sig1(k,[nea,neb],[nea,neb]));
        if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 || det(my_Sigma)> test_sig/dicoe/dicoe%%4e-8/dicoe/dicoe
            continue;
        end
        inv_Sigma=inv(my_Sigma);
        tpx=xla-read_cy1(k,nea);
        tpy=yla-read_cy1(k,neb);
        temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
        temp2=exp(-temp1/2);
        temp_use1=temp_use1+1;
        temp3=temp2/2/pi/sqrt(det(my_Sigma));
        zla1=zla1+temp3;
    end
    zla1=zla1/temp_use1;
    if coe3>0
        j=0;
        temp_use2=0;
        for k=winl2:1:winr2
            my_Sigma=squeeze(readall_Sig2(k,[nea,neb],[nea,neb]));
            if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1/test_sigu2 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>test_sigu2 || det(my_Sigma)/(my_Sigma(2,2)*my_Sigma(2,2))<1/test_sigu2 || det(my_Sigma)/(my_Sigma(2,2)*my_Sigma(2,2))>test_sigu2 || det(my_Sigma)> test_sig2/dicoe/dicoe%%4e-8/dicoe/dicoe
                continue;
            end
            inv_Sigma=inv(my_Sigma);
            tpx=xla-read_cy2(k,nea);
            tpy=yla-read_cy2(k,neb);
            temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
            temp2=exp(-temp1/2);
            temp_use2=temp_use2+1;
            temp3=temp2/2/pi/sqrt(det(my_Sigma));
            zla3=zla3+temp3;
        end
        zla3=zla3/temp_use2;
    end
    j=0;
    temp_use3=0;
    for k=1:1:all_N
%         if k~=4
%             continue
%         end
        my_Sigma=squeeze(read_simsigma(k,[nea,neb],[nea,neb]));
        if det(my_Sigma)<0 ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 
            continue;
        end
%         if det(my_Sigma)>2e-7
%             my_Sigma=my_Sigma/2;
%         end
        inv_Sigma=inv(my_Sigma);
        tpx=xla-read_simpoi(k,nea);
        tpy=yla-read_simpoi(k,neb);
        temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
        temp2=exp(-temp1/2);
        temp_use3=temp_use3+1;
        temp3=temp2/2/pi/sqrt(det(my_Sigma));
        zla2=zla2+temp3;
    end
    zla2=zla2/temp_use3;
end

if coe2==-1
    zla=-log(zla1);
elseif coe2==0
    zla=-log(zla2);
elseif coe3==0
    zla=(-log((coe2*zla1+zla2)/(1+coe2)));
else
%      zla=(-log(zla1)-coe*log(zla2))/(1+coe);
      zla=(-log((coe2*zla1+zla2+coe3*zla3)/(1+coe2+coe3)));
%       zla=(-log((coe2*zla1+zla2)/(1+coe2)));
%       zla=-log(zla3);
end
surf(xla,yla,min(top_g,zla),'DisplayName','');shading interp;
hold on






if have_tr==1
    xtr1=read_cy1(winl:winr,nea);
    ytr1=read_cy1(winl:winr,neb);
    ztr11=0*xtr1;
    ztr12=0*xtr1;
    ztr13=0*xtr1;
    if re_la==1
        j=0;
        temp_use=0;
        for k=winl:1:winr
            my_Sigma=squeeze(readall_Sig1(k,[nea,neb],[nea,neb]));
            if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 || det(my_Sigma)> test_sig/dicoe/dicoe%%4e-8/dicoe/dicoe
                continue;
            end
            inv_Sigma=inv(my_Sigma);
            tpx=xtr1-read_cy1(k,nea);
            tpy=ytr1-read_cy1(k,neb);
            temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
            temp2=exp(-temp1/2);
            temp_use=temp_use+1;
            temp3=temp2/2/pi/sqrt(det(my_Sigma));
            ztr11=ztr11+temp3;
        end
        ztr11=ztr11/temp_use;
        j=0;
        temp_use=0;
        for k=winl2:1:winr2
            my_Sigma=squeeze(readall_Sig2(k,[nea,neb],[nea,neb]));
             if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1/test_sigu2 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>test_sigu2 || det(my_Sigma)/(my_Sigma(2,2)*my_Sigma(2,2))<1/test_sigu2 || det(my_Sigma)/(my_Sigma(2,2)*my_Sigma(2,2))>test_sigu2 || det(my_Sigma)> test_sig2/dicoe/dicoe%%4e-8/dicoe/dicoe
                continue;
            end
            inv_Sigma=inv(my_Sigma);
            tpx=xtr1-read_cy2(k,nea);
            tpy=ytr1-read_cy2(k,neb);
            temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
            temp2=exp(-temp1/2);
            temp_use=temp_use+1;
            temp3=temp2/2/pi/sqrt(det(my_Sigma));
            ztr13=ztr13+temp3;
        end
        ztr13=ztr13/temp_use;
        j=0;
        temp_use=0;
        for k=1:1:all_N
            if k~=4
                continue
            end
            my_Sigma=squeeze(read_simsigma(k,[nea,neb],[nea,neb]));
            if det(my_Sigma)<0 ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1
                continue;
            end
            inv_Sigma=inv(my_Sigma);
            tpx=xtr1-read_simpoi(k,nea);
            tpy=ytr1-read_simpoi(k,neb);
            temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
            temp2=exp(-temp1/2);
            temp_use=temp_use+1;
            temp3=temp2/2/pi/sqrt(det(my_Sigma));
            ztr12=ztr12+temp3;
        end
        ztr12=ztr12/temp_use;
    end
    
    if coe2==-1
        ztr1=-log(ztr11);
    elseif coe2==0
        ztr1=-log(ztr12);
    else
        ztr1=(-log((coe2*ztr11+ztr12+coe3*ztr13)/(1+coe2+coe3)))+1;
%         ztr1=(-log((coe2*ztr11+ztr12+ztr13)/(1+coe2)))+1;
    end
    % % pic3f
%     ztr1(10)=-1;    
%     ztr1(11)=-1;    
%     ztr1(12)=-0.2;
%     ztr1(9)=-1;
%     ztr1(8)=-0.2;
%     ztr1(7)=-0.2;
%     ztr1(13)=-0.1;    
    
    line1=plot3(xtr1,ytr1,ztr1,"black");
    
    
    
    xtr2=read_cy2(winl2:winr2,nea);
    ytr2=read_cy2(winl2:winr2,neb);
    ztr21=0*xtr2;
    ztr22=0*xtr2;
    ztr23=0*xtr2;
    if re_la==1
        j=0;
        temp_use=0;
        for k=winl:1:winr
            my_Sigma=squeeze(readall_Sig1(k,[nea,neb],[nea,neb]));
            if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 || det(my_Sigma)> test_sig/dicoe/dicoe%%4e-8/dicoe/dicoe
                continue;
            end
            inv_Sigma=inv(my_Sigma);
            tpx=xtr2-read_cy1(k,nea);
            tpy=ytr2-read_cy1(k,neb);
            temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
            temp2=exp(-temp1/2);
            temp_use=temp_use+1;
            temp3=temp2/2/pi/sqrt(det(my_Sigma));
            ztr21=ztr21+temp3;
        end
        ztr21=ztr21/temp_use;
        j=0;
        temp_use=0;
        for k=winl2:1:winr2
            my_Sigma=squeeze(readall_Sig2(k,[nea,neb],[nea,neb]));
             if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1/test_sigu2 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>test_sigu2 || det(my_Sigma)/(my_Sigma(2,2)*my_Sigma(2,2))<1/test_sigu2 || det(my_Sigma)/(my_Sigma(2,2)*my_Sigma(2,2))>test_sigu2 || det(my_Sigma)> test_sig2/dicoe/dicoe%%4e-8/dicoe/dicoe
                continue;
            end
            inv_Sigma=inv(my_Sigma);
            tpx=xtr2-read_cy2(k,nea);
            tpy=ytr2-read_cy2(k,neb);
            temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
            temp2=exp(-temp1/2);
            temp_use=temp_use+1;
            temp3=temp2/2/pi/sqrt(det(my_Sigma));
            ztr23=ztr23+temp3;
        end
        ztr23=ztr23/temp_use;
        j=0;
        temp_use=0;
        for k=1:1:all_N
            if k~=4
                continue
            end
            my_Sigma=squeeze(read_simsigma(k,[nea,neb],[nea,neb]));
            if det(my_Sigma)<0 ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1
                continue;
            end
            inv_Sigma=inv(my_Sigma);
            tpx=xtr2-read_simpoi(k,nea);
            tpy=ytr2-read_simpoi(k,neb);
            temp1=tpx.^2*inv_Sigma(1,1)+tpy.^2*inv_Sigma(2,2)+2*tpx.*tpy*inv_Sigma(1,2);
            temp2=exp(-temp1/2);
            temp_use=temp_use+1;
            temp3=temp2/2/pi/sqrt(det(my_Sigma));
            ztr22=ztr22+temp3;
        end
        ztr22=ztr22/temp_use;
    end
    
    if coe2==-1
        ztr2=-log(ztr21);
    elseif coe2==0
        ztr2=-log(ztr22);
    else
        ztr2=(-log((coe2*ztr21+ztr22+coe3*ztr23)/(1+coe2+coe3)))+1;
    end

    % % pic3f
%     ztr2(2178)=10;
%     ztr2(2177)=11;
%     ztr2(2176)=12;
%     ztr2(2175)=13;
%     ztr2(2172)=13; 
%     ztr2(2173)=15; 
%     ztr2(2174)=14; 
    line2=plot3(xtr2,ytr2,ztr2,"red");
    
end





name_of_ne=["AIBL","AIBR",	"ALA",	"AVAL",	"AVAR",	"AVBL",	"AVBR",	"AVER",	"RIBL",	"RID",	"RIML",	"RIMR",	"RIVL",	"RIVR",	"RMED",	"RMEL",	"RMER", 	"SMDDR",	"SMDVR",	"VB01",	"VB02"];
if pc==1
xlabel("PC1",'FontSize',12);
ylabel("PC2",'FontSize',12);
end
if pc==0
    xlabel(name_of_ne(nea),'FontSize',12);
    ylabel(name_of_ne(neb),'FontSize',12);
    zlabel("");
end
if have_legend==1
    legend("landscape",cy_n,'Location','Northeast','FontSize',12);
end
if have_tr==1
    legend([line1,line2],"cycle1","cycle2",'FontSize',12);
%     legend("landscaoe","cycle1");
end
cbar=colorbar();
cbar.Position(1)=0.923;
set(gca,'FontSize',12);
set(h,'Units','Inches');
pos=get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
% zlabel("PC3");
% hold on
hold off
%%
%1
axis([xsidel,xsider,ysidel,ysider]);
set(gcf,'position',[3,2,6.6,4.4]);
set(gcf,'PaperSize',[6.6,4.4]);
%%
%2
% if have_tr==0
%     cy_n="noc/landscape";
% end 
% if pc==1 
%     save(qianzhui+"pc/"+cy_n+".mat");
% end
% if pc==0 
%     save(qianzhui+nea+"-"+neb+"/"+cy_n+".mat")
% end



