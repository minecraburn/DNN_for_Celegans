del=1;del2=4;my_D=0.000206754921237;
m=3137;
x=zeros(340*200/del/del,1);
y=zeros(340*200/del/del,1);
ua=zeros(340*200/del/del,1);
va=zeros(340*200/del/del,1);
ub=zeros(340*200/del/del,1);
vb=zeros(340*200/del/del,1);
B0=zeros(2,14);
B0(:,3:14)=surface_B(2:3,:);
B0(1,1)=1;B0(2,2)=1;
v1=B0(1,:);
v1=v1/sqrt(v1*v1');
v2=B0(2,:);
v2=v2-v2*v1'*v1;
v2=v2/sqrt(v2*v2');
the_v=[v1;v2];
temp_use=0;
all_N=size(simpoi);
all_N=all_N(1);
D_sigma=0.0055483;
for i = 1:del:340
    for j =1:del:200
%         x_flux=-0.8+0.005*(i-1);
%         y_flux=-0.6+0.005*(j-1);
        x((j-1)/del+1+(i-1)/del*200/del)=x1(i,j,1);
        y((j-1)/del+1+(i-1)/del*200/del)=x1(i,j,2);
        fv=squeeze(fx(i,j,:));
        fv=v1'*(v1*fv)+v2'*(v2*fv);
        ua((j-1)/del+1+(i-1)/del*200/del)=fv(1);
        va((j-1)/del+1+(i-1)/del*200/del)=fv(2);
    end
end

ub=0;
vb=0;
fp1=0;fp2=0;
j=0;sumland1=0;
for k=winl:1:winr
    my_Sigma=squeeze(readall_Sig1(k,[nea,neb],[nea,neb]));
    if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 || det(my_Sigma)> test_sig/dicoe/dicoe%%4e-8/dicoe/dicoe
        continue;
    end
    my_Sigma=the_v*squeeze(readall_Sig1(k,:,:))*the_v';
    if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 || det(my_Sigma)> test_sig/dicoe/dicoe/4
        continue;
    end
    tempsf=[1,read_cy1(k,1:2)]*surface_B-read_cy1(k,3:14);
    tempsf=[0,0,tempsf];
    tempsfno=tempsf/sqrt(tempsf*tempsf');
    temp5=tempsfno*squeeze(readall_Sig1(k,:,:))*tempsfno';
    temp5=tempsf*tempsf'/temp5;
    temp5=max(0,temp5);
    tpx=x-read_cy1(k,1);
    tpy=y-read_cy1(k,2);
    invsig=inv(my_Sigma);
    tempv1=-invsig(1,1)*tpx-invsig(1,2)*tpy;
    tempv2=-invsig(2,1)*tpx-invsig(2,2)*tpy;
    temp3=-tpx.*tempv1-tpy.*tempv2;
    temp4=exp(-temp3/2-temp5/200)/2/pi/sqrt(det(my_Sigma))/temp_use1;
    fp1=fp1+tempv1.*temp4;
    fp2=fp2+tempv2.*temp4;
    sumland1=sumland1+temp4;
end
fp1=fp1./sumland1;fp2=fp2./sumland1;
fp1(isnan(fp1))=0;fp2(isnan(fp2))=0;
ub=ub+fp1/(1+coe2+coe3);
vb=vb+fp2/(1+coe2+coe3);
fp1=0;fp2=0;
j=0;sumland2=0;
for k=winl2:1:winr2
    my_Sigma=squeeze(readall_Sig2(k,[nea,neb],[nea,neb]));
    if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 || det(my_Sigma)> test_sig2/dicoe/dicoe%%4e-8/dicoe/dicoe
        continue;
    end
    my_Sigma=the_v*squeeze(readall_Sig2(k,:,:))*the_v';
    if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 || det(my_Sigma)> test_sig2/dicoe/dicoe/4
        continue;
    end
    tempsf=[1,read_cy2(k,1:2)]*surface_B-read_cy2(k,3:14);
    tempsf=[0,0,tempsf];
    tempsfno=tempsf/sqrt(tempsf*tempsf');
    temp5=tempsfno*squeeze(readall_Sig2(k,:,:))*tempsfno';
    temp5=tempsf*tempsf'/temp5;
    temp5=max(0,temp5);
    tpx=x-read_cy2(k,1);
    tpy=y-read_cy2(k,2);
    invsig=inv(my_Sigma);
    tempv1=-invsig(1,1)*tpx-invsig(1,2)*tpy;
    tempv2=-invsig(2,1)*tpx-invsig(2,2)*tpy;
    temp3=-tpx.*tempv1-tpy.*tempv2;
    temp4=exp(-temp3/2-temp5/2)/2/pi/sqrt(det(my_Sigma))/temp_use2;
    fp1=fp1+tempv1.*temp4;
    fp2=fp2+tempv2.*temp4;
    sumland2=sumland2+temp4;
end
fp1=fp1./sumland2;fp2=fp2./sumland2;
fp1(isnan(fp1))=0;fp2(isnan(fp2))=0;
ub=ub+fp1*coe3/(1+coe2+coe3);
vb=vb+fp2*coe3/(1+coe2+coe3);
fp1=0;fp2=0;
j=0;sumland3=0;
for k=1:1:all_N
    my_Sigma=squeeze(read_simsigma(k,[nea,neb],[nea,neb]));
    if det(my_Sigma)<0 ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1
        continue;
    end
    my_Sigma=the_v*squeeze(read_simsigma(k,:,:))*the_v';
    if det(my_Sigma)<0  ||my_Sigma(1,1)<0 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))<1e-1 || det(my_Sigma)/(my_Sigma(1,1)*my_Sigma(1,1))>1e1 
        continue;
    end
    tempsf=[1,read_simpoi(k,1:2)]*surface_B-read_simpoi(k,3:14);
    tempsf=[0,0,tempsf];
    tempsfno=tempsf/sqrt(tempsf*tempsf');
    temp5=tempsfno*squeeze(read_simsigma(k,:,:))*tempsfno';
    temp5=tempsf*tempsf'/temp5;
    temp5=max(0,temp5);
    tpx=x-read_simpoi(k,1);
    tpy=y-read_simpoi(k,2);
    invsig=inv(my_Sigma);
    tempv1=-invsig(1,1)*tpx-invsig(1,2)*tpy;
    tempv2=-invsig(2,1)*tpx-invsig(2,2)*tpy;
    temp3=-tpx.*tempv1-tpy.*tempv2;
    temp4=exp(-temp3/2-temp5/2)/2/pi/sqrt(det(my_Sigma))/temp_use2;
    fp1=fp1+tempv1.*temp4;
    fp2=fp2+tempv2.*temp4;
    sumland3=sumland3+temp4;
end
fp1=fp1./sumland3;fp2=fp2./sumland3;
fp1(isnan(fp1))=0;fp2(isnan(fp2))=0;
ub=ub+fp1*coe2/(1+coe2+coe3);
vb=vb+fp2*coe2/(1+coe2+coe3);
ub=ub;
vb=vb;
sumland=-log((sumland1+sumland2*coe3+sumland3*coe2)/(1+coe2+coe3));

%% 
% xst=1;xend=340;
% yst=1;yend=200;
xst=40;xend=265;
yst=35;yend=170;
coe0=0.00015
% coe0=0.00015;
% coe0=0;
h=figure;
xsize=size(x);
xsize=xsize(1);
% quiver(x,y,coe0*ub,coe0*vb,15);

reshape_x=reshape(x,200/del,340/del)';
reshape_y=reshape(y,200/del,340/del)';
reshape_ua=reshape(ua,200/del,340/del)';
reshape_ub=reshape(ub,200/del,340/del)';
reshape_va=reshape(va,200/del,340/del)';
reshape_vb=reshape(vb,200/del,340/del)';
% qu1=quiver(reshape_x(1:del2:340/del,1:del2:200/del),reshape_y(1:del2:340/del,1:del2:200/del), ...
% reshape_ua(1:del2:340/del,1:del2:200/del)-coe0*reshape_ub(1:del2:340/del,1:del2:200/del), ...
% reshape_va(1:del2:340/del,1:del2:200/del)-coe0*reshape_vb(1:del2:340/del,1:del2:200/del),3,"red");


hold on
reshape_land=reshape(sumland,200/del,340/del)';
% qu1=quiver3(reshape_x(1:del2:340/del,1:del2:200/del),reshape_y(1:del2:340/del,1:del2:200/del), ...
% 0*reshape_x(1:del2:340/del,1:del2:200/del)+100,...
% reshape_ua(1:del2:340/del,1:del2:200/del)-coe0*reshape_ub(1:del2:340/del,1:del2:200/del), ...
% reshape_va(1:del2:340/del,1:del2:200/del)-coe0*reshape_vb(1:del2:340/del,1:del2:200/del),...
% 0*reshape_x(1:del2:340/del,1:del2:200/del),...
% 2,"red");
[itop,jtop]=size(reshape_x);
for i=1:1:itop
    for j=1:1:jtop
        if reshape_land(i,j)>80
            reshape_ua(i,j)=0;
            reshape_va(i,j)=0;
            reshape_ub(i,j)=0;
            reshape_vb(i,j)=0;
        end
    end
end

reshape_uc=reshape_ua-coe0*reshape_ub;
reshape_vc=reshape_va-coe0*reshape_vb;
reshape_ud=reshape_uc./sqrt(reshape_uc.^2+reshape_vc.^2);
reshape_vd=reshape_vc./sqrt(reshape_uc.^2+reshape_vc.^2);

% qu1=quiver3(reshape_x(xst:del2:xend,yst:del2:yend),reshape_y(xst:del2:xend,yst:del2:yend), ...
% 0*reshape_xxst:del2:xend,yst:del2:yend)+100,...
% reshape_ua(xst:del2:xend,yst:del2:yend), ...
% reshape_va(xst:del2:xend,yst:del2:yend),...
% 0*reshape_x(xst:del2:xend,yst:del2:yend),...
% 2,"red");
% qu1=quiver3(reshape_x(xst:del2:xend,yst:del2:yend),reshape_y(xst:del2:xend,yst:del2:yend), ...
% 0*reshape_x(xst:del2:xend,yst:del2:yend)+100,...
% reshape_uc(xst:del2:xend,yst:del2:yend), ...
% reshape_vc(xst:del2:xend,yst:del2:yend),...
% 0*reshape_x(xst:del2:xend,yst:del2:yend),...
% 1.6,"red");
qu1=quiver3(reshape_x(xst:del2:xend,yst:del2:yend),reshape_y(xst:del2:xend,yst:del2:yend), ...
0*reshape_x(xst:del2:xend,yst:del2:yend)+100,...
reshape_ud(xst:del2:xend,yst:del2:yend), ...
reshape_vd(xst:del2:xend,yst:del2:yend),...
0*reshape_x(xst:del2:xend,yst:del2:yend),...
0.5,"red");
qu1.Color(4) = 0.1;

plot3(cy1(1000:11000,1),cy1(1000:11000,2),100+0*cy1(1000:11000,2),"black");
% surf(reshape_x,reshape_y,min(100,reshape_land)-100,'DisplayName','');shading interp;
surf(reshape_x(xst:xend,yst:yend),reshape_y(xst:xend,yst:yend),min(100,reshape_land(xst:xend,yst:yend)),'DisplayName','');shading interp;
cbar=colorbar();
cbar.Position(1)=0.923;
hold off
legend('flux',"cycle1",'FontSize',12)
%legend("flux","cycle1",'FontSize',12)
xlabel('PC1','FontSize',12);
ylabel('PC2','FontSize',12);
zlabel("");
axis([reshape_x(xst,yst),reshape_x(xend,yend),reshape_y(xst,yst),reshape_y(xend,yend)]);
set(h,'Units','Inches');
pos=get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
set(gcf,'position',[3,2,6.6,4.4]);
set(gcf,'PaperSize',[6.6,4.4]);