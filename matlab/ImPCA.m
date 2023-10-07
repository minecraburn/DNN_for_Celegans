copyer1=WT_NoStim(1).traces(:,d1(1:1,1:dim)); %%%%%%%%%
datasize=zeros(6,1);
datasize(1)=0;
temp=size(copyer1);
datasize(2)=temp(1);
for k=2:1:n0
 copyer1=[copyer1;WT_NoStim(k).traces(:,d1(k:k,1:dim))];%%%%%%%%%%%%%%%%%%%
 temp=size(copyer1);
 datasize(k+1)=temp(1);
end
[n,~]=size(copyer1);
s=max(copyer1)-min(copyer1);
gj=(copyer1-repmat(min(copyer1),n,1))./repmat(s,n,1);

mysee=1;
winl=1;winr=3137;line_num=0;
% winl=274;
% line_num=0;
% winr=1350;
[coeff,score,latent]=pca(gj); 
cor=cumsum(latent)/sum(latent);
x1=score(datasize(mysee)+1+winl:datasize(mysee+line_num)+winr,1)';
x2=score(datasize(mysee)+1+winl:datasize(mysee+line_num)+winr,2)';
x3=score(datasize(mysee)+1+winl:datasize(mysee+line_num)+winr,3)';
% x1=copyer1(1:3137,1)';
% x2=copyer1(1:3137,4)';
% x3=copyer1(1:3137,3)';
plot3(x1,x2,x3);
xlabel("PC1");
ylabel("PC2");
zlabel("PC3");
figure(2);
plot(x1,x2);
xlabel("PC1");
ylabel("PC2");
invcoe=inv(coeff);

cc2=min(copyer1);
cc3=s;
www=copyer1-score*invcoe.*cc3+cc2;
cc1=www(1:1,:);
cc0=cc1-cc2;

traA=invcoe.*cc3;
traB=coeff./cc3';
temp11=copyer1-score*traA-cc0;
temp22=(copyer1-cc0)*traB-score;
fpsarr=zeros(5,1);
for k=1:1:n0
 fpsarr(k)=WT_NoStim(k).fps;%%%%%%%%%%%%%%%%%%%
end