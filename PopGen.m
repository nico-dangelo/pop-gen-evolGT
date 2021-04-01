% %s= selection strength
 s=0.01;
w0=[1.02,1.03,0.99;1.00,1.01,0.98;1.01,1.04,1.02];
 p0=[0.1111,0.1111,0.1111;0.1111,0.1111,0.1111;0.1111,0.1111,0.1111];
 gens=100;
 p=p0;
w=w0;
wbar=sum(p0,'all');
x=[0.3333,0.3333,0.3333];
y=[0.3333;0.3333;0.3333];
for t=2:gens+1
p(:,:,t)=p(:,:,t-1).*w;
wbar(t)=sum(p(:,:,t),'all'); %sum of elements to normalize frequencies and ensure probability distribution 
p(:,:,t)=p(:,:,t)./wbar(t);%normalize frequencies
x(:,:,t)=sum(p(:,:,t-1),1);
y(:,:,t)=sum(p(:,:,t-1),2);
wbarnew(t)=sum(p(:,:,t),'all'); %Check that probabilities still add to 1
sumx(t)=sum(x(:,:,t));
sumy(t)=sum((y(:,:,t)));     
end
% for t2=1:gens
% for i=1:3
% for j=1:3
% D(i,j,t)=p(i,j,t)-(x(:,i,t)*y(j,:,t));
% end
% end
% end
%disp(p(:,:,gens+1));

% for j=1:t
% 
% A=[1,2,3];
% B=[1;2;3];
% surf(A,B,p(:,:,j));
% 
% colormap hot;
% colorbar;
% xlabel('Locus B');
% ylabel('Locus A')
% zlabel('Genotype frequencies');
% F(j)=getframe(gcf);
%  end

    
   
  

%%%%%%%%%%%%%%%%%%%%%%%% Game Theoretic 
%Payoff matrix is differential fitness landscape
Q=(w0-1)/s;
for t3=2:gens+1
xg(:,:,t3)=x(:,:,t3-1)+x(:,:,t3-1)*s*sum((y(:,:,t3-1).*Q),'all');
qbar1(t3)=sum(x(:,:,t3),'all');
qbar2(t3)=sum(y(:,:,t3),'all');
xg(:,:,t3)=xg(:,:,t3)./qbar1(t3);
yg(:,:,t3)=y(:,:,t3-1)+sum((s*(x(:,:,t3-1).*Q)),'all')*y(:,:,t3-1);
yg(:,:,t3)=yg(:,:,t3)./qbar2(t3);
end
% v = VideoWriter('freqs_strat.avi');
% open(v);
% for k=1:gens
% subplot(3,1,1);
% plot([1,2,3],yg(:,:,k));
% xlim([0 4]);
% xlabel('Allele');
% ylabel('Strategy probability');
% title("Mixed strategy at locus A after " + k + " generations");
% axis equal
% subplot(3,1,2);
% plot([1,2,3],xg(:,:,k));
% xlabel('Allele');
% ylabel('Strategy probability');
% title("Mixed strategy at locus B after " + k + " generations");
% xlim([0 4]);
% axis equal;
% subplot(3,1,3);
% surf([1,2,3],[1,2,3],p(:,:,k));
% colormap hot;
% colorbar;
% xlabel('Locus B');
% ylabel('Locus A')
% zlabel('Frequency');
% title("Genotype frequencies after " + k +" generations");
% 
% F(k)=getframe(gcf);
% writeVideo(v,F(k));
% end
% close(v);

% surf([1,2,3],[1,2,3],p(:,:,101));
% colormap hot;
% colorbar;
% xlabel('Locus B');
% ylabel('Locus A')
% zlabel('Frequency');
% title("Genotype frequencies after " + 100 +" generations");
% saveas(gcf,'Genfreq100.jpg');