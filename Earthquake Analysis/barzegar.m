clc;
clear;
close all;

%Negar Sangari 950212443
%Mohammad Barzegar 950212402


[FileName,PathName] = uigetfile({'*.*','All Files (*.*)'}, ...
'Pick a File');

path=strcat(PathName,FileName);

Eq=load(path);

Eq=Eq';

Eq=Eq(:);

Eq=Eq*9.81;

Eq=Eq';

m=1;
e=input('Enter Damping ratio (0.05): ');


tf=input('Enter final T: ');
dt=input('Time step (0.01s for Elcentro): ');





type=input('Choose and type linear or constant: ', 's');

if (type=="constant") || (type=="Constant")
    gamma=0.5;
    beta=0.25;
else
    gamma=0.5;
    beta=1/6;
end


p=-m*Eq;

u=zeros(1,length(Eq));
ud=zeros(1,length(Eq));
udd=zeros(1,length(Eq));

u(1)=0;
ud(1)=0;
udd(1)=p(1);


it=(tf/dt);

Sa=zeros(1,400);    
Sd=zeros(1,400);    




for j=0.01:0.01:4
    
    T=0.01+j;
    o=2*3.14/T;
    k=m*(o^2);
    c=2*m*e*o;

    %udd=(p(1)-c*ud-k*u)/m;
    a1=(m/(beta*(dt^2)))+((gamma*c)/(beta*dt));
    a2=(m/(beta*dt))+((gamma/beta)-1)*c;
    a3=(((1/(2*beta))-1)*m)+dt*c*((gamma/(2*beta))-1);
    khat=k+a1;


        for i=2:it
        
        ph(i)=p(i)+a1*u(i-1)+a2*ud(i-1)+a3*udd(i-1);
    
        u(i)=ph(i)/khat;
    
        ud(i)=(((gamma)*(u(i)-u(i-1)))/(beta*dt))+(1-(gamma/beta))*ud(i-1)+(1-(gamma/(2*beta)))*dt*udd(i-1);
    
        udd(i)=((u(i)-u(i-1))/(beta*dt^2))-(ud(i-1)/(beta*dt))-((1/(2*beta))-1)*udd(i-1);
        
        end
        
u2=u(isfinite(u));
u3=max(abs(u2));
Sd(round(j*100))=u3;
Sa(round(j*100))=Sd(round(j*100))*o^2;       
        
end
 


time=0.01:0.01:4;
figure;
subplot(1,2,1);
plot(time,Sa,'linewidth',2);
xlabel('Period [sec]');
ylabel('Sa [M/s^2]');

subplot(1,2,2);
plot(time,Sd,'linewidth',2);
xlabel('Period [sec]');
ylabel('Sd [Meters]');
grid on;


