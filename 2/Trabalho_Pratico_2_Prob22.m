close all
clear all
clc

% Alínea a)

% Definição de Variáveis

GMs=4*(pi^2);

tempo_inicial=0;
delta_t=0.0001;
tempo_final=1;

t=tempo_inicial:delta_t:tempo_final;

% Começamos no sitio mais afastado do Sol
posicao_inicial_X=0.47;
posicao_inicial_Y=0;

velocidade_inicial_X=0;
velocidade_inicial_Y=8.2;

N=length(t);

posicao_X=zeros(1,N);
posicao_X(1)=posicao_inicial_X;
posicao_Y=zeros(1,N);
posicao_Y(1)=posicao_inicial_Y;

velocidade_X=zeros(1,N);
velocidade_X(1)=velocidade_inicial_X;
velocidade_Y=zeros(1,N);
velocidade_Y(1)=velocidade_inicial_Y;


for i=1:N-1
    
    r=sqrt(posicao_X(i)^2 + posicao_Y(i)^2);
    
    aceleracao_X=-((GMs*posicao_X(i))/(r^3));
    aceleracao_Y=-((GMs*posicao_Y(i))/(r^3));
    
    velocidade_X(i+1)=velocidade_X(i) + aceleracao_X*delta_t;
    velocidade_Y(i+1)=velocidade_Y(i) + aceleracao_Y*delta_t;
    
    posicao_X(i+1)=posicao_X(i)+velocidade_X(i+1)*delta_t;
    posicao_Y(i+1)=posicao_Y(i)+velocidade_Y(i+1)*delta_t;
    
end

figure(1)

plot(posicao_X,posicao_Y,'-k')
hold on
plot(0,0,'oc','LineWidth',3)
title('Órbita de Mercúrio')
legend('Órb. Mercurio','Sol')
grid on
axis([-0.5 0.5 -0.5 0.5])
set(gca,'PlotBoxAspectRatio',[1 1 1])

figure(2)

plot(t,posicao_X)
title('Período da Órbita de Mercúrio')
grid on
ylabel('Posição (m)')
xlabel('Tempo (s)')


% Alínea b)

% Vamos identificar os indices que sao máximos, vamos ajustar uma 
%parabola aos 3 pontos de cima da função do periodo

t_max=[];
% posicao_X_max=[];

for i=2:N-1
    
    if posicao_X(i) > posicao_X(i-1) && posicao_X(i) > posicao_X(i+1)
        
        aux=lagr(t(i-1:i+1),posicao_X(i-1:i+1));
        t_max=[t_max, aux(1)];
%         posicao_X_max=[posicao_X_max, aux(2)];
        
    end
    
end

N_max=length(t_max);
T_med=0; % ---> período médio

for i=1:N_max-1
    
    T_med=T_med +t_max(i+1)-t_max(i);
    
end

T_med = T_med/(N_max-1);

% OU

T=((t_max(N_max)-t_max(1))/(N_max-1)); % ---> Período

fprintf('O período da Orbita de Mercúrio em torno do Sol é de %f (rad/s).\n',T)

