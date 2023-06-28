close all
clear all
clc

% Alínea a) -- Euler Explicito

% Definição de Variáveis (S.I.)

posicao_inicial=1;
velocidade_inicial=0;
K=1; % Constante da Mola
m=1;

tempo_inicial=0;
delta_t=0.01;
tempo_final=10;

t=tempo_inicial:delta_t:tempo_final;

N=length(t);

posicao=zeros(1,N);
posicao(1)=posicao_inicial;

velocidade=zeros(1,N);
velocidade(1)=velocidade_inicial;

for i=1:N-1
    
    aceleracao=-((K*posicao(i))/m);
    velocidade(i+1)=velocidade(i)+aceleracao*delta_t;
    posicao(i+1)=posicao(i)+velocidade(i)*delta_t;
    
end

figure(1)

plot(t,posicao,'-c')
hold on
plot(t,velocidade,'-r')
ylabel('Posição (m)')
xlabel('Tempo (s)')
title('a) EULER')
legend('Posicao/Tempo','Velocidade/Tempo')
grid on


% Alínea b) -- Euler-Cromer

clear all

% Definição de Variáveis (S.I.)

posicao_inicial=1;
velocidade_inicial=0;
K=1; % Constante da Mola
m=1;

tempo_inicial=0;
delta_t=0.01;
tempo_final=10;

t=tempo_inicial:delta_t:tempo_final;

N=length(t);

posicao=zeros(1,N);
posicao(1)=posicao_inicial;

velocidade=zeros(1,N);
velocidade(1)=velocidade_inicial;

for i=1:N-1
    
    aceleracao=-((K*posicao(i))/m);
    velocidade(i+1)=velocidade(i)+aceleracao*delta_t;
    posicao(i+1)=posicao(i)+velocidade(i+1)*delta_t;
    
end

figure(2)

plot(t,posicao,'-c')
hold on
plot(t,velocidade,'-r')
ylabel('Posição (m)')
xlabel('Tempo (s)')
title('b) EULER-CROMER')
legend('Posicao/Tempo','Velocidade/Tempo')
grid on


% Alínea C) -- Euler Implícito

clear all

% Definição de Variáveis (S.I.)

posicao_inicial=1;
velocidade_inicial=0;
K=1; % Constante da Mola
m=1;

tempo_inicial=0;
delta_t=0.01;
tempo_final=10;

t=tempo_inicial:delta_t:tempo_final;

N=length(t);

posicao=zeros(1,N);
posicao(1)=posicao_inicial;

velocidade=zeros(1,N);
velocidade(1)=velocidade_inicial;

% a*Z=b

a=[1,(K/m)*delta_t;-delta_t,1];

for i=1:N-1
    
%     posicao(i+1)=posicao(i)+velocidade(i+1)*delta_t;
%     aceleracao=-((K*posicao(i+1))/m);
%     velocidade(i+1)=velocidade(i)+aceleracao*delta_t;
    
%     Z=[velocidade(i+1);posicao(i+1)];

    b=[velocidade(i);posicao(i)];
    
    Z=linsolve(a,b);
    
    velocidade(i+1)=Z(1);
    posicao(i+1)=Z(2);
    
end

figure(3)

plot(t,posicao,'-c')
hold on
plot(t,velocidade,'-r')
ylabel('Posição (m)')
xlabel('Tempo (s)')
title('c) EULER IMPLÍCITO')
legend('Posicao/Tempo','Velocidade/Tempo')
grid on


% Alínea D) -- Crank-Nicolson

clear all

% Definição de Variáveis (S.I.)

posicao_inicial=1;
velocidade_inicial=0;
K=1; % Constante da Mola
m=1;

tempo_inicial=0;
delta_t=0.01;
tempo_final=10;

t=tempo_inicial:delta_t:tempo_final;

N=length(t);

posicao=zeros(1,N);
posicao(1)=posicao_inicial;

velocidade=zeros(1,N);
velocidade(1)=velocidade_inicial;

% a*Z=b

a=[1,(K/(2*m))*delta_t;-(delta_t/2),1];

for i=1:N-1
    
%     posicao(i+1)=posicao(i)+velocidade(i+1)*delta_t;
%     aceleracao=-((K*posicao(i+1))/m);
%     velocidade(i+1)=velocidade(i)+aceleracao*delta_t;
    
%     Z=[velocidade(i+1);posicao(i+1)];

    b=[-(K/(2*m))*delta_t*posicao(i)+velocidade(i);posicao(i)+velocidade(i)*(delta_t/2)];
    
    Z=linsolve(a,b);
    
    velocidade(i+1)=Z(1);
    posicao(i+1)=Z(2);
    
end

figure(4)

plot(t,posicao,'-c')
hold on
plot(t,velocidade,'-r')
ylabel('Posição (m)')
xlabel('Tempo (s)')
title('d) CRANK-NICOLSON')
legend('Posicao/Tempo','Velocidade/Tempo')
grid on

% Alínea E) usando Crank-Nicolson

aux=lagr();





