clear
close

real_I = [340991 612248 1032066 1461474 2100273 2817199 2442853 2058896 1459822 972307 589666 380604 267036 237754 175780 129371];


% parameters
d = 1000;
dt = 1 / d;
T = 20;
N = T * d;
total = 50000000;
min_loss = Inf;
good_beta = 0;
good_gamma = 0;

% 가중치 설정
w = zeros(1,20);
for x = 1:20 
    w(x) = log10(x) + 0.3;
end


% beta, gamma 탐색
for temp_beta = 0.1: 0.01 : 3                 
    for temp_gamma = 0.1 : 0.01 : 3
         
        S = zeros(1, N);
        I = zeros(1, N);
        R = zeros(1, N);
   
        S(1) = 0.983;
        I(1) = 0.004;
        R(1) = 0.013;

        for t = 1:N-1
            S(t) = 1 - I(t) - R(t);
            I(t+1) = I(t) + temp_beta * S(t) * I(t) *dt - temp_gamma * I(t) * dt ;
            R(t+1) = R(t) +  temp_gamma * I(t) * dt;
        end
        
        I = I * total;
        
        for i = 1 : T
            I(i) = I(i * d);
        end

        current_loss = 0;
        real_I_idx = 1;

% loss 계산
        for x = 1:size(real_I,2)                     
            current_loss = current_loss + (I(x) - real_I(real_I_idx))^2 * w(x);     
            real_I_idx = real_I_idx + 1;
        end

% loss 및 parameter 업데이트
        if min_loss > current_loss          
            min_loss = current_loss; 
            good_beta = temp_beta;
            good_gamma = temp_gamma;       
        end       
    end
end


S = zeros(1, N);
I = zeros(1, N);
R = zeros(1, N);

S(1) = 0.983;
I(1) = 0.004;
R(1) = 0.013;

for t = 1:N-1
    if t<13*d
        S(t) = 1 - I(t) - R(t);
        I(t+1) = I(t) + good_beta * S(t) * I(t) *dt - good_gamma * I(t) * dt ;
        R(t+1) = R(t) +  good_gamma * I(t) * dt;
    else
        S(t) = 1 - I(t) - R(t);
        I(t+1) = I(t) + good_beta*1.2 * S(t) * I(t) *dt - good_gamma * I(t) * dt ;
        R(t+1) = R(t) +  good_gamma * I(t) * dt;
    end
end

I = I * total;

for i = 1 : T
    I(i) = I(i * d);
end

good_beta
good_gamma


% plot
figure
plot(real_I, '-x');
hold on;

plot(I(1:20), '-o');
hold on;

xlabel('weeks')
ylabel('# of Infectors')
legend("real", "pred");
grid on