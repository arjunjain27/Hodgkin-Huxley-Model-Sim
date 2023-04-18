%%
% graphing the voltage-dependent functions a, B, n, tau for gating variable n
% using the steady state activation rate, opening/closing rate, and time constant function equations from HH

V = -120:0.1:0;
a_n = zeros(1, length(V));
B_n = zeros(1, length(V));
tau_n = zeros(1, length(V));
n_inf = zeros(1, length(V));

for i = 1:length(V)
    a_n(i) = (0.01*(V(i)+55))/(1-exp(-0.1*(V(i)+55)));
    B_n(i) = 0.125*exp(-0.0125*(V(i)+65));
    
    tau_n(i) = 1/(a_n(i)+B_n(i));
    n_inf(i) = tau_n(i)*a_n(i);
end

figure();
subplot(1, 3, 1)
hold on
plot(V, a_n)
plot(V, B_n)
hold off
xlabel("V (mV)")
ylabel("a or B (1/ms)")
legend({'a', 'B'})

subplot(1, 3, 2)
plot(V, n_inf)
xlabel("V (mV)")
ylabel("n (limiting value)")

subplot(1, 3, 3)
plot(V, tau_n)
xlabel("V (mV)")
ylabel("tau (ms)")

%%
% graphing the voltage-dependent functions n, m, h and voltage-dependent time constants tau_n, tau_m, tau_h
% using the rate function equations from HH

V = -120:0.1:20;

n_inf = zeros(1, length(V));
a_n = zeros(1, length(V));
tau_n = zeros(1, length(V));
B_n = zeros(1, length(V));

m_inf = zeros(1, length(V));
a_m = zeros(1, length(V));
tau_m = zeros(1, length(V));
B_m = zeros(1, length(V));

h_inf = zeros(1, length(V));
a_h = zeros(1, length(V));
tau_h = zeros(1, length(V));
B_h = zeros(1, length(V));

for i = 1:length(V)
    a_n(i) = (0.01*(V(i)+55))/(1-exp(-0.1*(V(i)+55)));
    B_n(i) = 0.125*exp(-0.0125*(V(i)+65));
    tau_n(i) = 1/(a_n(i)+B_n(i));
    n_inf(i) = tau_n(i)*a_n(i);
    
    a_m(i) = (0.1*(V(i)+40))/(1-exp(-0.1*(V(i)+40)));
    B_m(i) = 4*exp(-0.0556*(V(i)+65));
    tau_m(i) = 1/(a_m(i)+B_m(i));
    m_inf(i) = tau_m(i)*a_m(i);
    
    a_h(i) = 0.07*exp(-0.05*(V(i)+65));
    B_h(i) = 1/(1+exp(-0.1*(V(i)+35)));
    tau_h(i) = 1/(a_h(i)+B_h(i));
    h_inf(i) = tau_h(i)*a_h(i);
end

figure();
hold on
plot(V, n_inf)
plot(V, m_inf)
plot(V, h_inf)
hold off
xlabel("V (mv)")
legend({'n', 'm', 'h'})

figure();
hold on
plot(V, tau_n)
plot(V, tau_m)
plot(V, tau_h)
hold off
xlabel("V (mv)")
ylabel("tau (ms)")
legend({'n', 'm', 'h'})
