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

%%
% full simulation of Hodgkin-Huxley model using standard parameters (Dayan and Abbott)
% simulating a short injection current

c_m = 0.01; %nF/mm^2
A = 0.1; %mm^2

%maximal conductances (mS/mm^2)
g_L = 0.003;
g_K = 0.36;
g_Na = 1.2;

%reversal potentials (mV)
E_L = -54.387;
E_K = -77;
E_Na = 50;

dt = 0.0001; %(ms)
time = 0:dt:15;

V = zeros(1, length(time));
V(1) = -65; %(mV)

i_m = zeros(1, length(V));

%initializing gating variables
n = zeros(1, length(V));
a_n = zeros(1, length(V));
B_n = zeros(1, length(V));
tau_n = zeros(1, length(V));
n_inf = zeros(1, length(V));

m = zeros(1, length(V));
a_m = zeros(1, length(V));
B_m = zeros(1, length(V));
tau_m = zeros(1, length(V));
m_inf = zeros(1, length(V));

h = zeros(1, length(V));
a_h = zeros(1, length(V));
B_h = zeros(1, length(V));
tau_h = zeros(1, length(V));
h_inf = zeros(1, length(V));

a_n(1) = (0.01*(V(1)+55))/(1-exp(-0.1*(V(1)+55)));
B_n(1) = 0.125*exp(-0.0125*(V(1)+65));
tau_n(1) = 1/(a_n(1)+B_n(1));
n_inf(1) = tau_n(1)*a_n(1);
n(1) = n_inf(1);
    
a_m(1) = (0.1*(V(1)+40))/(1-exp(-0.1*(V(1)+40)));
B_m(1) = 4*exp(-0.0556*(V(1)+65));
tau_m(1) = 1/(a_m(1)+B_m(1));
m_inf(1) = tau_m(1)*a_m(1);
m(1) = m_inf(1);
    
a_h(1) = 0.07*exp(-0.05*(V(1)+65));
B_h(1) = 1/(1+exp(-0.1*(V(1)+35)));
tau_h(1) = 1/(a_h(1)+B_h(1));
h_inf(1) = tau_h(1)*a_h(1);
h(1) = h_inf(1);

%injection current
start = 5;
fin = 8;
I_e = zeros(1, length(V));
I_e(start/dt:fin/dt) = 0.015;

for i = 1:length(time)-1
    i_m(i) = g_L*(V(i)-E_L) + g_K*n(i)^4*(V(i)-E_K) + g_Na*m(i)^3*h(i)*(V(i)-E_Na);
    dV = 1/c_m * (I_e(i)/A - i_m(i)) * dt;
    
    V(i+1) = V(i) + dV;
    
    a_n(i+1) = (0.01*(V(i+1)+55))/(1-exp(-0.1*(V(i+1)+55)));
    B_n(i+1) = 0.125*exp(-0.0125*(V(i+1)+65));
    tau_n(i+1) = 1/(a_n(i+1)+B_n(i+1));
    n_inf(i+1) = tau_n(i+1)*a_n(i+1);

    a_m(i+1) = (0.1*(V(i+1)+40))/(1-exp(-0.1*(V(i+1)+40)));
    B_m(i+1) = 4*exp(-0.0556*(V(i+1)+65));
    tau_m(i+1) = 1/(a_m(i+1)+B_m(i+1));
    m_inf(i+1) = tau_m(i+1)*a_m(i+1);

    a_h(i+1) = 0.07*exp(-0.05*(V(i+1)+65));
    B_h(i+1) = 1/(1+exp(-0.1*(V(i+1)+35)));
    tau_h(i+1) = 1/(a_h(i+1)+B_h(i+1));
    h_inf(i+1) = tau_h(i+1)*a_h(i+1);
    
    dn = (a_n(i)*(1-n(i)) - B_n(i)*n(i)) * dt;
    dm = (a_m(i)*(1-m(i)) - B_m(i)*m(i)) * dt;
    dh = (a_h(i)*(1-h(i)) - B_h(i)*h(i)) * dt;
    
    n(i+1) = n(i)+dn;
    m(i+1) = m(i)+dm;
    h(i+1) = h(i)+dh;
end

figure();
subplot(5, 1, 1)
plot(time, V)
ylabel("V (mV)")

subplot(5, 1, 2)
plot(time, i_m)
ylabel("i_m (uA/mm^2)")

subplot(5, 1, 3)
plot(time, m)
ylabel("m")
ylim([0,1])

subplot(5, 1, 4)
plot(time, h)
ylabel("h")
ylim([0,1])

subplot(5, 1, 5)
plot(time, n)
ylabel("n")
xlabel("Time (ms)")
ylim([0,1])

%%
% full simulation of HH model using standard parameters but blocking sodium channels

V = zeros(1, length(time));
V(1) = -65; %(mV)

i_m = zeros(1, length(V));

%initializing gating variables
n = zeros(1, length(V));
a_n = zeros(1, length(V));
B_n = zeros(1, length(V));
tau_n = zeros(1, length(V));
n_inf = zeros(1, length(V));

m = zeros(1, length(V));
a_m = zeros(1, length(V));
B_m = zeros(1, length(V));
tau_m = zeros(1, length(V));
m_inf = zeros(1, length(V));

h = zeros(1, length(V));
a_h = zeros(1, length(V));
B_h = zeros(1, length(V));
tau_h = zeros(1, length(V));
h_inf = zeros(1, length(V));

a_n(1) = (0.01*(V(1)+55))/(1-exp(-0.1*(V(1)+55)));
B_n(1) = 0.125*exp(-0.0125*(V(1)+65));
tau_n(1) = 1/(a_n(1)+B_n(1));
n_inf(1) = tau_n(1)*a_n(1);
n(1) = n_inf(1);
    
a_m(1) = (0.1*(V(1)+40))/(1-exp(-0.1*(V(1)+40)));
B_m(1) = 4*exp(-0.0556*(V(1)+65));
tau_m(1) = 1/(a_m(1)+B_m(1));
m_inf(1) = tau_m(1)*a_m(1);
m(1) = m_inf(1);
    
a_h(1) = 0.07*exp(-0.05*(V(1)+65));
B_h(1) = 1/(1+exp(-0.1*(V(1)+35)));
tau_h(1) = 1/(a_h(1)+B_h(1));
h_inf(1) = tau_h(1)*a_h(1);
h(1) = h_inf(1);

%injection current
start = 5;
fin = 8;
I_e = zeros(1, length(V));
I_e(start/dt:fin/dt) = 0.015;

for i = 1:length(time)-1
    i_m(i) = g_L*(V(i)-E_L) + g_K*n(i)^4*(V(i)-E_K) + g_Na/10*m(i)^3*h(i)*(V(i)-E_Na);
    dV = 1/c_m * (I_e(i)/A - i_m(i)) * dt;
    
    V(i+1) = V(i) + dV;
    
    a_n(i+1) = (0.01*(V(i+1)+55))/(1-exp(-0.1*(V(i+1)+55)));
    B_n(i+1) = 0.125*exp(-0.0125*(V(i+1)+65));
    tau_n(i+1) = 1/(a_n(i+1)+B_n(i+1));
    n_inf(i+1) = tau_n(i+1)*a_n(i+1);

    a_m(i+1) = (0.1*(V(i+1)+40))/(1-exp(-0.1*(V(i+1)+40)));
    B_m(i+1) = 4*exp(-0.0556*(V(i+1)+65));
    tau_m(i+1) = 1/(a_m(i+1)+B_m(i+1));
    m_inf(i+1) = tau_m(i+1)*a_m(i+1);

    a_h(i+1) = 0.07*exp(-0.05*(V(i+1)+65));
    B_h(i+1) = 1/(1+exp(-0.1*(V(i+1)+35)));
    tau_h(i+1) = 1/(a_h(i+1)+B_h(i+1));
    h_inf(i+1) = tau_h(i+1)*a_h(i+1);
    
    dn = (a_n(i)*(1-n(i)) - B_n(i)*n(i)) * dt;
    dm = (a_m(i)*(1-m(i)) - B_m(i)*m(i)) * dt;
    dh = (a_h(i)*(1-h(i)) - B_h(i)*h(i)) * dt;
    
    n(i+1) = n(i)+dn;
    m(i+1) = m(i)+dm;
    h(i+1) = h(i)+dh;
end

figure();
subplot(5, 1, 1)
plot(time, V)
ylabel("V (mV)")
title("Membrane Potential with Blocked Na Channels");

subplot(5, 1, 2)
plot(time, i_m)
ylabel("i_m (uA/mm^2)")

subplot(5, 1, 3)
plot(time, m)
ylabel("m")
ylim([0,1])

subplot(5, 1, 4)
plot(time, h)
ylabel("h")
ylim([0,1])

subplot(5, 1, 5)
plot(time, n)
ylabel("n")
xlabel("Time (ms)")
ylim([0,1])

%%
% full simulation of HH model using standard parameters but blocking sodium channels

V = zeros(1, length(time));
V(1) = -65; %(mV)

i_m = zeros(1, length(V));

%initializing gating variables
n = zeros(1, length(V));
a_n = zeros(1, length(V));
B_n = zeros(1, length(V));
tau_n = zeros(1, length(V));
n_inf = zeros(1, length(V));

m = zeros(1, length(V));
a_m = zeros(1, length(V));
B_m = zeros(1, length(V));
tau_m = zeros(1, length(V));
m_inf = zeros(1, length(V));

h = zeros(1, length(V));
a_h = zeros(1, length(V));
B_h = zeros(1, length(V));
tau_h = zeros(1, length(V));
h_inf = zeros(1, length(V));

a_n(1) = (0.01*(V(1)+55))/(1-exp(-0.1*(V(1)+55)));
B_n(1) = 0.125*exp(-0.0125*(V(1)+65));
tau_n(1) = 1/(a_n(1)+B_n(1));
n_inf(1) = tau_n(1)*a_n(1);
n(1) = n_inf(1);
    
a_m(1) = (0.1*(V(1)+40))/(1-exp(-0.1*(V(1)+40)));
B_m(1) = 4*exp(-0.0556*(V(1)+65));
tau_m(1) = 1/(a_m(1)+B_m(1));
m_inf(1) = tau_m(1)*a_m(1);
m(1) = m_inf(1);
    
a_h(1) = 0.07*exp(-0.05*(V(1)+65));
B_h(1) = 1/(1+exp(-0.1*(V(1)+35)));
tau_h(1) = 1/(a_h(1)+B_h(1));
h_inf(1) = tau_h(1)*a_h(1);
h(1) = h_inf(1);

%injection current
start = 5;
fin = 8;
I_e = zeros(1, length(V));
I_e(start/dt:fin/dt) = 0.015;

for i = 1:length(time)-1
    i_m(i) = g_L*(V(i)-E_L) + g_K/10*n(i)^4*(V(i)-E_K) + g_Na*m(i)^3*h(i)*(V(i)-E_Na);
    dV = 1/c_m * (I_e(i)/A - i_m(i)) * dt;
    
    V(i+1) = V(i) + dV;
    
    a_n(i+1) = (0.01*(V(i+1)+55))/(1-exp(-0.1*(V(i+1)+55)));
    B_n(i+1) = 0.125*exp(-0.0125*(V(i+1)+65));
    tau_n(i+1) = 1/(a_n(i+1)+B_n(i+1));
    n_inf(i+1) = tau_n(i+1)*a_n(i+1);

    a_m(i+1) = (0.1*(V(i+1)+40))/(1-exp(-0.1*(V(i+1)+40)));
    B_m(i+1) = 4*exp(-0.0556*(V(i+1)+65));
    tau_m(i+1) = 1/(a_m(i+1)+B_m(i+1));
    m_inf(i+1) = tau_m(i+1)*a_m(i+1);

    a_h(i+1) = 0.07*exp(-0.05*(V(i+1)+65));
    B_h(i+1) = 1/(1+exp(-0.1*(V(i+1)+35)));
    tau_h(i+1) = 1/(a_h(i+1)+B_h(i+1));
    h_inf(i+1) = tau_h(i+1)*a_h(i+1);
    
    dn = (a_n(i)*(1-n(i)) - B_n(i)*n(i)) * dt;
    dm = (a_m(i)*(1-m(i)) - B_m(i)*m(i)) * dt;
    dh = (a_h(i)*(1-h(i)) - B_h(i)*h(i)) * dt;
    
    n(i+1) = n(i)+dn;
    m(i+1) = m(i)+dm;
    h(i+1) = h(i)+dh;
end

figure();
subplot(5, 1, 1)
plot(time, V)
ylabel("V (mV)")
title("Membrane Potential with Blocked K+ Channels");

subplot(5, 1, 2)
plot(time, i_m)
ylabel("i_m (uA/mm^2)")

subplot(5, 1, 3)
plot(time, m)
ylabel("m")
ylim([0,1])

subplot(5, 1, 4)
plot(time, h)
ylabel("h")
ylim([0,1])

subplot(5, 1, 5)
plot(time, n)
ylabel("n")
xlabel("Time (ms)")
ylim([0,1])

%%
% full simulation of HH model using standard parameters but with persistent sodium channels

V = zeros(1, length(time));
V(1) = -65; %(mV)

i_m = zeros(1, length(V));

%initializing gating variables
n = zeros(1, length(V));
a_n = zeros(1, length(V));
B_n = zeros(1, length(V));
tau_n = zeros(1, length(V));
n_inf = zeros(1, length(V));

m = zeros(1, length(V));
a_m = zeros(1, length(V));
B_m = zeros(1, length(V));
tau_m = zeros(1, length(V));
m_inf = zeros(1, length(V));

h = zeros(1, length(V));
a_h = zeros(1, length(V));
B_h = zeros(1, length(V));
tau_h = zeros(1, length(V));
h_inf = zeros(1, length(V));

a_n(1) = (0.01*(V(1)+55))/(1-exp(-0.1*(V(1)+55)));
B_n(1) = 0.125*exp(-0.0125*(V(1)+65));
tau_n(1) = 1/(a_n(1)+B_n(1));
n_inf(1) = tau_n(1)*a_n(1);
n(1) = n_inf(1);
    
a_m(1) = (0.1*(V(1)+40))/(1-exp(-0.1*(V(1)+40)));
B_m(1) = 4*exp(-0.0556*(V(1)+65));
tau_m(1) = 1/(a_m(1)+B_m(1));
m_inf(1) = tau_m(1)*a_m(1);
m(1) = m_inf(1);
    
a_h(1) = 0.07*exp(-0.05*(V(1)+65));
B_h(1) = 1/(1+exp(-0.1*(V(1)+35)));
tau_h(1) = 1/(a_h(1)+B_h(1));
h_inf(1) = tau_h(1)*a_h(1);
h(1) = h_inf(1);

%injection current
start = 5;
fin = 8;
I_e = zeros(1, length(V));
I_e(start/dt:fin/dt) = 0.015;

for i = 1:length(time)-1
    i_m(i) = g_L*(V(i)-E_L) + g_K*n(i)^4*(V(i)-E_K) + g_Na*m(i)^4*(V(i)-E_Na);
    dV = 1/c_m * (I_e(i)/A - i_m(i)) * dt;
    
    V(i+1) = V(i) + dV;
    
    a_n(i+1) = (0.01*(V(i+1)+55))/(1-exp(-0.1*(V(i+1)+55)));
    B_n(i+1) = 0.125*exp(-0.0125*(V(i+1)+65));
    tau_n(i+1) = 1/(a_n(i+1)+B_n(i+1));
    n_inf(i+1) = tau_n(i+1)*a_n(i+1);

    a_m(i+1) = (0.1*(V(i+1)+40))/(1-exp(-0.1*(V(i+1)+40)));
    B_m(i+1) = 4*exp(-0.0556*(V(i+1)+65));
    tau_m(i+1) = 1/(a_m(i+1)+B_m(i+1));
    m_inf(i+1) = tau_m(i+1)*a_m(i+1);

    a_h(i+1) = 0.07*exp(-0.05*(V(i+1)+65));
    B_h(i+1) = 1/(1+exp(-0.1*(V(i+1)+35)));
    tau_h(i+1) = 1/(a_h(i+1)+B_h(i+1));
    h_inf(i+1) = tau_h(i+1)*a_h(i+1);
    
    dn = (a_n(i)*(1-n(i)) - B_n(i)*n(i)) * dt;
    dm = (a_m(i)*(1-m(i)) - B_m(i)*m(i)) * dt;
    dh = (a_h(i)*(1-h(i)) - B_h(i)*h(i)) * dt;
    
    n(i+1) = n(i)+dn;
    m(i+1) = m(i)+dm;
    h(i+1) = h(i)+dh;
end

figure();
subplot(5, 1, 1)
plot(time, V)
ylabel("V (mV)")
title("Membrane Potential with Persistent Na+ Channels");

subplot(5, 1, 2)
plot(time, i_m)
ylabel("i_m (uA/mm^2)")

subplot(5, 1, 3)
plot(time, m)
ylabel("m")
ylim([0,1])

subplot(5, 1, 4)
plot(time, h)
ylabel("h")
ylim([0,1])

subplot(5, 1, 5)
plot(time, n)
ylabel("n")
xlabel("Time (ms)")
ylim([0,1])
