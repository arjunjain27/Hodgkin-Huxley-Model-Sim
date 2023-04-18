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
