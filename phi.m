function f = phi(mu,sigma)
tau_s = 0.0005; % s
tau_r = 0.002; % s
tau_m = 0.01; % s
V_leak = -65; % mV
V_reset = -65; % mV
V_thresh = -50; % mV
gamma = abs(zeta(0.5))/sqrt(2);
s_tau = sqrt(tau_s/tau_m);
N = numel(mu);
f = zeros(N,1);
for n=1:N
    lower = (V_reset - V_leak - mu(n)) / sigma(n) + gamma * s_tau;
    upper = (V_thresh -V_leak - mu(n)) / sigma(n) + gamma * s_tau;
    x = linspace(lower,upper,100);
    y = exp(x.^2) .* (1 + erf(x));
    f(n) = (tau_r + tau_m * sqrt(pi) * trapz(x,y))^-1; % Hz
end
