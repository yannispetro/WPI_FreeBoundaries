function [ tt,dt,N ] = Dim_time( Fs,Tot )
% Time vector (*w_n due to nondimensionalization tau -> w_n*t)

%     Fs = 100;                % Sampling frequency  
%     Tot  = 20;               % Total time  (sec)
    dt = 1/Fs;               % Sampling period (Original - sec)

    T = Tot;             % Total time
    
    N  = floor(T/dt);      % Length of signal
    tt = (0:N-1)*dt;      % Time vector

end

