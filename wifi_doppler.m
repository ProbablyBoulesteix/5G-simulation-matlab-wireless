f_0 = 2.4 * 10^9;   % Base frequency (Hz)
v = 100;            % Vehicle speed relative to base station (m/s)   
x_val = -50:1:50;
f_x=zeros(1,length(x_val));
f_x(x_val<0)=-v;
f_x(x_val>0)=v;
shift=zeros(1,length(f_x));
for i = 1:length(f_x)
    dopp(i)=doppler_delta(f_x(i),f_0);
    shift(i)=phase_shift(1/f_0,f_x(i),f_0);
end
max(shift)-min(shift) %TODO calculate by hand if we're supposed to have such small results, if yes : we're fucked
close all
figure(1);
subplot(2,1,1);
plot(x_val,dopp);
ylabel('Doppler shift');
subplot(2,1,2);
plot(x_val,shift);
ylabel('Phase shift');

% Calculate Doppler shift delta in Hz given relative speed and base frequency
function delta_f = doppler_delta(v_relative, f_0)
    c = physconst('LightSpeed');
    %delta_f = f_0 * ((v_relative) / c);                    % Speed of light (m/s)
    delta_f = f_0 * ((c + v_relative) / c) - f_0;   % New frequency after doppler shift (Hz)
end

% Calculate phase shift as a result of Doppler shift
function delta_phi = phase_shift(t, v_relative, f_0)
    delta_phi = 2 * pi * (f_0 + doppler_delta(v_relative, f_0)) * t;
end

%fplot(@(t) phase_shift(t, v, f_0), [-100, 100]); 