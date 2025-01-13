f_0 = 2.4 * 10^9;   % Base frequency (Hz)
v = 5:5:50;            % Vehicle speed relative to base station (m/s)   
x_val = -50:1:50;   %set up path for computation of doppler shift along it
f_x=zeros(1,length(x_val));
in_sim=Simulink.SimulationInput("wifi_1");
close all;
shift=zeros(1,length(f_x));
res_rms=zeros(1,length(v));
for j = 1:length(v) %iterate the test with different speed
    f_x(x_val<0)=v(j); %relative speed along the path
    f_x(x_val>0)=-(j);
    for i = 1:length(f_x)
        shift(i)=phase_shift(1/f_0,f_x(i),f_0); %calculte phase shift along the path
        if f_x(i)>0
            shift(i)=-shift(i); %inverse phase shift once we passed the origin (base station)
        end
    end
    in_sim=in_sim.setVariable('shift',shift);
    out = sim(in_sim); %run simulink sim with the previously calculated phase shift
    SER = out.yout.get('ErrorVec').Values.Data(:,1); %get SER data from the sim
    res_rms(j)=rms(SER); %calculate RMS of SER for each speed
end
plot(v,res_rms,'ro');
ylabel("SER");
xlabel("Speed (m/s)");
yline(mean(res_rms));
legend("RMS SER","Mean value")
%uncomment this part for a look at phase shift along the path
%close all
%figure(1);
%plot(x_val,shift);
%ylabel('Phase shift');
%xlabel("Distance from base station")

% Calculate Doppler shift delta in Hz given relative speed and base frequency
function delta_f = doppler_delta(v_relative, f_0)
    c = physconst('LightSpeed'); % Speed of light (m/s)
    delta_f = f_0 * ((v_relative) / c);                    
    %delta_f = f_0 * ((c + v_relative) / c) - f_0; % New frequency after doppler shift (Hz)
end

% Calculate phase shift as a result of Doppler shift
function delta_phi = phase_shift(t, v_relative, f_0)
    delta_phi = 2 * pi * (f_0 + doppler_delta(v_relative, f_0)) * t;
end