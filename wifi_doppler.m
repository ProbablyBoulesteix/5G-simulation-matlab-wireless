open_system('wifi_1');
f_0 = 2.4 * 10^9;   % Base frequency (Hz)
v = 5:1:50;            % Vehicle speed relative to base station (m/s)   
x_val = -50:1:50;   %set up path for computation of doppler shift along it
symbol_rate = 250e3;
SNR = 15;
f_x=zeros(1,length(x_val));
in_sim=Simulink.SimulationInput("wifi_1");
shift=zeros(1,length(f_x));
res_rms=zeros(1,length(v));
res_dataRate=zeros(1,length(v));
for j = 1:length(v) %iterate the test with different speed
    f_x(x_val<0)=v(j); %relative speed along the path
    f_x(x_val>0)=-v(j);
    for i = 1:length(f_x)
        shift(i)=phase_shift(1/f_0,f_x(i),f_0); %calculte phase shift along the path
        if f_x(i)>0
            shift(i)=-shift(i); %inverse phase shift once we passed the origin (base station)
        end
    end
    set_param('wifi_1/Constant','value',mat2str(shift(1,:))); %input the calculated phase shift into the Simulink model
    %get_param('wifi_1/Constant','Value');
    out = sim(in_sim); %run simulink sim with the previously calculated phase shift
    SER = out.yout.get('ErrorVec').Values.Data(:,1); %get SER data from the sim
    res_rms(j)=rms(SER); %calculate RMS of SER for each speed
    res_dataRate(j)=(52*symbol_rate*6*res_rms(j));% data rate in bits = symbol rate*number of bit per symbol * SER (52 because there are 52 channel in Wifi4)
end
close all;
figure(1)
subplot(3,1,1);
plot(v,res_rms,'ro');
ylabel("RMS SER");
xlabel("Speed (m/s)");
subplot(3,1,2);
plot(v,res_dataRate,'bo');
ylabel("Data rate (Mbit/s)");
xlabel("Speed (m/s)")
title("Average data rate over the path at different speed");
subplot(3,1,3)
histogram(res_dataRate,'Normalization','cdf');
xlabel("Data rate (Mbit/s)");
ylabel("CDF")
txt='Probability that the datarate is below or equal a certain value';
text(5e4,0.6,txt,'FontSize',8);
title("Data rate CDF");

% Calculate Doppler shift delta in Hz given relative speed and base frequency
function delta_f = doppler_delta(v_relative, f_0)
    c = physconst('LightSpeed'); % Speed of light (m/s)
    %delta_f = f_0 * ((v_relative) / c);                    
    delta_f = f_0 * ((c + v_relative) / c) - f_0; % New frequency after doppler shift (Hz)
end

% Calculate phase shift as a result of Doppler shift
function delta_phi = phase_shift(t, v_relative, f_0)
    delta_phi = 2 * pi * (f_0 + doppler_delta(v_relative, f_0)) * t;
end