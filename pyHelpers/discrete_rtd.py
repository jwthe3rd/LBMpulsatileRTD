import numpy as np
import matplotlib.pyplot as plt
import math

num_particles = 5000
num_pulses = 9
deltaT = 0.005
data_file = "run_log"
#comp_data_file = "dataFVM/2ml.dat"
plt.rcParams['figure.figsize'] = [10, 7]
plt.rcParams['figure.dpi'] = 200
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.size'] = 16.0
def javad_velocity(iT):

    t = ((iT)*1.111-math.floor((iT)*1.111))/1.111

    if (t <= 0.15):
      vel_char =1.322515215*0.722378598771807*((0.349999999999992+9.39592352092547*t-191.515151515187*(t**2))+871.933621933773*((t**3)))
    else:
      vel_char =1.322515215*0.722378598771807*((-8.87872100015116+158.853244453899*t-1009.813663592*((t**2))+3083.59699791101*((t**3))-4840.042878547*((t**4))+3765.53701950272*((t**5))-1142.95591854388*((t**6))-8.95498219229757*((t**7))))

    return 0.5*vel_char;

with open(data_file,"r") as f:

    lines = f.readlines()

stats_line = []
time_offset_line = []

for line in lines:
    line_split = line.split()

    if line_split[0] == '[captureStatistics]' and len(line_split) == 4:
        stats_line.append(int(line_split[3].split('=')[1]))
    elif line_split[0] == '[LatticeStatistics]' and len(line_split) == 6:
        if float(line_split[2][2]) > 0:
            time_offset_line.append(float(line_split[2][2])-1)

e_t = np.zeros(shape=np.asarray(stats_line).shape)

for i, line in enumerate(stats_line):
    if i == 0:
        summer = 0
    else:
        summer+= line
        e_t[i] = ((line - stats_line[i-1])/(num_particles*deltaT))

t = np.arange(0, deltaT*len(e_t)-0.5*deltaT, deltaT)
time_offset = np.asarray(time_offset_line)
#t = np.arange(0, deltaT*len(e_t), deltaT)
print(t)
print(e_t)

# with open(comp_data_file, "r") as f:

#     lines = f.readlines()

# C_t_comp = []
# delta_Ts = []

# for line in lines:
#     C_t_comp.append(float(line.split()[1]))
#     delta_Ts.append(float(line.split()[2]))

# C_t_comp = np.asarray(C_t_comp)
# e_t_comp = C_t_comp / np.trapz(C_t_comp, delta_Ts)

# e_t_comp = np.zeros(shape=C_t_comp.shape)

# for i, line in enumerate(C_t_comp):
#     if i!= 0:
#         e_t_comp[i] = ((line - C_t_comp[i-1])/(delta_Ts[i] - delta_Ts[i-1]))

print(f"Area under the E curve  -- > {np.sum(e_t*deltaT)}")
print(f"Mean Residence Time LBM -- > {np.trapz(t*e_t, t, deltaT)}")
# print(f"Mean Residence Time FVM -- > {np.trapz(delta_Ts*e_t_comp, delta_Ts)}")
# print(f"Space time (V/Q)  -- > {1000000*0.02*(np.pi*(0.005**2)) / 1.985}")


pulse_separated_t = []
pulse_separated_e_t = []
mean_rts = []

plt.figure()

for i in range(num_pulses):
    pulse_separated_t.append((t-time_offset)[i*200 : (i+1)*200])
    pulse_separated_e_t.append(e_t[i*200 : (i+1)*200])

pulse_separated_t = np.asarray(pulse_separated_t)
pulse_separated_e_t = np.asarray(pulse_separated_e_t)

for i in range(num_pulses):
    plt.plot(pulse_separated_t[i],pulse_separated_e_t[i], label=f'Pulse #{i+1}')
    mean_rts.append(np.trapz(pulse_separated_t*pulse_separated_e_t, pulse_separated_t, deltaT))

plt.xlabel("Time (s)")
plt.ylabel("E(t)")
plt.legend()
plt.savefig("layered_rtd.png")   

vels = np.zeros(shape=t.shape)

for i, iT in enumerate(t):
    vels[i] = javad_velocity(iT)


prev_val = -1

for idx, val in enumerate(time_offset):

    if val == prev_val:
        continue
    else:
        print(f"Velocity of pulse {int(idx/200)+1} is {vels[idx]}")
        prev_val = val


plt.figure()

fig, axs = plt.subplots(2)
fig.supxlabel('Time (s)')
axs[0].plot(t, e_t)
for i in range(num_pulses):
    axs[0].plot([i, i],[min(e_t), max(e_t)], 'r-')
axs[1].plot(t, vels)
for i in range(num_pulses):
    axs[1].plot([i, i],[min(vels), max(vels)], 'r-')
axs[0].set(ylabel='E(t)')
axs[1].set(ylabel='Velocity (m/s)')
plt.savefig('vel_plot.png')