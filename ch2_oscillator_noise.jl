using GNSSTools
using PyPlot
pygui(true)

t_length = 1.  # second
f_s = 5e6  # Hz
v_0 = L1_freq  # Hz
h₋₂, h₋₁, h₀, h₁, h₂ = h_parms_tcxo[4]

noise_h₋₂ = generate_phase_noise(t_length, f_s, v_0, [h₋₂, 0,  0,  0,  0])
noise_h₋₁ = generate_phase_noise(t_length, f_s, v_0, [0, h₋₁,  0,  0,  0])
noise_h₀  = generate_phase_noise(t_length, f_s, v_0, [0,   0, h₀,  0,  0])
noise_h₁  = generate_phase_noise(t_length, f_s, v_0, [0,   0,  0, h₁,  0])
noise_h₂  = generate_phase_noise(t_length, f_s, v_0, [0,   0,  0,  0, h₂])

t = range(0, t_length, length=length(noise_h₋₂))
title_font_size = 10

fig = figure(figsize=(6.5,1.2))
ax1 = fig.add_subplot(1,5,1)
ax1.plot(t, real.(noise_h₋₂), "k")
title("h₋₂", fontsize=title_font_size)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax2 = fig.add_subplot(1,5,2)
ax2.plot(t, real.(noise_h₋₁), "k")
title("h₋₁", fontsize=title_font_size)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax3 = fig.add_subplot(1,5,3)
ax3.plot(t, real.(noise_h₀), "k")
title("h₀", fontsize=title_font_size)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax4 = fig.add_subplot(1,5,4)
ax4.plot(t, real.(noise_h₁), "k")
title("h₁", fontsize=title_font_size)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax5 = fig.add_subplot(1,5,5)
ax5.plot(t, real.(noise_h₂), "k")
title("h₂", fontsize=title_font_size)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
subplots_adjust(hspace=0.3, wspace=0.4, top=0.81, left=0.08, right=0.93, bottom=0.06)
savefig("figures/ch2_oscillator_noise.pdf", dpi=300)