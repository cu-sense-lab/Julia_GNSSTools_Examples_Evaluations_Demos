using GNSSTools
using PyPlot
pygui(true)
PyPlot.matplotlib[:rc]("text", usetex=false) # allow tex rendering

σ_n = 1
SNR = 9  # dB
v = Array(0:0.01:10)
pfa = false_alarm_pdf.(v, σ_n)
pd = detection_pdf.(v, σ_n; SNR=SNR)
V_T = 1.9
V_T_idx = findmin(abs.(v .- V_T))[2]
pfa_color = "grey"
pd_color = "lightgrey"

ymax = max(maximum(pfa), maximum(pd))
fig = figure(figsize=(6.5, 1.5))
ax = fig.add_subplot(1,1,1)
ax.plot(v, pfa, pfa_color, label=L"p_{fa}(v)") 
ax.plot(v, pd, pd_color, label=L"p_{d}(v)") 
ax.vlines(x=1.9, ymin=0, ymax=ymax, color="k", linestyle="--")#, label=L"V_T")
fill_between(v[V_T_idx:end], pd[V_T_idx:end], color=pd_color)
fill_between(v[V_T_idx:end], pfa[V_T_idx:end], color=pfa_color)
ylim(0, ymax+0.05)
xlabel("v")
ylabel("pᵥ(v)")
legend()
subplots_adjust(bottom=0.3, left=0.09, right=0.91)
savefig("figures/ch3_probability_detection_false_alarm.svg", dpi=300)