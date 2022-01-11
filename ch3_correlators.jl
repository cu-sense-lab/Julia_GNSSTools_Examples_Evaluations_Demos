using PyPlot
pygui(true)


M = 3
f_s = 11
d = floor(Int, (f_s-1)/2)
N = f_s*M
center_idx = floor(Int, (N+1)/2)
t = Array(range(-M/2, M/2, length=N))
δt = 3
s = zeros(Int, N)
s[center_idx-d:center_idx+d] .= 1
zp = circshift(s, δt)
ze = circshift(zp, -d)
zl = circshift(zp, d)
S = float.(deepcopy(s))
S[center_idx-2*d:center_idx-1] .= Array(0:2*d-1)/(2*d)
S[center_idx+1:center_idx+2*d] .= reverse(S[center_idx-2*d:center_idx-1])

fig = figure(figsize=(7.5, 3))
ax1 = subplot2grid((4,2), (0,0), rowspan=1, colspan=1) 
ax1.step(t, ze, "b-")
ylim([0, 1.2])
ax1.axis("off")
ax2 = subplot2grid((4,2), (1,0), rowspan=1, colspan=1) 
ax2.step(t, zp, color="grey", linestyle="-")
ylim([0, 1.2])
ax2.axis("off")
ax3 = subplot2grid((4,2), (2,0), rowspan=1, colspan=1)
ax3.step(t, zl, "r-")
ylim([0, 1.2])
ax3.axis("off")
ax4 = subplot2grid((4,2), (3,0), rowspan=1, colspan=1)
ax4.step(t, s, "k-")
ax4.arrow(t[1], 0, t[end]-t[1], 0., fc="k", ec="k", lw = 1., 
             head_width=8/80*(t[end]-t[1]), head_length=1.75/20*1.2, 
             overhang=0.3, length_includes_head=true, clip_on=false)
ax4.arrow(0, 0, 0, 11, fc="k", ec="k", lw = 1., 
             head_width=2/30*1.2, head_length=9/60*(t[end]-t[1]), 
             overhang=0.3, length_includes_head=true, clip_on=false) 
ylim([0, 1.2])
ax4.axis("off")

ax5 = subplot2grid((4,2), (0,1), rowspan=4, colspan=1)
ax5.plot(t, S, "k-", label=L"s[n]")
ax5.vlines(x=t[center_idx+δt-d], ymin=0, ymax=S[center_idx+δt-d], 
           color="grey", linestyle=":")
ax5.vlines(x=t[center_idx+δt], ymin=0, ymax=S[center_idx+δt], 
           color="grey", linestyle=":")
ax5.vlines(x=t[center_idx+δt+d], ymin=0, ymax=S[center_idx+δt+d], 
           color="grey", linestyle=":")
ax5.scatter(t[center_idx+δt-d], S[center_idx+δt-d], c="b", label=L"Z_E")
ax5.scatter(t[center_idx+δt], S[center_idx+δt], c="grey", label=L"Z_P")
ax5.scatter(t[center_idx+δt+d], S[center_idx+δt+d], c="r", label=L"Z_L")
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax5.spines["top"].set_visible(false)
ax5.spines["right"].set_visible(false)
ax5.spines["bottom"].set_visible(false)
ax5.spines["left"].set_visible(false)
# ax5.spines["left"].set_position("center")
ax5.arrow(t[1], 0, t[end]-t[1], 0., fc="k", ec="k", lw = 1., 
             head_width=1/80*(t[end]-t[1]), head_length=1.5/20*1.2, 
             overhang=0.3, length_includes_head=true, clip_on=false)
ax5.arrow(0, 0, 0, 1.285, fc="k", ec="k", lw = 1., 
             head_width=2/30*1.2, head_length=1/60*(t[end]-t[1]), 
             overhang=0.3, length_includes_head=true, clip_on=false) 
# xlim(t[1], t[end])
ylim([0, 1.2])
legend()

subplots_adjust(wspace=0.075, hspace=1.5, left=0.01, right=0.99)
savefig("figures/ch3_correlators_for_figure.svg", dpi=300)