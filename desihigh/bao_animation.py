"""
Baryon Acoustic Oscillation (BAO) animation.

Run:
    python bao_animation.py
    python bao_animation.py --save bao.mp4
    python bao_animation.py --save bao.gif
"""

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.widgets import Button, Slider

# ---------- constants ----------
T_RECOMB = 0.55
Z_START, Z_RECOMB = 6000, 1100
DM_FRACTION = 0.83

W, H = 800, 576
CX, CY = W / 2, H / 2
MAX_R = min(W, H) * 0.42

PHASES = [
    (0.00, 0.05, "Radiation era — DM (83%) and baryons (17%) concentrated at centre"),
    (0.05, 0.50, "Sound wave sweeps outward — only baryons inside the wave front are carried"),
    (0.50, 0.58, "Recombination (z ≈ 1100) — wave freezes; shell of baryons + DM centre"),
    (0.58, 1.01, "Matter era — DM dominates central clump; baryon shell at ~150 Mpc"),
]

DM_RGB = np.array([90/255, 150/255, 255/255])
BAR_RGB = np.array([250/255, 165/255, 55/255])


def t_to_z(t):
    if t <= T_RECOMB:
        f = t / T_RECOMB
        return int(round(np.exp(np.log(Z_START) + f * (np.log(Z_RECOMB) - np.log(Z_START)))))
    return int(round(Z_RECOMB * (1 - (t - T_RECOMB) / (1 - T_RECOMB))))


def format_z(z):
    if z >= 1000: return f"z = {z:,}"
    if z >= 10:   return f"z = {z}"
    return f"z = {z:.1f}"


def sound_radius(t):
    return min(t * MAX_R * 0.85, MAX_R * 0.72)


def baryon_r(r0, shell_flag, t):
    rs_now = sound_radius(t)
    rs_end = sound_radius(T_RECOMB)
    out = np.array(r0, copy=True)

    non_shell = shell_flag == 0
    infall = 1.0 - min((t - T_RECOMB) * 0.05, 0.07) if t > T_RECOMB else 1.0
    out[non_shell] = r0[non_shell] * infall

    shell = ~non_shell
    if t < T_RECOMB:
        reached = shell & (r0 < rs_now)
        lag = r0[reached] * 0.15
        out[reached] = np.minimum(rs_now - lag, rs_now * 0.98)
    else:
        reached = shell & (r0 < rs_end)
        settle = min((t - T_RECOMB) * 0.04, 0.03)
        lag = r0[reached] * 0.12
        out[reached] = (rs_end - lag) * (1.0 - settle)
    return out


def phase_text(t):
    for a, b, msg in PHASES:
        if a <= t < b:
            return msg
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--particles", type=int, default=3000)
    ap.add_argument("--speed", type=float, default=0.75)  # 4x slower than prior default
    ap.add_argument("--frames", type=int, default=500)
    ap.add_argument("--save", type=str, default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n = args.particles
    angle = rng.uniform(0, 2 * np.pi, n)
    jitter = rng.uniform(0, 2 * np.pi, n)
    is_dm = rng.random(n) < DM_FRACTION
    r0 = MAX_R * rng.random(n) ** np.where(is_dm, 0.55, 0.58)
    shell = np.where(is_dm, 0.0, (rng.random(n) < 0.35).astype(float))
    dm_idx = np.where(is_dm)[0]
    bar_idx = np.where(~is_dm)[0]

    # stars
    si = np.arange(180)
    sx = (np.sin(si * 127.1 + 3) * 0.5 + 0.5) * W
    sy = (np.sin(si * 311.7 + 7) * 0.5 + 0.5) * H

    fig = plt.figure(figsize=(8, 6.3), facecolor="#080c18")
    ax = fig.add_axes([0, 0.07, 1, 0.82])
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_aspect("equal"); ax.set_facecolor("#080c18"); ax.axis("off")

    star_colors = np.ones((180, 4)); star_colors[:, 3] = 0.15
    stars = ax.scatter(sx, sy, s=3, c=star_colors, zorder=1)

    theta = np.linspace(0, 2 * np.pi, 181)
    wave_ring, = ax.plot([], [], color=(30/255, 220/255, 180/255), lw=1.8,
                         dashes=(6, 4), zorder=2)
    shell_ring, = ax.plot([], [], color=(220/255, 130/255, 60/255), lw=1, zorder=3)

    dm_colors = np.tile(np.append(DM_RGB, 0.5), (len(dm_idx), 1))
    bar_colors = np.tile(np.append(BAR_RGB, 0.8), (len(bar_idx), 1))
    dm_scatter = ax.scatter(np.zeros(len(dm_idx)), np.zeros(len(dm_idx)),
                            s=5, c=dm_colors, edgecolors="none", zorder=4)
    bar_scatter = ax.scatter(np.zeros(len(bar_idx)), np.zeros(len(bar_idx)),
                             s=8, c=bar_colors, edgecolors="none", zorder=5)

    glow = ax.scatter([CX], [CY], s=900, c=[[1, 0.82, 0.43, 0.25]],
                      edgecolors="none", zorder=6)
    ax.scatter([CX], [CY], s=50, c="#fff5e0", edgecolors="none", zorder=7)

    z_text = ax.text(0.02, 0.97, "", transform=ax.transAxes,
                     color="white", fontsize=13, fontweight="bold", va="top")
    phase = fig.text(0.5, 0.03, "", ha="center", color="#9aa4b4", fontsize=10)

    legend = [
        mpatches.Patch(color="#5a96ff", label="Dark matter ~83%"),
        mpatches.Patch(color="#faa537", label="Baryons ~17%"),
        mpatches.Patch(color=(30/255, 220/255, 180/255), label="Sound wave front"),
    ]
    ax.legend(handles=legend, loc="upper right", frameon=False,
              labelcolor="white", fontsize=9)

    state = {"t": 0.0, "playing": True, "speed": args.speed}

    def update(frame):
        if state["playing"]:
            state["t"] = (state["t"] + 0.0008 * state["speed"] * 8) % 1.0
        t = state["t"]
        rs = sound_radius(t)

        if t < T_RECOMB:
            wave_ring.set_data(CX + rs * np.cos(theta), CY + rs * np.sin(theta))
            wave_ring.set_alpha(0.65); wave_ring.set_linewidth(1.8)
            shell_ring.set_data([], [])
        else:
            rs_f = sound_radius(T_RECOMB)
            fade = max(0.0, 1 - (t - T_RECOMB) * 2)
            wave_ring.set_data(CX + rs_f * np.cos(theta), CY + rs_f * np.sin(theta))
            wave_ring.set_alpha(fade * 0.4); wave_ring.set_linewidth(1.2)
            sh_r = rs_f * (1 + (t - T_RECOMB) * 0.02)
            shell_ring.set_data(CX + sh_r * np.cos(theta), CY + sh_r * np.sin(theta))
            shell_ring.set_alpha(min((t - T_RECOMB) * 1.5, 0.6))
            shell_ring.set_linewidth(min((t - T_RECOMB) * 20 + 1, 8))

        infall_dm = 1.0 - min((t - T_RECOMB) * 0.06, 0.08) if t > T_RECOMB else 1.0
        r_dm = r0[dm_idx] * infall_dm
        jx = np.sin(jitter[dm_idx] * 3.1 + t * 1.4) * 1.1
        jy = np.cos(jitter[dm_idx] * 2.3 + t * 1.8) * 1.1
        dm_scatter.set_offsets(np.column_stack([
            CX + np.cos(angle[dm_idx]) * r_dm + jx,
            CY + np.sin(angle[dm_idx]) * r_dm + jy,
        ]))
        a_dm = np.clip(0.70 - (r_dm / MAX_R) * 0.80, 0.06, 1)
        dm_colors[:, 3] = a_dm
        dm_scatter.set_facecolors(dm_colors)

        r_b = np.minimum(baryon_r(r0[bar_idx], shell[bar_idx], t), MAX_R * 0.97)
        jx = np.sin(jitter[bar_idx] * 3.1 + t * 1.6) * 1.2
        jy = np.cos(jitter[bar_idx] * 2.3 + t * 2.0) * 1.2
        bar_scatter.set_offsets(np.column_stack([
            CX + np.cos(angle[bar_idx]) * r_b + jx,
            CY + np.sin(angle[bar_idx]) * r_b + jy,
        ]))
        a_b = np.clip(0.95 - (r_b / MAX_R) * 0.70, 0.14, 1)
        bar_colors[:, 3] = a_b
        bar_scatter.set_facecolors(bar_colors)

        glow_r = 22 + (t - T_RECOMB) * 6 if t > T_RECOMB else 18
        glow.set_sizes([glow_r * glow_r * 2])

        z_text.set_text(format_z(t_to_z(t)))
        phase.set_text(phase_text(t))
        fig.canvas.draw_idle()
        return stars, wave_ring, shell_ring, dm_scatter, bar_scatter, glow, z_text, phase

    anim = FuncAnimation(fig, update, frames=args.frames,
                         interval=33, blit=False, repeat=True, cache_frame_data=False)

    # Speed slider (top-left)
    spd_ax = fig.add_axes([0.08, 0.94, 0.18, 0.025])
    spd_ax.set_facecolor("#1a2238")
    spd = Slider(spd_ax, "Speed", 0.1, 5.0, valinit=args.speed,
                 valstep=0.05, color="#4a79d4")
    spd.label.set_color("white"); spd.label.set_fontsize(9)
    spd.valtext.set_color("white"); spd.valtext.set_fontsize(9)
    spd.on_changed(lambda v: state.__setitem__("speed", float(v)))

    # Play/Pause button
    btn_ax = fig.add_axes([0.30, 0.925, 0.12, 0.05])
    btn = Button(btn_ax, "Pause", color="#1a2238", hovercolor="#2a3454")
    btn.label.set_color("white")

    def toggle(event):
        state["playing"] = not state["playing"]
        btn.label.set_text("Play" if not state["playing"] else "Pause")
        fig.canvas.draw_idle()

    btn.on_clicked(toggle)

    # Save PNG (current frame)
    btn_png_ax = fig.add_axes([0.44, 0.925, 0.12, 0.05])
    btn_png = Button(btn_png_ax, "Save PNG", color="#1a2238", hovercolor="#2a3454")
    btn_png.label.set_color("white")

    def save_png(event):
        from datetime import datetime
        fname = f"bao_frame_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        fig.savefig(fname, facecolor=fig.get_facecolor(), dpi=150)
        print(f"saved → {fname}")

    btn_png.on_clicked(save_png)

    # Save MP4 (full animation)
    btn_mp4_ax = fig.add_axes([0.58, 0.925, 0.14, 0.05])
    btn_mp4 = Button(btn_mp4_ax, "Save MP4", color="#1a2238", hovercolor="#2a3454")
    btn_mp4.label.set_color("white")

    def save_mp4(event):
        from datetime import datetime
        was_playing = state["playing"]
        state["playing"] = False
        saved_t = state["t"]
        fname = f"bao_{datetime.now().strftime('%Y%m%d_%H%M%S')}.mp4"
        print(f"writing {fname} … ({args.frames} frames)")
        try:
            state["t"] = 0.0
            FFMpegWriter(fps=30, bitrate=4000)
            anim.save(fname, writer=FFMpegWriter(fps=30, bitrate=4000))
            print(f"saved → {fname}")
        except Exception as e:
            # fall back to GIF if ffmpeg missing
            gif = fname.replace(".mp4", ".gif")
            print(f"ffmpeg failed ({e}); writing GIF instead → {gif}")
            state["t"] = 0.0
            anim.save(gif, writer=PillowWriter(fps=30))
            print(f"saved → {gif}")
        state["t"] = saved_t
        state["playing"] = was_playing

    btn_mp4.on_clicked(save_mp4)

    fig._btns = (btn, btn_png, btn_mp4, spd)  # keep refs

    if args.save:
        writer = (PillowWriter(fps=30) if args.save.endswith(".gif")
                  else FFMpegWriter(fps=30, bitrate=4000))
        anim.save(args.save, writer=writer)
        print(f"saved → {args.save}")
    else:
        # keep anim reference on figure to prevent GC on some backends
        fig._anim = anim
        plt.show()


if __name__ == "__main__":
    main()
