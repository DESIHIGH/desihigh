"""
desi_sky_subtraction.py — animated DESI-style sky subtraction.

Builds up a galaxy-fiber spectrum and 10 sky-fiber spectra photon by photon
(Poisson increments) over a 600s exposure, recomputing the sky model
(median of sky fibers) and the sky-subtracted spectrum every frame.

Four panels (top to bottom):
  1. Raw galaxy fiber — dominated by OH sky lines, [OII] doublet buried in it.
  2. Sky model (median of 10 sky fibers) ± 1σ band.
  3. Sky-subtracted — [OII] emerges as exposure grows; true signal overlaid.
  4. Residual (sky-subtracted − true signal) vs ±1σ noise envelope.

Controls:  Speed slider · Shutter / Pause / Reset · Save PNG.

Run:
    python desi_sky_subtraction.py
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation

# ---------- constants ----------
LAM_LO, LAM_HI = 738.0, 756.0
N_PIX = 720
lam = np.linspace(LAM_LO, LAM_HI, N_PIX)

T_EXP = 600.0
READ_NOISE = 3.0
LSF_SIGMA = 0.11
N_SKY_FIBERS = 10

OII_LINES = [
    {"lam": 746.2, "rate": 8.0, "label": "[OII] 746.2"},
    {"lam": 746.8, "rate": 8.0, "label": "[OII] 746.8"},
]
OH_LINES = [
    {"lam": 739.2, "rate": 55.0},
    {"lam": 747.2, "rate": 70.0},
    {"lam": 752.3, "rate": 48.0},
    {"lam": 754.5, "rate": 38.0},
]
SKY_CONTINUUM = 8.0


def gauss(lam_grid, lam0, sigma):
    return np.exp(-0.5 * ((lam_grid - lam0) / sigma) ** 2)


def make_sky_rate():
    r = np.full(N_PIX, SKY_CONTINUUM)
    for line in OH_LINES:
        r += line["rate"] * gauss(lam, line["lam"], LSF_SIGMA)
    return r


def make_obj_rate():
    r = np.zeros(N_PIX)
    for line in OII_LINES:
        r += line["rate"] * gauss(lam, line["lam"], LSF_SIGMA)
    return r


SKY_RATE = make_sky_rate()
OBJ_RATE = make_obj_rate()
TRUE_OBJ = OBJ_RATE * T_EXP  # for reference line only (noiseless, at full exposure)


# ---------- simulation state ----------
class Sim:
    def __init__(self, seed=42):
        self.rng = np.random.default_rng(seed)
        self.speed = 20.0   # sim seconds per real second
        self.running = False
        self.reset()

    def reset(self):
        self.rng = np.random.default_rng(self.rng.integers(1 << 31))
        self.t = 0.0
        self.galaxy = self.rng.normal(0, READ_NOISE, N_PIX)
        self.sky_fibers = self.rng.normal(0, READ_NOISE, (N_SKY_FIBERS, N_PIX))
        self.running = False

    def step(self, dt):
        if self.t >= T_EXP or dt <= 0:
            return
        dt = min(dt, T_EXP - self.t)
        self.galaxy += self.rng.poisson((SKY_RATE + OBJ_RATE) * dt).astype(float)
        self.sky_fibers += self.rng.poisson(
            np.broadcast_to(SKY_RATE * dt, (N_SKY_FIBERS, N_PIX))).astype(float)
        self.t += dt
        if self.t >= T_EXP:
            self.running = False

    # derived products
    def sky_model(self):
        return np.median(self.sky_fibers, axis=0)

    def sky_model_err(self):
        # median of N Poisson draws: σ ≈ 1.253 σ_mean / √N
        if self.t <= 0:
            return np.zeros(N_PIX)
        return (1.253 * np.sqrt(SKY_RATE * self.t + READ_NOISE ** 2)
                / np.sqrt(N_SKY_FIBERS))

    def noise_theory(self):
        if self.t <= 0:
            return np.full(N_PIX, READ_NOISE)
        return np.sqrt(OBJ_RATE * self.t + SKY_RATE * self.t
                       + READ_NOISE ** 2 + self.sky_model_err() ** 2)

    def snr(self, lam0, nsig=3):
        mask = np.abs(lam - lam0) < nsig * LSF_SIGMA
        sub = self.galaxy - self.sky_model()
        sig = sub[mask].sum()
        n = np.sqrt((self.noise_theory()[mask] ** 2).sum())
        return sig / n if n > 0 else 0.0


# ---------- UI ----------
def main():
    sim = Sim()

    fig = plt.figure(figsize=(13, 9), facecolor="#07070f")
    fig.canvas.manager.set_window_title("DESI sky subtraction — animated")

    gs = gridspec.GridSpec(4, 1, hspace=0.12,
                           top=0.87, bottom=0.09, left=0.08, right=0.94)

    PANEL_BG = "#0d0d1a"
    SPINE = "#333355"
    TICK = "#888"
    TEXT = "#ccccdd"
    GRID = "#1a1a2e"

    def style(ax, ylabel):
        ax.set_facecolor(PANEL_BG)
        for s in ax.spines.values(): s.set_color(SPINE)
        ax.tick_params(colors=TICK, labelcolor=TEXT, labelsize=8)
        ax.set_ylabel(ylabel, color=TEXT, fontsize=9)
        ax.yaxis.set_label_coords(-0.06, 0.5)
        ax.grid(True, color=GRID, linewidth=0.4)
        ax.set_xlim(LAM_LO, LAM_HI)

    # Panel 1: raw galaxy fiber
    ax1 = fig.add_subplot(gs[0])
    l1_raw, = ax1.plot(lam, sim.galaxy, color="#4488ff", lw=0.7, alpha=0.85,
                       label="raw galaxy fiber")
    l1_sky, = ax1.plot(lam, np.zeros_like(lam), color="#ff8833", lw=1.0,
                       ls="--", alpha=0.7, label="true sky × t")
    for line in OII_LINES:
        ax1.axvline(line["lam"], color="#ffee44", lw=0.8, ls=":", alpha=0.5)
    style(ax1, "counts (e⁻)")
    ax1.legend(loc="upper left", fontsize=7, facecolor=PANEL_BG,
               labelcolor=TEXT, framealpha=0.7)
    title = ax1.set_title("", color=TEXT, fontsize=10, pad=6)
    ax1.tick_params(labelbottom=False)

    # Panel 2: sky fibers + model
    ax2 = fig.add_subplot(gs[1])
    sky_lines = []
    for i in range(N_SKY_FIBERS):
        l, = ax2.plot(lam, sim.sky_fibers[i], color="#336688",
                      lw=0.4, alpha=0.35,
                      label="sky fiber" if i == 0 else None)
        sky_lines.append(l)
    l2_model, = ax2.plot(lam, sim.sky_model(), color="#55aaff", lw=1.2,
                         label=f"sky model (median of {N_SKY_FIBERS})")
    l2_err = ax2.fill_between(lam, sim.sky_model() - 1, sim.sky_model() + 1,
                              color="#55aaff", alpha=0.15, label="model ±1σ")
    style(ax2, "counts (e⁻)")
    ax2.legend(loc="upper left", fontsize=7, facecolor=PANEL_BG,
               labelcolor=TEXT, framealpha=0.7)
    ax2.tick_params(labelbottom=False)

    # Panel 3: sky-subtracted
    ax3 = fig.add_subplot(gs[2])
    sub0 = sim.galaxy - sim.sky_model()
    l3_fill = ax3.fill_between(lam, sub0, alpha=0.25, color="#ffee44")
    l3_sub, = ax3.plot(lam, sub0, color="#ffee44", lw=0.8, label="sky-subtracted")
    l3_true, = ax3.plot(lam, TRUE_OBJ, color="#ff4444", lw=1.2, ls="--",
                        alpha=0.8, label="true [OII] × T_EXP (noiseless)")
    ax3.axhline(0, color=SPINE, lw=0.6)
    for line in OII_LINES:
        ax3.axvline(line["lam"], color="#ffee44", lw=0.8, ls=":", alpha=0.6)
    style(ax3, "counts (e⁻)")
    ax3.legend(loc="upper left", fontsize=7, facecolor=PANEL_BG,
               labelcolor=TEXT, framealpha=0.7)
    snr_txt = ax3.text(0.98, 0.90, "", transform=ax3.transAxes,
                       ha="right", va="top", color="#ffee44", fontsize=8,
                       bbox=dict(facecolor=PANEL_BG, edgecolor=SPINE, alpha=0.8))
    ax3.tick_params(labelbottom=False)

    # Panel 4: residual
    ax4 = fig.add_subplot(gs[3])
    zeros = np.zeros_like(lam)
    l4_env = ax4.fill_between(lam, zeros, zeros, color="#334466", alpha=0.4,
                              label="±1σ noise envelope")
    l4_res, = ax4.plot(lam, zeros, color="#aaaacc", lw=0.6, alpha=0.8,
                       label="residual = subtracted − true")
    ax4.axhline(0, color=SPINE, lw=0.6)
    for line in OH_LINES:
        ax4.axvline(line["lam"], color="#ff8833", lw=0.5, ls=":", alpha=0.4)
    style(ax4, "residual (e⁻)")
    ax4.set_xlabel("Observed wavelength (nm)", color=TEXT, fontsize=9)
    ax4.legend(loc="upper left", fontsize=7, facecolor=PANEL_BG,
               labelcolor=TEXT, framealpha=0.7)

    # --- controls (top row) ---
    sax = fig.add_axes([0.10, 0.93, 0.25, 0.025], facecolor="#1a2238")
    spd = Slider(sax, "Speed (sim s / real s)", 1.0, 100.0, valinit=sim.speed,
                 valstep=1, color="#4a79d4")
    spd.label.set_color("white"); spd.label.set_fontsize(8)
    spd.valtext.set_color("white"); spd.valtext.set_fontsize(8)
    spd.on_changed(lambda v: setattr(sim, "speed", float(v)))

    def mkbtn(x, y, w, h, label):
        a = fig.add_axes([x, y, w, h], facecolor="#2a3454")
        b = Button(a, label, color="#2a3454", hovercolor="#3d4a70")
        b.label.set_color("white"); b.label.set_fontsize(9)
        return b

    btn_go    = mkbtn(0.45, 0.93, 0.10, 0.035, "Shutter")
    btn_pause = mkbtn(0.56, 0.93, 0.08, 0.035, "Pause")
    btn_reset = mkbtn(0.65, 0.93, 0.08, 0.035, "Reset")
    btn_save  = mkbtn(0.84, 0.93, 0.10, 0.035, "Save PNG")

    def on_go(_):    sim.running = (sim.t < T_EXP);
    def on_pause(_): sim.running = False
    def on_reset(_):
        sim.reset()
        rerender()

    def on_save(_):
        from datetime import datetime
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           f"desi_sky_subtraction_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png")
        fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
        print(f"saved → {out}")

    btn_go.on_clicked(on_go)
    btn_pause.on_clicked(on_pause)
    btn_reset.on_clicked(on_reset)
    btn_save.on_clicked(on_save)

    # --- rendering ---
    def rerender():
        model = sim.sky_model()
        sub = sim.galaxy - model
        merr = sim.sky_model_err()
        nth = sim.noise_theory()
        res = sub - TRUE_OBJ * (sim.t / T_EXP)  # compare to true signal *so far*

        l1_raw.set_ydata(sim.galaxy)
        l1_sky.set_ydata(SKY_RATE * sim.t)
        ymax = max(sim.galaxy.max() * 1.1, 10)
        ax1.set_ylim(min(sim.galaxy.min() - 20, -50), ymax)

        for i, ln in enumerate(sky_lines):
            ln.set_ydata(sim.sky_fibers[i])
        l2_model.set_ydata(model)
        # rebuild fill
        nonlocal l2_err
        l2_err.remove()
        l2_err = ax2.fill_between(lam, model - merr, model + merr,
                                  color="#55aaff", alpha=0.15)
        ax2.set_ylim(min(0, sim.sky_fibers.min() - 20),
                     max(sim.sky_fibers.max() * 1.1, 10))

        l3_sub.set_ydata(sub)
        nonlocal l3_fill
        l3_fill.remove()
        l3_fill = ax3.fill_between(lam, sub, 0, alpha=0.25, color="#ffee44")
        l3_true.set_ydata(TRUE_OBJ * (sim.t / T_EXP))
        ymax3 = max(sub.max() * 1.2, TRUE_OBJ.max() * (sim.t / T_EXP) * 1.2, 50)
        ymin3 = min(sub.min() * 1.1, -50)
        ax3.set_ylim(ymin3, ymax3)

        l4_res.set_ydata(res)
        nonlocal l4_env
        l4_env.remove()
        l4_env = ax4.fill_between(lam, nth, -nth, color="#334466", alpha=0.4)
        ye = max(nth.max() * 1.3, abs(res).max() * 1.1, 30)
        ax4.set_ylim(-ye, ye)

        snr_a = sim.snr(OII_LINES[0]["lam"])
        snr_b = sim.snr(OII_LINES[1]["lam"])
        snr_txt.set_text(f"t = {sim.t:6.1f} s\n"
                         f"SNR 746.2 = {snr_a:5.1f}\n"
                         f"SNR 746.8 = {snr_b:5.1f}")
        title.set_text(
            f"DESI sky subtraction — t = {sim.t:.1f}/{int(T_EXP)} s  "
            f"({N_SKY_FIBERS} sky fibers)  "
            f"{'[EXPOSING]' if sim.running else '[paused]' if sim.t < T_EXP else '[DONE]'}"
        )

    state = {"last": None}

    def tick(_):
        if sim.running:
            now = time.time()
            wall_dt = 0.033 if state["last"] is None else min(now - state["last"], 0.1)
            state["last"] = now
            sim.step(wall_dt * sim.speed)
        else:
            state["last"] = None
        rerender()
        return []

    anim = FuncAnimation(fig, tick, interval=33, blit=False,
                         cache_frame_data=False)

    fig._keep = (spd, btn_go, btn_pause, btn_reset, btn_save, anim)

    rerender()
    plt.show()


if __name__ == "__main__":
    main()
