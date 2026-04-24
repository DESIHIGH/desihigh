"""
CCD / multi-fiber spectrograph mapping simulator.

20 fibers × 120 wavelength rows, three arms (blue/red/NIR), with an expose
button that adds sky + object + photon + read noise, cosmic-ray injection,
and a "clock" operation that shifts columns into a serial register.

Run:
    python ccd_mapping.py

Fiber selector: slider 0..19 below the CCD, or click a fiber column on the image.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.patches import Rectangle

# ---------- constants ----------
NF, NW, FW = 20, 120, 100000.0

ARMS = {
    "blue": {"lo": 360, "hi": 550, "name": "Blue 360–550nm", "skysc": 0.4,
             "sky_lines": [487, 520]},
    "red":  {"lo": 550, "hi": 720, "name": "Red 550–720nm",  "skysc": 1.0,
             "sky_lines": [589, 630, 666, 700]},
    "nir":  {"lo": 720, "hi": 980, "name": "NIR 720–980nm",  "skysc": 2.5,
             "sky_lines": [742, 780, 812, 855, 904, 940]},
}

FIBERS = [
    ("sky",    0.00), ("galaxy", 0.08), ("galaxy", 0.31), ("qso",    0.95),
    ("galaxy", 0.17), ("star",   0.00), ("galaxy", 0.52), ("galaxy", 0.44),
    ("sky",    0.00), ("galaxy", 0.23), ("qso",    1.28), ("galaxy", 0.67),
    ("star",   0.00), ("galaxy", 0.88), ("galaxy", 0.11), ("qso",    0.41),
    ("galaxy", 0.76), ("sky",    0.00), ("galaxy", 0.35), ("galaxy", 0.59),
]

# name, rest-frame nm, equivalent-width weight, RGB hint
ALL_LINES = [
    ("[OII]",  372.7, 40), ("Hδ",    410.2,  6), ("Hγ",    434.0,  9),
    ("Hβ",     486.1, 22), ("[OIII]",495.9, 15), ("[OIII]",500.7, 45),
    ("Hα",     656.3, 65), ("[NII]", 658.4, 22), ("[SII]", 671.6, 14),
    ("Paβ",    821.0,  9), ("CaII",  854.2, 12), ("CaII",  866.2,  8),
    ("[SIII]", 907.0,  7),
]


def lam_rgb(lam):
    """Approximate visible-light wavelength → RGB."""
    if lam < 380:   return (55,  0,  90)
    if lam < 440:   t = (440 - lam) / 60;  return (int(t * 80), 0, 255)
    if lam < 490:   t = (lam - 440) / 50;  return (0, int(t * 210), 255)
    if lam < 510:   t = (510 - lam) / 20;  return (0, 255, int(t * 255))
    if lam < 580:   t = (lam - 510) / 70;  return (int(t * 255), 255, 0)
    if lam < 645:   t = (645 - lam) / 65;  return (255, int(t * 210), 0)
    return (255, 0, 0)


# ---------- model state ----------
class Sim:
    def __init__(self):
        self.arm_key = "red"
        self.lstr = 8; self.cont = 2; self.sky = 3; self.rn = 2
        self.grid = np.zeros((NF, NW), dtype=np.float32)
        self.cr = set()
        self.ck_col = 0
        self.sel_f = 1

    @property
    def arm(self):
        return ARMS[self.arm_key]

    def w_to_lam(self, w):
        a = self.arm
        return a["hi"] - w / NW * (a["hi"] - a["lo"])

    def sky_at(self, lam):
        a = self.arm
        s = 600 * self.sky * a["skysc"]
        for ll in a["sky_lines"]:
            d = lam - ll
            s += 600 * 8 * self.sky * a["skysc"] * np.exp(-d * d / 10.0)
        return s

    def obj_at(self, lam, fib_type, z):
        if fib_type == "sky":
            return 0.0
        a = self.arm
        s = 0.0
        f = (lam - a["lo"]) / (a["hi"] - a["lo"])
        if fib_type == "galaxy":
            s += 700 * self.cont * (0.7 + 0.5 * f)
        elif fib_type == "qso":
            s += 700 * self.cont * 2 * (0.3 + 1.5 * max(0, f) ** 1.2)
        elif fib_type == "star":
            s += 700 * self.cont * 3 * np.exp(-((lam - 580) ** 2) / 11000)
        for _, rest, eqw in ALL_LINES:
            obs = rest * (1 + z)
            d = lam - obs
            if abs(d) < 10:
                s += 3500 * self.lstr * (eqw / 40) * np.exp(-d * d / 0.7)
        return s

    def expose(self):
        rng = np.random.default_rng()
        a = self.arm
        ws = np.arange(NW)
        lams = a["hi"] - ws / NW * (a["hi"] - a["lo"])
        for fi in range(NF):
            ftype, z = FIBERS[fi]
            sig = np.array([self.sky_at(l) + self.obj_at(l, ftype, z) for l in lams])
            photon = rng.normal(0, np.sqrt(np.maximum(sig, 1)))
            read = rng.normal(0, self.rn, size=NW)
            self.grid[fi] = np.clip(sig + photon + read, 0, FW)
        for (fi, w) in self.cr:
            if fi < NF and w < NW:
                self.grid[fi, w] = min(FW, self.grid[fi, w] + FW * 0.7)

    def add_cr(self):
        rng = np.random.default_rng()
        fi0 = rng.integers(NF); w0 = rng.integers(NW)
        length = 2 + rng.integers(5)
        dfi = 1 if rng.random() > 0.5 else 0
        dw  = 1 if rng.random() > 0.4 else 0
        for i in range(length):
            fi = min(NF - 1, fi0 + i * dfi)
            w  = min(NW - 1, w0  + i * dw)
            self.cr.add((fi, w))
            self.grid[fi, w] = min(FW, self.grid[fi, w] + FW * 0.7)

    def reset(self):
        self.grid[:] = 0; self.cr.clear(); self.ck_col = 0

    def clock_one(self):
        if self.ck_col >= NF:
            return False
        self.grid[:-1] = self.grid[1:]
        self.grid[-1] = 0
        self.ck_col += 1
        return True

    def lines_in_arm(self, fi):
        ftype, z = FIBERS[fi]
        if ftype in ("sky", "star"):
            return 0
        a = self.arm
        n = 0
        for _, rest, _ in ALL_LINES:
            obs = rest * (1 + z)
            if a["lo"] <= obs <= a["hi"]:
                n += 1
        return n


# ---------- rendering ----------
def ccd_image(sim: Sim):
    a = sim.arm
    dmax = max(sim.grid.max(), 1)
    ws = np.arange(NW)
    lams = a["hi"] - ws / NW * (a["hi"] - a["lo"])
    wl_cols = np.array([lam_rgb(l) for l in lams], dtype=np.float32) / 255.0

    img = np.zeros((NW, NF, 3), dtype=np.float32)
    f = np.sqrt(np.clip(sim.grid.T / dmax, 0, 1))
    dim = f < 0.025
    col = wl_cols[:, None, :] * (f[..., None] * 1.5)
    highlight = np.clip((f - 0.78) / 0.22, 0, 1)[..., None]
    col = col + highlight * (1 - col)
    col = np.clip(col, 0, 1)
    base = np.array([4, 4, 14], dtype=np.float32) / 255.0
    img[:] = np.where(dim[..., None], base, col)

    for (fi, w) in sim.cr:
        if fi < NF and w < NW:
            img[w, fi] = np.array([1.0, 0.51, 0.08])
    return img


def spectrum_points(sim: Sim):
    a = sim.arm
    ws = np.arange(NW)
    lams = a["hi"] - ws / NW * (a["hi"] - a["lo"])
    row = sim.grid[sim.sel_f]
    order = np.argsort(lams)
    return lams[order], row[order]


# ---------- UI ----------
def main():
    sim = Sim()
    sim.expose()

    fig = plt.figure(figsize=(13, 6.8), facecolor="#0b0d14")
    fig.canvas.manager.set_window_title("CCD multi-fiber spectrograph")

    # ---- CCD image axis
    ax_ccd = fig.add_axes([0.09, 0.32, 0.48, 0.60])
    ax_ccd.set_facecolor("#04040e")
    im = ax_ccd.imshow(ccd_image(sim), aspect="auto", origin="upper",
                       extent=[-0.5, NF - 0.5, NW, 0], interpolation="nearest")
    ax_ccd.set_xlabel("fiber", color="#bbb", fontsize=9)
    ax_ccd.set_xticks([0, 5, 10, 15, 19])
    ax_ccd.tick_params(colors="#999", labelsize=8)
    arm_title = ax_ccd.set_title(sim.arm["name"], color="white", fontsize=11)
    # Selected-fiber highlight — use Rectangle directly (compatible with
    # matplotlib 3.9+, where axvspan returns a Rectangle whose set_xy
    # signature differs from the older Polygon's).
    sel_rect = Rectangle((sim.sel_f - 0.5, 0), 1.0, NW,
                         ec=(1, 0.84, 0.2), fc="none", lw=1.6)
    ax_ccd.add_patch(sel_rect)

    # wavelength colour bar (left of CCD)
    ax_wb = fig.add_axes([0.065, 0.32, 0.015, 0.60])
    a = sim.arm
    wbar = np.array([lam_rgb(a["hi"] - y / 300 * (a["hi"] - a["lo"]))
                     for y in range(300)], dtype=np.float32) / 255.0
    wb_im = ax_wb.imshow(wbar[:, None, :], aspect="auto", origin="upper")
    ax_wb.set_xticks([]); ax_wb.tick_params(colors="#999", labelsize=7)
    ax_wb.set_yticks([0, 299]); ax_wb.set_yticklabels([a["hi"], a["lo"]])

    # ---- 1D spectrum axis
    ax_sp = fig.add_axes([0.065, 0.08, 0.50, 0.18])
    ax_sp.set_facecolor("#04040e")
    ax_sp.tick_params(colors="#999", labelsize=8)
    for s in ax_sp.spines.values(): s.set_color("#333")
    lams, row = spectrum_points(sim)
    (sp_line,) = ax_sp.plot(lams, row, color="white", lw=1)
    line_markers = []
    sp_title = ax_sp.set_title("", color="white", fontsize=9, loc="left")

    # ---- right-side controls
    panel_bg = "#1a1f2e"
    def mkax(x, y, w, h, facecolor=panel_bg):
        a = fig.add_axes([x, y, w, h])
        a.set_facecolor(facecolor); return a

    # arm radio
    ax_arm = mkax(0.62, 0.78, 0.11, 0.14)
    arm_radio = RadioButtons(ax_arm, ("blue", "red", "nir"), active=1,
                             activecolor="#ffaa44")
    for lbl in arm_radio.labels: lbl.set_color("white"); lbl.set_fontsize(9)
    ax_arm.set_title("Arm", color="#bbb", fontsize=9)

    # sliders
    slider_axes = {}
    sliders = {}
    slider_defs = [
        ("lstr", "Line str", 1, 20, 8),
        ("cont", "Continuum", 0, 8, 2),
        ("sky",  "Sky level", 0, 10, 3),
        ("rn",   "Read noise", 0, 100, 2),
    ]
    for i, (key, label, lo, hi, val) in enumerate(slider_defs):
        sax = mkax(0.76, 0.87 - i * 0.05, 0.21, 0.025, "#0b0d14")
        slider_axes[key] = sax
        sl = Slider(sax, label, lo, hi, valinit=val, valstep=1,
                    color="#4a79d4")
        sl.label.set_color("#bbb"); sl.label.set_fontsize(8)
        sl.valtext.set_color("white"); sl.valtext.set_fontsize(8)
        sliders[key] = sl

    # buttons
    def mkbtn(x, y, w, h, label):
        a = mkax(x, y, w, h, "#2a3454")
        b = Button(a, label, color="#2a3454", hovercolor="#3d4a70")
        b.label.set_color("white"); b.label.set_fontsize(9)
        return b

    btn_expose = mkbtn(0.62, 0.66, 0.08, 0.05, "Expose")
    btn_cr     = mkbtn(0.71, 0.66, 0.07, 0.05, "+CR")
    btn_reset  = mkbtn(0.79, 0.66, 0.09, 0.05, "Reset")
    btn_ck1    = mkbtn(0.62, 0.60, 0.12, 0.045, "Clock 1")
    btn_ckall  = mkbtn(0.75, 0.60, 0.13, 0.045, "Clock all")
    btn_png    = mkbtn(0.62, 0.54, 0.17, 0.045, "Save PNG")
    btn_data   = mkbtn(0.80, 0.54, 0.18, 0.045, "Save data")

    # summary stats (compact)
    stats_ax = mkax(0.62, 0.42, 0.36, 0.11, "#0f121b")
    stats_ax.set_xticks([]); stats_ax.set_yticks([])
    stats_txt = stats_ax.text(0.04, 0.95, "", color="white",
                              fontsize=9, va="top", family="monospace")

    # info bar
    info_ax = mkax(0.62, 0.36, 0.36, 0.05, "#0f121b")
    info_ax.set_xticks([]); info_ax.set_yticks([])
    info_txt = info_ax.text(0.03, 0.5, "Red arm exposed. Use slider to pick fiber.",
                            color="#bbb", fontsize=8, va="center")

    # fiber table, two columns (0-9 left, 10-19 right) — all 20 visible
    ftable_ax = mkax(0.62, 0.03, 0.36, 0.32, "#0f121b")
    ftable_ax.set_xticks([]); ftable_ax.set_yticks([])
    ftable_ax.text(0.5, 0.97, "Fiber table  (# type  z   lines)",
                   color="#bbb", fontsize=8, ha="center", va="top",
                   transform=ftable_ax.transAxes)
    ftable_left  = ftable_ax.text(0.03, 0.88, "", color="white",
                                  fontsize=8, va="top", family="monospace",
                                  transform=ftable_ax.transAxes)
    ftable_right = ftable_ax.text(0.53, 0.88, "", color="white",
                                  fontsize=8, va="top", family="monospace",
                                  transform=ftable_ax.transAxes)

    # clock-all timer
    timer = [None]

    # ---------- redraw helpers ----------
    def update_arm_bar():
        a = sim.arm
        wbar = np.array([lam_rgb(a["hi"] - y / 300 * (a["hi"] - a["lo"]))
                         for y in range(300)], dtype=np.float32) / 255.0
        wb_im.set_data(wbar[:, None, :])
        ax_wb.set_yticklabels([a["hi"], a["lo"]])
        arm_title.set_text(a["name"])

    def redraw_ccd():
        im.set_data(ccd_image(sim))
        sel_rect.set_bounds(sim.sel_f - 0.5, 0, 1.0, NW)

    def redraw_spectrum():
        nonlocal line_markers
        for m in line_markers:
            m.remove()
        line_markers = []

        a = sim.arm
        lams, row = spectrum_points(sim)
        sp_line.set_data(lams, row)
        ax_sp.set_xlim(a["lo"], a["hi"])
        ymax = max(row.max(), 1)
        ax_sp.set_ylim(0, ymax * 1.15)

        ftype, z = FIBERS[sim.sel_f]
        for ll in a["sky_lines"]:
            m = ax_sp.axvline(ll, color="#ffff99", alpha=0.25, lw=0.8)
            line_markers.append(m)
        if ftype not in ("sky", "star"):
            for name, rest, _ in ALL_LINES:
                obs = rest * (1 + z)
                if a["lo"] <= obs <= a["hi"]:
                    r, g, b = lam_rgb(obs)
                    c = (r/255, g/255, b/255)
                    m = ax_sp.axvline(obs, color=c, alpha=0.7, lw=1, ls="--")
                    line_markers.append(m)
                    t = ax_sp.text(obs, ymax * 1.08, name, color=c,
                                   fontsize=7, ha="center")
                    line_markers.append(t)
        sp_title.set_text(f"fiber {sim.sel_f} · {ftype}"
                          + (f"  z = {z:.2f}" if z > 0 else ""))

    def redraw_stats():
        ftype, z = FIBERS[sim.sel_f]
        row = sim.grid[sim.sel_f]
        peak = max(row.max(), 1)
        a = sim.arm
        skyE = 600 * sim.sky * a["skysc"]
        snr = peak / np.sqrt(peak + skyE + sim.rn ** 2)
        lines = sim.lines_in_arm(sim.sel_f)
        stats_txt.set_text(
            f"fiber {sim.sel_f}  {ftype}  z={z:.2f}\n"
            f"SNR {snr:.1f}   peak/FW {int(peak/FW*100)}%\n"
            f"lines in arm {lines}   clocked {sim.ck_col}/{NF}"
        )
        left, right = [], []
        for i, (ft, zz) in enumerate(FIBERS):
            mk = "→" if i == sim.sel_f else " "
            ln = sim.lines_in_arm(i)
            stars = "★" * min(ln, 3) if ln else "—"
            line = f"{mk}{i:2d} {ft:6s} {zz:4.2f} {stars}"
            (left if i < 10 else right).append(line)
        ftable_left.set_text("\n".join(left))
        ftable_right.set_text("\n".join(right))

    def redraw_all():
        update_arm_bar()
        redraw_ccd()
        redraw_spectrum()
        redraw_stats()
        fig.canvas.draw()

    # ---------- callbacks ----------
    def on_arm(label):
        sim.arm_key = label
        sim.cr.clear(); sim.ck_col = 0
        sim.expose()
        redraw_all()
        info_txt.set_text(f"{sim.arm['name']} exposed.")

    def on_slider(_):
        sim.lstr = int(sliders["lstr"].val)
        sim.cont = int(sliders["cont"].val)
        sim.sky  = int(sliders["sky"].val)
        sim.rn   = int(sliders["rn"].val)

    def on_click(event):
        if event.inaxes is ax_ccd and event.xdata is not None and event.button == 1:
            fi = int(round(event.xdata))
            if 0 <= fi < NF and fi != sim.sel_f:
                sim.sel_f = fi
                redraw_ccd(); redraw_spectrum(); redraw_stats()
                fig.canvas.draw()

    def do_expose(_):
        sim.expose(); redraw_all()
        info_txt.set_text(f"{sim.arm['name']} exposed. Lines shift with each galaxy's z.")

    def do_cr(_):
        sim.add_cr(); redraw_ccd(); fig.canvas.draw()
        info_txt.set_text("Cosmic-ray streak added.")

    def do_reset(_):
        sim.reset(); redraw_all()
        info_txt.set_text("Reset — click Expose.")

    def do_ck1(_):
        sim.clock_one(); redraw_all()
        info_txt.set_text(f"Clocked {sim.ck_col}/{NF} columns.")

    def do_ckall(_):
        if timer[0] is not None:
            timer[0].stop(); timer[0] = None
        t = fig.canvas.new_timer(interval=160)
        def step():
            if not sim.clock_one():
                t.stop(); timer[0] = None
                return
            redraw_all()
            info_txt.set_text(f"Clocking… {sim.ck_col}/{NF}")
        t.add_callback(step)
        timer[0] = t
        t.start()

    def do_save_png(_):
        from datetime import datetime
        stem = f"ccd_{sim.arm_key}_f{sim.sel_f}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        fname = stem + ".png"
        fig.savefig(fname, facecolor=fig.get_facecolor(), dpi=150)
        info_txt.set_text(f"PNG saved → {fname}")
        print(f"saved → {fname}")

    def do_save_data(_):
        from datetime import datetime
        a = sim.arm
        ws = np.arange(NW)
        lams = a["hi"] - ws / NW * (a["hi"] - a["lo"])
        stem = f"ccd_{sim.arm_key}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        # try FITS, fall back to npz + csv
        try:
            from astropy.io import fits
            hdu = fits.PrimaryHDU(sim.grid.astype(np.float32))
            hdu.header["ARM"]   = sim.arm_key
            hdu.header["WMIN"]  = a["lo"]
            hdu.header["WMAX"]  = a["hi"]
            hdu.header["LSTR"]  = sim.lstr
            hdu.header["CONT"]  = sim.cont
            hdu.header["SKY"]   = sim.sky
            hdu.header["RN"]    = sim.rn
            hdu.header["CKCOL"] = sim.ck_col
            fname = stem + ".fits"
            hdu.writeto(fname, overwrite=True)
            info_txt.set_text(f"FITS saved → {fname}")
            print(f"saved → {fname}")
            return
        except ImportError:
            pass
        # fallback: NPZ of grid + CSV of selected spectrum
        npz = stem + ".npz"
        np.savez(npz, grid=sim.grid, wavelengths=lams,
                 arm=sim.arm_key, lstr=sim.lstr, cont=sim.cont,
                 sky=sim.sky, rn=sim.rn, ck_col=sim.ck_col)
        csv = f"{stem}_fiber{sim.sel_f}.csv"
        np.savetxt(csv, np.column_stack([lams, sim.grid[sim.sel_f]]),
                   header="wavelength_nm,counts", delimiter=",", comments="")
        info_txt.set_text(f"saved → {npz}, {csv}")
        print(f"saved → {npz}\nsaved → {csv}")

    arm_radio.on_clicked(on_arm)
    for s in sliders.values(): s.on_changed(on_slider)
    btn_expose.on_clicked(do_expose)
    btn_cr.on_clicked(do_cr)
    btn_reset.on_clicked(do_reset)
    btn_ck1.on_clicked(do_ck1)
    btn_ckall.on_clicked(do_ckall)
    btn_png.on_clicked(do_save_png)
    btn_data.on_clicked(do_save_data)
    fig.canvas.mpl_connect("button_press_event", on_click)

    # keep refs so widgets don't get GC'd
    fig._widgets = (arm_radio, btn_expose, btn_cr, btn_reset,
                    btn_ck1, btn_ckall, btn_png, btn_data, sliders)

    redraw_all()
    plt.show()


if __name__ == "__main__":
    main()
