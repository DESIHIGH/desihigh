"""
desi_focal_plane.py — interactive DESI focal plane / fiber positioner sim.

Three panels:
  1. Full focal plane: 5000 positioners, targets, assignments, collisions
  2. Zoom into a selectable petal — patrol disks + two-arm robot schematics
  3. Completeness-by-class bar chart

Controls (top row):
  Patrol mm · Collision mm · Density × · Petal  · Re-run · Animate
  Save PNG  · Save CSV

"Animate" reveals the greedy assignment one priority class at a time
(QSO → ELG → LRG → BGS → STAR). "Re-run" re-randomises targets and
re-runs the full assignment with current slider values.
"""

import os
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button
from matplotlib.collections import LineCollection, EllipseCollection
from scipy.spatial import KDTree

# ---------- instrument constants ----------
FOCAL_RADIUS_MM = 406.0
FOV_DEG = 3.2
N_PETALS = 10
N_FIBERS = 5000
R_INNER = 1.4
R_OUTER = 1.4
DEFAULT_PATROL = R_INNER + R_OUTER        # 2.8 mm
DEFAULT_COLLISION = 8.0                   # mm
POSITIONER_PITCH = 10.4
ARCSEC_PER_MM = 67.7

TARGET_TYPES = {
    "QSO":  {"priority": 0, "color": "#ff8800", "density": 120},
    "ELG":  {"priority": 1, "color": "#4488ff", "density": 2400},
    "LRG":  {"priority": 2, "color": "#ff4444", "density": 350},
    "BGS":  {"priority": 3, "color": "#88cc44", "density": 800},
    "STAR": {"priority": 4, "color": "#aaaaaa", "density": 200},
}
TYPE_NAMES = list(TARGET_TYPES.keys())


# ---------- geometry / targets ----------
def build_positioner_grid(pitch=POSITIONER_PITCH, n_fibers=N_FIBERS,
                          focal_radius=FOCAL_RADIUS_MM, n_petals=N_PETALS):
    nx = int(2 * focal_radius / pitch) + 4
    xs, ys = [], []
    for row in range(-nx, nx + 1):
        for col in range(-nx, nx + 1):
            x = col * pitch + (row % 2) * pitch * 0.5
            y = row * pitch * np.sqrt(3) / 2
            r = np.hypot(x, y)
            if 50 < r < focal_radius * 0.985:
                xs.append(x); ys.append(y)
    xy = np.array(list(zip(xs, ys)))
    if len(xy) > n_fibers:
        idx = np.argsort(np.hypot(xy[:, 0], xy[:, 1]))[:n_fibers]
        xy = xy[idx]
    angles = np.degrees(np.arctan2(xy[:, 1], xy[:, 0])) % 360.0
    petal = (angles / (360.0 / n_petals)).astype(int)
    return xy, petal


def generate_targets(density_mul, rng, focal_radius=FOCAL_RADIUS_MM):
    fov_area = np.pi * (FOV_DEG / 2) ** 2
    all_xy, types, pri = [], [], []
    for name, tp in TARGET_TYPES.items():
        n = int(tp["density"] * density_mul * fov_area)
        r = focal_radius * np.sqrt(rng.uniform(0, 1, n))
        th = rng.uniform(0, 2 * np.pi, n)
        x = r * np.cos(th); y = r * np.sin(th)
        mask = r > 50
        xy = np.column_stack([x[mask], y[mask]])
        all_xy.append(xy)
        types.extend([name] * len(xy))
        pri.extend([tp["priority"]] * len(xy))
    tgt_xy = np.vstack(all_xy) if all_xy else np.zeros((0, 2))
    return tgt_xy, np.array(types), np.array(pri)


def assign_fibers(pos_xy, tgt_xy, tgt_pri, patrol, collision_dist, rng,
                  priority_filter=None):
    """
    Greedy assignment by priority. If priority_filter is a set of ints,
    only targets with priority in that set are considered.
    Returns (assignment dict, collisions list).
    """
    pos_tree = KDTree(pos_xy)
    assigned = {}
    order = np.lexsort((rng.uniform(size=len(tgt_pri)), tgt_pri))
    for ti in order:
        if priority_filter is not None and tgt_pri[ti] not in priority_filter:
            continue
        tx, ty = tgt_xy[ti]
        reachable = pos_tree.query_ball_point([tx, ty], patrol)
        if not reachable:
            continue
        best, best_d = None, np.inf
        for pi in reachable:
            if pi in assigned:
                continue
            d = np.hypot(pos_xy[pi, 0] - tx, pos_xy[pi, 1] - ty)
            if d < best_d:
                best_d = d; best = pi
        if best is not None:
            assigned[best] = ti

    collisions = []
    idxs = list(assigned.keys())
    if len(idxs) > 1:
        tree = KDTree(pos_xy[idxs])
        for i, j in tree.query_pairs(POSITIONER_PITCH * 2.5):
            pi, pj = idxs[i], idxs[j]
            ti, tj = assigned[pi], assigned[pj]
            if np.hypot(tgt_xy[ti, 0] - tgt_xy[tj, 0],
                        tgt_xy[ti, 1] - tgt_xy[tj, 1]) < collision_dist:
                collisions.append((pi, pj))
    return assigned, collisions


# ---------- state ----------
class Sim:
    def __init__(self):
        self.patrol = DEFAULT_PATROL
        self.collision = DEFAULT_COLLISION
        self.density_mul = 1.0
        self.zoom_petal = 2
        self.rng = np.random.default_rng(1234)
        self.pos_xy, self.pos_petal = build_positioner_grid()
        self.rebuild_targets()
        self.assignment = {}
        self.collisions = []
        self.priority_reveal = set(range(len(TARGET_TYPES)))  # all by default

    def rebuild_targets(self):
        self.tgt_xy, self.tgt_type, self.tgt_pri = generate_targets(
            self.density_mul, self.rng)

    def run_assignment(self, priority_filter=None):
        self.assignment, self.collisions = assign_fibers(
            self.pos_xy, self.tgt_xy, self.tgt_pri,
            self.patrol, self.collision, self.rng,
            priority_filter=priority_filter)


# ---------- plotting ----------
PANEL_BG = "#0a0a18"
TEXT = "#ccccdd"
TICK = "#666688"
SPINE = "#222244"


def style(ax, title=""):
    ax.set_facecolor(PANEL_BG)
    for sp in ax.spines.values(): sp.set_color(SPINE)
    ax.tick_params(colors=TICK, labelcolor=TEXT, labelsize=7)
    if title:
        ax.set_title(title, color=TEXT, fontsize=9, pad=5)
    ax.set_aspect("equal")


def main():
    sim = Sim()
    sim.run_assignment()

    fig = plt.figure(figsize=(16, 8), facecolor="#07070f")
    fig.canvas.manager.set_window_title("DESI focal plane (interactive)")
    gs = GridSpec(1, 3, figure=fig, wspace=0.10,
                  left=0.03, right=0.985, top=0.80, bottom=0.06,
                  width_ratios=[1.0, 1.0, 1.05])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    # ── top-row controls ──
    def mkax(x, y, w, h, fc=PANEL_BG):
        a = fig.add_axes([x, y, w, h]); a.set_facecolor(fc); return a

    def mkslider(x, y, w, h, label, lo, hi, val, step):
        a = mkax(x, y, w, h, "#1a2238")
        s = Slider(a, label, lo, hi, valinit=val, valstep=step, color="#4a79d4")
        s.label.set_color("white"); s.label.set_fontsize(8)
        s.valtext.set_color("white"); s.valtext.set_fontsize(8)
        return s

    def mkbtn(x, y, w, h, label):
        a = mkax(x, y, w, h, "#2a3454")
        b = Button(a, label, color="#2a3454", hovercolor="#3d4a70")
        b.label.set_color("white"); b.label.set_fontsize(8)
        return b

    spd_patrol   = mkslider(0.07, 0.94, 0.08, 0.022, "Patrol",  1.0, 5.0, sim.patrol, 0.1)
    spd_collide  = mkslider(0.21, 0.94, 0.08, 0.022, "Collide", 2.0, 15.0, sim.collision, 0.5)
    spd_density  = mkslider(0.35, 0.94, 0.08, 0.022, "Density", 0.3, 3.0, sim.density_mul, 0.1)
    spd_petal    = mkslider(0.49, 0.94, 0.06, 0.022, "Petal",   0, N_PETALS - 1, sim.zoom_petal, 1)

    btn_rerun    = mkbtn(0.60, 0.935, 0.07, 0.035, "Re-run")
    btn_anim     = mkbtn(0.68, 0.935, 0.07, 0.035, "Animate")
    btn_png      = mkbtn(0.78, 0.935, 0.09, 0.035, "Save PNG")
    btn_csv      = mkbtn(0.88, 0.935, 0.09, 0.035, "Save CSV")

    # status bar
    info_ax = mkax(0.05, 0.905, 0.94, 0.022, "#0f121b")
    info_ax.set_xticks([]); info_ax.set_yticks([])
    info_txt = info_ax.text(0.01, 0.5, "", color="#aaa", fontsize=8, va="center")

    # persistent artists for panel 1 (created once, updated in place)
    p1 = {}

    def setup_panel1():
        ax1.clear()
        style(ax1, "Focal plane")
        for p in range(N_PETALS):
            rad = np.radians(p * 360 / N_PETALS)
            ax1.plot([0, FOCAL_RADIUS_MM * np.cos(rad)],
                     [0, FOCAL_RADIUS_MM * np.sin(rad)],
                     color="#222244", lw=0.5)
        th = np.linspace(0, 2 * np.pi, 361)
        ax1.plot(FOCAL_RADIUS_MM * np.cos(th), FOCAL_RADIUS_MM * np.sin(th),
                 color="#334466", lw=0.9)
        ax1.add_patch(plt.Circle((0, 0), 50, color=PANEL_BG, zorder=2))
        ax1.set_xlim(-FOCAL_RADIUS_MM * 1.05, FOCAL_RADIUS_MM * 1.05)
        ax1.set_ylim(-FOCAL_RADIUS_MM * 1.05, FOCAL_RADIUS_MM * 1.05)
        ax1.set_xlabel("focal x (mm)", color=TEXT, fontsize=8)
        ax1.set_ylabel("focal y (mm)", color=TEXT, fontsize=8)

        # one rasterised scatter per target type (updated in refresh_targets)
        p1["tgt_scatters"] = {}
        for name, tp in TARGET_TYPES.items():
            sc = ax1.scatter([], [], s=0.4, color=tp["color"], alpha=0.25,
                             rasterized=True, zorder=3)
            p1["tgt_scatters"][name] = sc

        # single scatter for all positioners; we'll recolor each frame
        p1["pos_scatter"] = ax1.scatter(sim.pos_xy[:, 0], sim.pos_xy[:, 1],
                                        s=1.4, c="#333355", alpha=0.85,
                                        rasterized=True, zorder=5)
        # two line collections: assignment lines + collision lines
        p1["asgn_lines"]  = LineCollection([], colors="#888", linewidths=0.3,
                                           alpha=0.45, zorder=4)
        p1["coll_lines"]  = LineCollection([], colors="#ff0000", linewidths=0.9,
                                           alpha=0.6, zorder=6)
        ax1.add_collection(p1["asgn_lines"])
        ax1.add_collection(p1["coll_lines"])
        p1["legend"] = None

    def update_panel1():
        # targets
        for name, sc in p1["tgt_scatters"].items():
            m = sim.tgt_type == name
            sc.set_offsets(sim.tgt_xy[m])

        # positioner colors
        colors = np.tile(np.array([0.2, 0.2, 0.33, 0.85]), (len(sim.pos_xy), 1))
        for pi, ti in sim.assignment.items():
            c = TARGET_TYPES[sim.tgt_type[ti]]["color"]
            r = int(c[1:3], 16) / 255; g = int(c[3:5], 16) / 255
            b = int(c[5:7], 16) / 255
            colors[pi] = [r, g, b, 0.9]
        p1["pos_scatter"].set_facecolor(colors)

        # assignment lines (as segments)
        segs = np.array([
            [[sim.pos_xy[pi, 0], sim.pos_xy[pi, 1]],
             [sim.tgt_xy[ti, 0], sim.tgt_xy[ti, 1]]]
            for pi, ti in sim.assignment.items()
        ]) if sim.assignment else np.zeros((0, 2, 2))
        line_colors = [TARGET_TYPES[sim.tgt_type[ti]]["color"]
                       for _, ti in sim.assignment.items()]
        p1["asgn_lines"].set_segments(segs)
        if line_colors:
            p1["asgn_lines"].set_color(line_colors)

        # collision lines (cap at 100 so we never explode the artist count)
        csegs = np.array([
            [[sim.pos_xy[pi, 0], sim.pos_xy[pi, 1]],
             [sim.pos_xy[pj, 0], sim.pos_xy[pj, 1]]]
            for pi, pj in sim.collisions[:100]
        ]) if sim.collisions else np.zeros((0, 2, 2))
        p1["coll_lines"].set_segments(csegs)

        ax1.set_title(f"Focal plane  —  {len(sim.pos_xy)} positioners, "
                      f"{len(sim.tgt_xy):,} targets", color=TEXT, fontsize=9,
                      pad=18)

        if p1["legend"] is not None:
            p1["legend"].remove()
        handles = [mpatches.Patch(color=tp["color"], label=tn)
                   for tn, tp in TARGET_TYPES.items()]
        handles += [mpatches.Patch(color="#333355", label="unassigned"),
                    mpatches.Patch(color="#ff0000",
                                   label=f"collision ({len(sim.collisions)})")]
        # legend above the focal plane, horizontal
        p1["legend"] = ax1.legend(
            handles=handles, loc="lower center",
            bbox_to_anchor=(0.5, 1.00), ncol=len(handles),
            fontsize=6, facecolor=PANEL_BG, labelcolor=TEXT,
            framealpha=0.85, handlelength=1.0, columnspacing=1.0,
            borderpad=0.3)

    # persistent artists for panel 2 (no EllipseCollection — fragile on py3.8 mpl)
    p2 = {"inited": False}

    def setup_panel2():
        ax2.clear()
        style(ax2, "Zoom — petal")
        ax2.set_xlabel("focal x (mm)", color=TEXT, fontsize=8)
        # arms + assignment lines
        p2["arms"] = LineCollection([], linewidths=0.8, alpha=0.7)
        ax2.add_collection(p2["arms"])
        # positioner scatter (colored by assignment)
        p2["pos"] = ax2.scatter([], [], s=6)
        # assigned target scatter
        p2["tgt_assigned"] = ax2.scatter([], [], s=10)
        # background targets per type
        p2["tgt_bg"] = {}
        for name, tp in TARGET_TYPES.items():
            p2["tgt_bg"][name] = ax2.scatter(
                [], [], s=3, color=tp["color"], alpha=0.3)
        p2["inited"] = True

    def update_panel2():
        if not p2["inited"]:
            setup_panel2()
        pmask = sim.pos_petal == sim.zoom_petal
        petal_pos = sim.pos_xy[pmask]
        if len(petal_pos) == 0:
            return
        petal_idx = np.where(pmask)[0]

        assigned_mask = np.array([pi in sim.assignment for pi in petal_idx])
        ua_pos = petal_pos[~assigned_mask]
        a_pos  = petal_pos[assigned_mask]
        a_idxs = petal_idx[assigned_mask]

        segs = []; arm_colors = []
        tgt_xys = []; tgt_cols = []
        for pi in a_idxs:
            ti = sim.assignment[pi]
            col = TARGET_TYPES[sim.tgt_type[ti]]["color"]
            px, py = sim.pos_xy[pi]; tx, ty = sim.tgt_xy[ti]
            ex = (px + tx) * 0.5; ey = (py + ty) * 0.5
            segs.append([[px, py], [ex, ey]])
            segs.append([[ex, ey], [tx, ty]])
            arm_colors.extend([col, col])
            tgt_xys.append([tx, ty]); tgt_cols.append(col)
        p2["arms"].set_segments(segs)
        if arm_colors:
            p2["arms"].set_color(arm_colors)

        pos_xys = np.vstack([ua_pos, a_pos]) if len(ua_pos) or len(a_pos) else np.zeros((0, 2))
        pos_fc = (["#333355"] * len(ua_pos)) + (["#ffffff"] * len(a_pos))
        p2["pos"].set_offsets(pos_xys)
        if pos_fc:
            p2["pos"].set_facecolor(pos_fc)

        p2["tgt_assigned"].set_offsets(
            np.asarray(tgt_xys) if tgt_xys else np.zeros((0, 2)))
        if tgt_cols:
            p2["tgt_assigned"].set_facecolor(tgt_cols)

        xmin = petal_pos[:, 0].min() - 20; xmax = petal_pos[:, 0].max() + 20
        ymin = petal_pos[:, 1].min() - 20; ymax = petal_pos[:, 1].max() + 20
        ax2.set_xlim(xmin, xmax); ax2.set_ylim(ymin, ymax)
        inview = ((sim.tgt_xy[:, 0] > xmin) & (sim.tgt_xy[:, 0] < xmax) &
                  (sim.tgt_xy[:, 1] > ymin) & (sim.tgt_xy[:, 1] < ymax))
        for name, sc in p2["tgt_bg"].items():
            m = inview & (sim.tgt_type == name)
            sc.set_offsets(sim.tgt_xy[m])

        ax2.set_title(f"Zoom — petal {sim.zoom_petal}",
                      color=TEXT, fontsize=9, pad=5)

    def draw_panel3():
        ax3.clear()
        ax3.set_facecolor(PANEL_BG)
        for sp in ax3.spines.values(): sp.set_color(SPINE)
        ax3.tick_params(colors=TICK, labelcolor=TEXT, labelsize=8)
        ax3.set_title("Assignment completeness", color=TEXT, fontsize=9, pad=5)

        avail = np.array([int((sim.tgt_type == t).sum()) for t in TYPE_NAMES])
        asgnd = np.zeros(len(TYPE_NAMES), int)
        for _, ti in sim.assignment.items():
            asgnd[TYPE_NAMES.index(sim.tgt_type[ti])] += 1
        colors = [TARGET_TYPES[t]["color"] for t in TYPE_NAMES]
        y = np.arange(len(TYPE_NAMES))
        ax3.barh(y, avail, height=0.55, color="#1a1a2e")
        ax3.barh(y, asgnd, height=0.55, color=colors, alpha=0.85)
        max_a = max(avail.max(), 1)
        for i, (a, v) in enumerate(zip(asgnd, avail)):
            pct = (a / v * 100) if v else 0
            ax3.text(a + max_a * 0.01, i, f"{a:,}  ({pct:.0f}%)",
                     va="center", ha="left", color=TEXT, fontsize=8)
        ax3.set_yticks(y); ax3.set_yticklabels(TYPE_NAMES, color=TEXT)
        ax3.set_xlabel("Fibers assigned", color=TEXT, fontsize=8)
        ax3.set_xlim(0, max_a * 1.35); ax3.invert_yaxis()
        ax3.grid(axis="x", color="#1a1a2e", lw=0.4)

        summary = (f"patrol {sim.patrol:.1f} mm · collide {sim.collision:.1f} mm · "
                   f"density {sim.density_mul:.1f}×\n"
                   f"assigned {len(sim.assignment):,}/{len(sim.pos_xy):,} "
                   f"({100*len(sim.assignment)/len(sim.pos_xy):.0f}%) · "
                   f"collisions {len(sim.collisions)}")
        ax3.text(0.97, 0.03, summary, transform=ax3.transAxes, ha="right",
                 va="bottom", color=TEXT, fontsize=7.5,
                 bbox=dict(facecolor=PANEL_BG, edgecolor=SPINE, alpha=0.9, pad=4))

    def redraw_all(use_idle=True):
        update_panel1(); update_panel2(); draw_panel3()
        info_txt.set_text(
            f"assigned {len(sim.assignment):,}/{len(sim.pos_xy):,}  "
            f"· collisions {len(sim.collisions)}  · targets {len(sim.tgt_xy):,}")
        fig.canvas.draw_idle()

    # --- callbacks ---
    def on_slider(_):
        sim.patrol       = float(spd_patrol.val)
        sim.collision    = float(spd_collide.val)
        sim.density_mul  = float(spd_density.val)
        sim.zoom_petal   = int(spd_petal.val)
    for s in (spd_patrol, spd_collide, spd_density, spd_petal):
        s.on_changed(on_slider)
    # petal change = just redraw panel 2 without re-running assignment
    def on_petal(_):
        sim.zoom_petal = int(spd_petal.val)
        update_panel2(); fig.canvas.draw_idle()
    spd_petal.on_changed(on_petal)

    def on_rerun(_):
        stop_anim()
        info_txt.set_text("re-running assignment…"); fig.canvas.draw_idle()
        sim.rebuild_targets()
        sim.run_assignment()
        redraw_all()

    def on_save_png(_):
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
            f"desi_focal_plane_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png")
        fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
        info_txt.set_text(f"saved → {out}"); fig.canvas.draw_idle()
        print(f"saved → {out}")

    def on_save_csv(_):
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
            f"desi_focal_plane_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
        with open(out, "w") as f:
            f.write("pos_idx,pos_x_mm,pos_y_mm,petal,assigned,target_idx,"
                    "target_type,target_x_mm,target_y_mm\n")
            for i in range(len(sim.pos_xy)):
                px, py = sim.pos_xy[i]; petal = sim.pos_petal[i]
                if i in sim.assignment:
                    ti = sim.assignment[i]; tx, ty = sim.tgt_xy[ti]
                    f.write(f"{i},{px:.3f},{py:.3f},{petal},1,{ti},"
                            f"{sim.tgt_type[ti]},{tx:.3f},{ty:.3f}\n")
                else:
                    f.write(f"{i},{px:.3f},{py:.3f},{petal},0,,,,\n")
        info_txt.set_text(f"saved → {out}"); fig.canvas.draw_idle()
        print(f"saved → {out}")

    # --- animation of assignment reveal (synchronous — safe on py3.8 macOS) ---
    anim_state = {"running": False}

    def stop_anim():
        anim_state["running"] = False

    def on_anim(_):
        if anim_state["running"]:
            return
        anim_state["running"] = True
        sim.assignment = {}; sim.collisions = []
        info_txt.set_text("animating assignment…")
        redraw_all()
        plt.pause(0.3)   # show empty state briefly
        for k in range(len(TYPE_NAMES)):
            if not anim_state["running"]:
                break
            sim.run_assignment(priority_filter=set(range(k + 1)))
            info_txt.set_text(
                f"revealed priority 0..{k} ({TYPE_NAMES[k]}) — "
                f"{len(sim.assignment):,} assigned")
            redraw_all()
            plt.pause(0.8)   # paints the figure AND sleeps — key difference
        anim_state["running"] = False
        info_txt.set_text(
            f"animation complete — {len(sim.assignment):,}/"
            f"{len(sim.pos_xy):,} assigned, "
            f"{len(sim.collisions)} collisions")
        fig.canvas.draw_idle()

    btn_rerun.on_clicked(on_rerun)
    btn_anim.on_clicked(on_anim)
    btn_png.on_clicked(on_save_png)
    btn_csv.on_clicked(on_save_csv)

    setup_panel1()
    setup_panel2()

    fig._keep = (spd_patrol, spd_collide, spd_density, spd_petal,
                 btn_rerun, btn_anim, btn_png, btn_csv, anim_state)

    redraw_all()
    plt.show()


if __name__ == "__main__":
    main()
