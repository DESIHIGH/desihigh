"""
Mayall 4m + DESI corrector + focal plane — interactive cartoon ray-trace.

Geometry is to-scale (proportions accurate to published specs) but the
ray propagation is geometric, not a physical trace through each lens
element.  For a faithful trace, use `batoid` with the DESI prescription.

What's shown:
  - Mayall primary mirror  (4.0 m dia, f/2.7 paraboloid)
  - Serrurier truss telescope tube + top/bottom end rings
  - Prime-focus cage
  - DESI corrector (6 lens elements + barrel)
  - Curved focal plane with fiber-positioner pencils
  - Equatorial horseshoe-yoke mount at Kitt Peak latitude (31.96°)
  - Five field-angle ray bundles fanning across the 3.2° FOV

Controls:
  - Zenith angle slider — tilt telescope in elevation (mount follows)
  - Focal plane zoom button — zoom into focal plane region
  - Reset view button — back to full telescope view
  - Toggle mount button — hide/show mount silhouette
  - Save PNG button — timestamped snapshot
  - Click any ray bundle in the legend to isolate

Run in Jupyter:
    %matplotlib tk
    %run /path/to/mayall_desi_raytrace.py

Or from a terminal:
    python mayall_desi_raytrace.py
"""

import os
from datetime import datetime

import numpy as np
import matplotlib
# The MacOSX backend's navigation toolbar can swallow widget click events
# on some OS versions.  Disable it before we create the figure.
matplotlib.rcParams['toolbar'] = 'None'
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Button, Slider
from matplotlib.lines import Line2D
from matplotlib.transforms import Affine2D


# ==========================================================================
# Geometry — all in meters, telescope axis along +y
# ==========================================================================

# Primary mirror (Mayall 4m, f/2.7 paraboloid)
PRIMARY_D   = 4.0                 # diameter (m)
PRIMARY_F   = 10.8                # focal length (m) — f/2.7 × 4m
PRIMARY_R   = PRIMARY_D / 2       # radius (m)
PRIMARY_Y0  = 0.0                 # vertex at origin

# Parabola: y = x² / (4f)
def primary_sag(x):
    return x**2 / (4 * PRIMARY_F)

# Effective system (after corrector): f/3.9 → EFL ≈ 15.6 m
EFF_F       = 15.6                # effective focal length after corrector
FOV_DEG     = 3.2                 # full field of view (deg)
FOV_HALF    = FOV_DEG / 2         # half-FOV (deg)

# DESI corrector barrel — converging beam from primary enters at the
# BOTTOM (C1, widest lens) and exits at the TOP where the focal plane
# sits, just above the last lens (C4).
CORR_BOT_Y  = PRIMARY_F - 2.0     # bottom (C1 — entrance, beam still wide)
CORR_TOP_Y  = PRIMARY_F + 0.1     # top (C4 — just below focal plane)
CORR_R_BOT  = 0.55                # C1 is biggest (~1.1 m dia, beam widest here)
CORR_R_TOP  = 0.42                # narrows as beam converges toward focus
FOCAL_Y     = CORR_TOP_Y + 0.18   # focal plane sits just above last lens (C4)
FOCAL_R     = 0.40                # 0.8 m dia curved focal surface
FOCAL_SAG   = 0.045               # ~45 mm sag on focal surface

# Six DESI lens elements — relative z positions from BOTTOM of barrel.
# Light enters at C1 (bottom), passes up through C2, ADC1, ADC2, C3, C4,
# and comes to focus on the focal plane just above C4.
LENS_Z_FRAC = [0.05, 0.22, 0.38, 0.50, 0.70, 0.88]  # from bottom (C1 → C4)
LENS_R_FRAC = [1.00, 0.92, 0.85, 0.85, 0.88, 0.95]  # relative lens radius
LENS_SHAPE  = ['biconvex', 'meniscus_up', 'biconvex',
               'biconvex', 'meniscus_down', 'meniscus_up']

# Serrurier truss tube
TUBE_TOP_Y     = PRIMARY_F + 1.35  # top end ring (above focal plane)
TUBE_BOT_Y     = 0.6               # bottom end ring (above primary cell)
TUBE_TOP_R     = 2.1
TUBE_BOT_R     = 2.25

# Prime-focus cage (houses corrector + focal plane + readout electronics)
CAGE_TOP_Y     = TUBE_TOP_Y + 0.1
CAGE_BOT_Y     = CORR_BOT_Y - 0.15
CAGE_R         = CORR_R_BOT + 0.08   # must clear the widest lens (C1)

# Mount (equatorial horseshoe-yoke at Kitt Peak, lat ≈ 31.96°)
LATITUDE_DEG   = 31.96
POLAR_TILT_RAD = np.deg2rad(LATITUDE_DEG)   # polar axis from horizontal


# ==========================================================================
# Palette — dark theme matching the H0 figure
# ==========================================================================

BG       = '#050612'       # same as C_SKY so the sky fills the whole page
PANEL_BG = '#0f101d'
TEXT     = '#d8dce8'
TEXT_DIM = '#9aa0b0'
SPINE    = '#3a3e52'
GRID     = '#1c1e30'

C_MIRROR   = '#8fb9d9'
C_TUBE     = '#4a5070'
C_CAGE     = '#6b7296'
C_LENS     = '#a8d8ff'
C_LENS_EDGE= '#d8f0ff'
C_FOCAL    = '#f4d03f'
C_FIBER    = '#ff6ec7'   # bright magenta — where the light lands
C_MOUNT    = '#3a3f5a'
C_PIER     = '#2a2e42'
C_SKY      = '#050612'

# Ray bundle colours — five distinct field angles
RAY_COLORS = ['#76D7C4', '#4E9FD1', '#F4D03F', '#F28E2B', '#E15759']
FIELD_ANGLES = np.array([-1.6, -0.8, 0.0, 0.8, 1.6])   # deg
N_RAYS_PER_BUNDLE = 4   # marginal rays per field angle


# ==========================================================================
# Ray tracing (cartoon geometry)
# ==========================================================================

def reflect_off_primary(x_hit, theta_in_rad):
    """
    Given a ray hitting the primary at x = x_hit with incoming direction
    (sin θ, -cos θ), return its reflected direction unit vector.
    Primary surface normal at (x, y=x²/4f) points toward (-x, 2f) normalised.
    """
    # surface normal (pointing up-and-inward)
    nx = -x_hit
    ny = 2 * PRIMARY_F
    n = np.array([nx, ny]) / np.hypot(nx, ny)

    # incoming ray direction
    d_in = np.array([np.sin(theta_in_rad), -np.cos(theta_in_rad)])

    # reflection: d_out = d_in - 2(d_in·n)n
    d_out = d_in - 2 * np.dot(d_in, n) * n
    return d_out

def trace_bundle(field_angle_deg, n_rays=4):
    """
    Return a list of (xs, ys) polylines — one per ray in the bundle —
    for a parallel beam entering at field_angle_deg from vertical.

    Cartoon model:
      1. Parallel rays come from above.
      2. Reflect off parabolic primary (true law of reflection).
      3. Converge to a point on the curved focal surface.  The focal
         spot position is set by the effective FL (f/3.9), not the bare
         prime focus, so the rays "bend" through the corrector region
         to reach their corrected focus.
    """
    theta = np.deg2rad(field_angle_deg)

    # Marginal-ray entry positions — spread across the aperture, but
    # pulled in a bit so that at the steepest field angle (±1.6°) the
    # outermost rays still enter the tube-top ring (radius TUBE_TOP_R).
    x_entries = np.linspace(-PRIMARY_R * 0.85, PRIMARY_R * 0.85, n_rays)

    # Final focus spot on the curved focal surface
    x_focus = EFF_F * np.tan(theta)
    # Curved focal surface: y = FOCAL_Y + FOCAL_SAG * (x/FOCAL_R)²  (concave up)
    y_focus = FOCAL_Y + FOCAL_SAG * (x_focus / FOCAL_R)**2

    rays = []
    for x0 in x_entries:
        # Position where this ray crosses the tube-TOP ring.  If it would
        # be outside the aperture the telescope never sees it — skip.
        y_hit = primary_sag(x0)
        x_at_top = x0 + np.tan(theta) * (TUBE_TOP_Y - y_hit)
        if abs(x_at_top) > TUBE_TOP_R:
            continue

        # Start high above the telescope
        y_start = TUBE_TOP_Y + 3.0
        # direction (sin θ, -cos θ): dx/dy = -tan θ  → dx = -tan(θ) * dy
        dy = y_hit - y_start
        dx = -np.tan(theta) * dy
        x_start = x0 - dx

        # Reflected direction from true parabolic reflection
        d_out = reflect_off_primary(x0, theta)

        # Propagate reflected ray UP until it enters the corrector at the
        # bottom (C1).  Beyond that we cartoon the corrector as a straight
        # shot to the focal spot — the actual bending happens across
        # 6 refractive surfaces, which we don't model.
        t_entry = (CORR_BOT_Y - y_hit) / d_out[1]
        x_corr_entry = x0 + t_entry * d_out[0]
        y_corr_entry = CORR_BOT_Y

        xs = [x_start, x0, x_corr_entry, x_focus]
        ys = [y_start, y_hit, y_corr_entry, y_focus]
        rays.append((np.array(xs), np.array(ys)))

    return rays, (x_focus, y_focus)


# ==========================================================================
# Figure
# ==========================================================================


# ======================================================================
# Entry point
# ======================================================================

def main():
    fig = plt.figure(figsize=(13, 8.5), facecolor=BG)
    try:
        fig.canvas.manager.set_window_title('Mayall + DESI — cartoon ray-trace')
    except Exception:
        pass

    gs = fig.add_gridspec(1, 1, left=0.06, right=0.78, top=0.92, bottom=0.27)
    ax = fig.add_subplot(gs[0])
    ax.set_facecolor(C_SKY)
    for s in ax.spines.values():
        s.set_color(SPINE)
    ax.tick_params(colors=TEXT_DIM)

    fig.suptitle('The Mayall 4 m with DESI — prime-focus corrector ray path',
                 fontsize=14, fontweight='bold', color=TEXT, y=0.955)


    # ==========================================================================
    # Mount (drawn first so it sits behind the telescope)
    # ==========================================================================

    mount_artists = []

    def draw_mount():
        """Equatorial horseshoe-yoke mount, drawn at Kitt Peak latitude.
           Polar axis tilted ~32° from horizontal.  Tube assumed vertical
           (telescope pointing near zenith) — we show the mount in the pose
           it would be in for a zenith pointing."""
        tilt = POLAR_TILT_RAD

        # Pier / concrete base
        pier = mpatches.Polygon(
            [[-3.5, -4.0], [3.5, -4.0], [2.8, -2.2], [-2.8, -2.2]],
            facecolor=C_PIER, edgecolor=SPINE, linewidth=1.0, zorder=0)
        ax.add_patch(pier)
        mount_artists.append(pier)

        # North pier bearing block
        np_block = mpatches.Rectangle((-3.2, -2.2), 2.0, 1.4,
                                      facecolor=C_MOUNT, edgecolor=SPINE,
                                      linewidth=1.0, zorder=0)
        ax.add_patch(np_block)
        mount_artists.append(np_block)

        # South pier bearing block
        sp_block = mpatches.Rectangle((1.2, -2.2), 2.0, 1.4,
                                      facecolor=C_MOUNT, edgecolor=SPINE,
                                      linewidth=1.0, zorder=0)
        ax.add_patch(sp_block)
        mount_artists.append(sp_block)

        # Polar axis — long shaft tilted at latitude
        polar_len = 8.5
        px0, py0 = -3.5, -0.5       # south end (near pier top, lower in view)
        px1 = px0 + polar_len * np.cos(tilt)
        py1 = py0 + polar_len * np.sin(tilt)
        polar = Line2D([px0, px1], [py0, py1], color=C_MOUNT, linewidth=6,
                       solid_capstyle='round', zorder=0)
        ax.add_line(polar)
        mount_artists.append(polar)

        # Horseshoe bearing at the north end of the polar axis
        # (large open ring)
        hs_cx, hs_cy = px1, py1
        hs_r = 2.6
        horseshoe = mpatches.Annulus(
            (hs_cx, hs_cy), hs_r, 0.35, angle=np.rad2deg(tilt),
            facecolor=C_MOUNT, edgecolor=SPINE, linewidth=1.2, zorder=0)
        ax.add_patch(horseshoe)
        mount_artists.append(horseshoe)

        # Yoke / fork arms straddling the tube — two beams from horseshoe
        # to the declination trunnions on the tube
        for side in (-1, 1):
            trunnion_x = side * TUBE_TOP_R * 1.05
            trunnion_y = (TUBE_BOT_Y + TUBE_TOP_Y) / 2 - 1.5
            # point on horseshoe circumference nearest trunnion
            arm = Line2D([hs_cx + side * hs_r * 0.6, trunnion_x],
                         [hs_cy + hs_r * 0.3, trunnion_y],
                         color=C_MOUNT, linewidth=5, solid_capstyle='round',
                         zorder=0)
            ax.add_line(arm)
            mount_artists.append(arm)
            # declination trunnion bearing
            tr = mpatches.Circle((trunnion_x, trunnion_y), 0.35,
                                 facecolor=C_MOUNT, edgecolor=SPINE,
                                 linewidth=1.0, zorder=0)
            ax.add_patch(tr)
            mount_artists.append(tr)

        # Latitude annotation — place it low-right, clear of the tube,
        # with a short leader line pointing at the polar axis
        axis_mid_x = (px0 + px1) / 2
        axis_mid_y = (py0 + py1) / 2
        label_x    = 6.2
        label_y    = -3.2
        leader = Line2D([axis_mid_x, label_x - 0.1],
                        [axis_mid_y, label_y + 0.1],
                        color=TEXT_DIM, linewidth=0.6, alpha=0.7, zorder=0)
        ax.add_line(leader)
        mount_artists.append(leader)
        lab = ax.text(label_x, label_y,
                      f'Polar axis\n(tilted {LATITUDE_DEG:.1f}° — Kitt Peak latitude)',
                      color=TEXT_DIM, fontsize=8.5,
                      ha='left', va='top', zorder=1)
        mount_artists.append(lab)


    # ==========================================================================
    # Telescope tube (Serrurier truss)
    # ==========================================================================

    def draw_tube():
        # Top end ring (supports prime focus cage)
        ring_top = mpatches.Rectangle(
            (-TUBE_TOP_R, TUBE_TOP_Y - 0.08), 2 * TUBE_TOP_R, 0.16,
            facecolor=C_TUBE, edgecolor=SPINE, linewidth=1.0, zorder=1)
        ax.add_patch(ring_top)

        # Bottom end ring (at top of primary cell)
        ring_bot = mpatches.Rectangle(
            (-TUBE_BOT_R, TUBE_BOT_Y - 0.08), 2 * TUBE_BOT_R, 0.16,
            facecolor=C_TUBE, edgecolor=SPINE, linewidth=1.0, zorder=1)
        ax.add_patch(ring_bot)

        # Center ring (mid-tube, where declination trunnions attach)
        mid_y = (TUBE_TOP_Y + TUBE_BOT_Y) / 2 - 1.5
        mid_r = (TUBE_TOP_R + TUBE_BOT_R) / 2 + 0.05
        ring_mid = mpatches.Rectangle(
            (-mid_r, mid_y - 0.1), 2 * mid_r, 0.2,
            facecolor=C_TUBE, edgecolor=SPINE, linewidth=1.0, zorder=1)
        ax.add_patch(ring_mid)

        # Serrurier truss struts — crossed pairs, upper (top ring to mid)
        # and lower (mid to bottom ring)
        for side in (-1, 1):
            # Upper truss (X pattern)
            ax.plot([side * TUBE_TOP_R, side * mid_r * 0.95],
                    [TUBE_TOP_Y - 0.08, mid_y + 0.1],
                    color=C_TUBE, linewidth=1.2, zorder=1)
            ax.plot([side * TUBE_TOP_R, -side * mid_r * 0.1],
                    [TUBE_TOP_Y - 0.08, mid_y + 0.1],
                    color=C_TUBE, linewidth=1.0, zorder=1, alpha=0.7)
            # Lower truss
            ax.plot([side * mid_r * 0.95, side * TUBE_BOT_R],
                    [mid_y - 0.1, TUBE_BOT_Y + 0.08],
                    color=C_TUBE, linewidth=1.2, zorder=1)
            ax.plot([side * mid_r * 0.95, -side * TUBE_BOT_R * 0.1],
                    [mid_y - 0.1, TUBE_BOT_Y + 0.08],
                    color=C_TUBE, linewidth=1.0, zorder=1, alpha=0.7)

        # Primary mirror cell (below bottom ring)
        cell = mpatches.Rectangle(
            (-TUBE_BOT_R - 0.15, -0.6), 2 * (TUBE_BOT_R + 0.15), 1.1,
            facecolor='#1a1d30', edgecolor=SPINE, linewidth=1.0, zorder=1)
        ax.add_patch(cell)


    # ==========================================================================
    # Prime-focus cage + corrector
    # ==========================================================================

    def draw_cage_and_corrector():
        # Cage outline
        cage = mpatches.FancyBboxPatch(
            (-CAGE_R, CAGE_BOT_Y), 2 * CAGE_R, CAGE_TOP_Y - CAGE_BOT_Y,
            boxstyle='round,pad=0.02,rounding_size=0.05',
            facecolor='#141626', edgecolor=C_CAGE, linewidth=1.2, zorder=2)
        ax.add_patch(cage)

        # Corrector barrel — tapered, wide at bottom (entrance from primary)
        # and narrow at top (exit to focal plane)
        barrel = mpatches.Polygon(
            [[-CORR_R_BOT, CORR_BOT_Y],
             [ CORR_R_BOT, CORR_BOT_Y],
             [ CORR_R_TOP, CORR_TOP_Y],
             [-CORR_R_TOP, CORR_TOP_Y]],
            facecolor='#1a1d30', edgecolor=C_CAGE, linewidth=1.0, zorder=3)
        ax.add_patch(barrel)

        # Six lens elements — frac_z is fraction from BOTTOM of barrel
        # (C1 near bottom, C4 near top)
        for frac_z, frac_r, shape in zip(LENS_Z_FRAC, LENS_R_FRAC, LENS_SHAPE):
            y = CORR_BOT_Y + frac_z * (CORR_TOP_Y - CORR_BOT_Y)
            # linearly interpolate barrel radius at this height (wide at bottom)
            t = (CORR_TOP_Y - y) / (CORR_TOP_Y - CORR_BOT_Y)   # 1 at bottom, 0 at top
            r_here = CORR_R_TOP * (1 - t) + CORR_R_BOT * t
            r_lens = r_here * frac_r * 0.92

            # Lens profile — cosmetic, just to hint shape
            if shape == 'biconvex':
                top = 0.04; bot = -0.04
            elif shape == 'meniscus_up':
                top = 0.05; bot = 0.02
            elif shape == 'meniscus_down':
                top = -0.02; bot = -0.05
            else:
                top = 0.03; bot = -0.03

            xs_lens = np.linspace(-r_lens, r_lens, 50)
            # top surface
            ytop = y + top * np.cos(np.pi * xs_lens / (2 * r_lens))
            # bottom surface
            ybot = y + bot * np.cos(np.pi * xs_lens / (2 * r_lens))
            verts = list(zip(xs_lens, ytop)) + list(zip(xs_lens[::-1], ybot[::-1]))
            lens = mpatches.Polygon(verts, facecolor=C_LENS, alpha=0.35,
                                    edgecolor=C_LENS_EDGE, linewidth=0.9,
                                    zorder=4)
            ax.add_patch(lens)

        # Curved focal plane (sits ABOVE the last corrector lens; edges curve
        # upward, away from the corrector)
        xs_f = np.linspace(-FOCAL_R, FOCAL_R, 80)
        ys_f = FOCAL_Y + FOCAL_SAG * (xs_f / FOCAL_R)**2
        ax.plot(xs_f, ys_f, color=C_FOCAL, linewidth=2.6, solid_capstyle='round',
                zorder=5)

        # Fiber tips — SHORT magenta stubs hanging down from the focal plane;
        # these mark where the light actually lands (the entrance of each
        # fiber).  DESI is nearly telecentric, so they follow the local
        # focal-surface normal.
        for x_fp in np.linspace(-FOCAL_R * 0.9, FOCAL_R * 0.9, 11):
            y_fp = FOCAL_Y + FOCAL_SAG * (x_fp / FOCAL_R)**2
            slope = 2 * FOCAL_SAG * x_fp / FOCAL_R**2
            pencil_len = 0.035                       # much shorter fiber tips
            dx = pencil_len * np.sin(np.arctan(slope))
            dy = -pencil_len * np.cos(np.arctan(slope))
            ax.plot([x_fp, x_fp + dx], [y_fp, y_fp + dy],
                    color=C_FIBER, linewidth=2.2, alpha=0.95,
                    solid_capstyle='round', zorder=6)
            # Bright dot at the tip — the light-entry point
            ax.plot(x_fp + dx, y_fp + dy, 'o',
                    color=C_FIBER, markersize=2.6, zorder=7)


    # ==========================================================================
    # Primary mirror
    # ==========================================================================

    def draw_primary():
        xs = np.linspace(-PRIMARY_R, PRIMARY_R, 120)
        ys = primary_sag(xs)
        # Fill below the parabola to suggest the mirror substrate
        verts = list(zip(xs, ys)) + [(PRIMARY_R, -0.5), (-PRIMARY_R, -0.5)]
        substrate = mpatches.Polygon(verts, facecolor='#1a1d30',
                                     edgecolor='none', zorder=1)
        ax.add_patch(substrate)
        # Reflective surface highlight
        ax.plot(xs, ys, color=C_MIRROR, linewidth=2.8, solid_capstyle='round',
                zorder=2)
        # Subtle inner gradient band
        ax.plot(xs, ys + 0.02, color=C_MIRROR, linewidth=0.8, alpha=0.4,
                zorder=2)
        # Central hole (Cassegrain hole, unused at prime focus but there)
        hole_r = 0.5
        xs_h = np.linspace(-hole_r, hole_r, 30)
        ys_h = primary_sag(xs_h)
        ax.plot(xs_h, ys_h, color=C_SKY, linewidth=3.2, solid_capstyle='round',
                zorder=2)


    # ==========================================================================
    # Background stars (just a few, for atmosphere)
    # ==========================================================================

    def draw_stars():
        rng = np.random.default_rng(42)
        # Small faint stars
        n = 140
        xs = rng.uniform(-11, 11, n)
        ys = rng.uniform(TUBE_TOP_Y + 0.5, TUBE_TOP_Y + 8.5, n)
        sizes = rng.uniform(0.3, 4.0, n)
        ax.scatter(xs, ys, s=sizes, color='#d8dce8', alpha=0.55,
                   zorder=0, linewidths=0)
        # A few brighter ones
        n2 = 18
        xs = rng.uniform(-11, 11, n2)
        ys = rng.uniform(TUBE_TOP_Y + 0.5, TUBE_TOP_Y + 8.5, n2)
        ax.scatter(xs, ys, s=rng.uniform(8, 22, n2), color='#fff4c8',
                   alpha=0.85, zorder=0, linewidths=0, marker='*')


    # ==========================================================================
    # Build the scene
    # ==========================================================================

    draw_stars()
    draw_mount()

    # snapshot everything that's been added so far — those artists stay FIXED
    # (stars + mount).  Everything drawn after this block is part of the
    # telescope assembly and will rotate with the "tilt" slider.
    _fixed_ids = {id(a) for lst in (ax.lines, ax.patches, ax.texts,
                                     ax.collections) for a in lst}

    draw_tube()
    draw_primary()
    draw_cage_and_corrector()


    # ==========================================================================
    # Ray bundles — 5 field angles fanning across the FOV
    # ==========================================================================

    ray_lines = {}    # field_angle_deg → list of Line2D
    focal_spots = {}  # field_angle_deg → Line2D marker

    incoming_rays = []     # l_in segments (sky → primary); hidden when zoomed
    reflected_rays = []    # l_mid segments (primary → corrector); hidden when zoomed


    def draw_rays():
        for ang, color in zip(FIELD_ANGLES, RAY_COLORS):
            rays, spot = trace_bundle(ang, n_rays=N_RAYS_PER_BUNDLE)
            lines = []
            for xs, ys in rays:
                l_in, = ax.plot(xs[:2], ys[:2], color=color, linewidth=1.0,
                                alpha=0.55, solid_capstyle='round', zorder=6)
                incoming_rays.append(l_in)
                l_mid, = ax.plot(xs[1:3], ys[1:3], color=color, linewidth=1.4,
                                 alpha=0.85, solid_capstyle='round', zorder=6)
                reflected_rays.append(l_mid)
                l_out, = ax.plot(xs[2:4], ys[2:4], color=color, linewidth=1.2,
                                 alpha=0.75, linestyle=(0, (3, 2)),
                                 solid_capstyle='round', zorder=6)
                lines.extend([l_in, l_mid, l_out])
            # Focal spot marker
            spot_marker, = ax.plot([spot[0]], [spot[1]], 'o', color=color,
                                   markersize=5, zorder=7, alpha=0.95,
                                   markeredgecolor='white', markeredgewidth=0.6)
            ray_lines[ang] = lines
            focal_spots[ang] = spot_marker

    draw_rays()


    # ==========================================================================
    # Telescope tilt — rotate the whole tube/primary/corrector/FP/rays
    # assembly around the declination axis (mount trunnions).  Default
    # position is Kitt Peak latitude (32°); slider below lets the user
    # swing the telescope like a real equatorial mount.
    # ==========================================================================

    TUBE_PIVOT = (0.0, (TUBE_BOT_Y + TUBE_TOP_Y) / 2 - 1.5)

    tube_artists = [a for lst in (ax.lines, ax.patches, ax.texts,
                                   ax.collections)
                    for a in lst if id(a) not in _fixed_ids]

    # One shared rotation transform — artists hold a reference to it and see
    # updates automatically when we mutate its matrix in place.  This is
    # much cheaper than rebuilding a transform per slider tick.
    _tube_rotation = Affine2D()
    _tube_transform = _tube_rotation + ax.transData
    for _a in tube_artists:
        _a.set_transform(_tube_transform)


    def apply_tube_tilt(deg):
        _tube_rotation.clear()
        _tube_rotation.rotate_deg_around(TUBE_PIVOT[0], TUBE_PIVOT[1], deg)


    # ==========================================================================
    # Legend (right side)
    # ==========================================================================

    legend_handles = []
    for ang, color in zip(FIELD_ANGLES, RAY_COLORS):
        legend_handles.append(Line2D([0], [0], color=color, linewidth=2.4,
                                     label=f'{ang:+.1f}°  field angle'))
    legend = ax.legend(handles=legend_handles,
                       loc='upper left', bbox_to_anchor=(1.02, 1.0),
                       frameon=True, fancybox=True, borderpad=0.8,
                       title='Ray bundles', title_fontsize=10, fontsize=9,
                       facecolor=PANEL_BG, edgecolor=SPINE, labelcolor=TEXT)
    legend.get_title().set_color(TEXT)
    legend_pick_to_angle = {}
    for ll, txt, ang in zip(legend.get_lines(), legend.get_texts(), FIELD_ANGLES):
        ll.set_linewidth(3.0)
        # Only the TEXT label is pickable — swatch picker was causing Button
        # widget clicks to misfire on the macOS backend.
        txt.set_picker(True)
        legend_pick_to_angle[txt] = ang

    # Annotation box: geometry summary
    info_text = (
        'Mayall 4 m  ·  f/2.7 paraboloid\n'
        'DESI corrector  ·  6 lenses\n'
        'Delivered f/3.9  ·  3.2° FOV\n'
        'Focal plane  ·  0.8 m ⌀ curved\n'
        '5 020 robotic fiber positioners\n'
        'Equatorial horseshoe-yoke mount'
    )
    ax.text(1.02, 0.50, info_text, transform=ax.transAxes,
            fontsize=9, color=TEXT_DIM, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.6', facecolor=PANEL_BG,
                      edgecolor=SPINE, linewidth=0.8))

    # Component callouts
    ax.annotate('Primary mirror\n(4.0 m, f/2.7)',
                xy=(0, -0.05), xytext=(-3.3, -1.6),
                fontsize=8.5, color=C_MIRROR, ha='left',
                arrowprops=dict(arrowstyle='-', color=C_MIRROR, linewidth=0.8))
    ax.annotate('DESI corrector\n(6 lenses, C1→C4)',
                xy=(CORR_R_BOT + 0.05, (CORR_TOP_Y + CORR_BOT_Y) / 2),
                xytext=(3.3, 9.6),
                fontsize=9, color=C_LENS, ha='left', fontweight='medium',
                arrowprops=dict(arrowstyle='-', color=C_LENS, linewidth=0.8))
    ax.annotate('Focal plane\n(fibers hang down)',
                xy=(FOCAL_R + 0.02, FOCAL_Y + FOCAL_SAG),
                xytext=(2.7, 12.4),
                fontsize=9, color=C_FOCAL, ha='left', fontweight='medium',
                arrowprops=dict(arrowstyle='-', color=C_FOCAL, linewidth=0.8))


    # ==========================================================================
    # Axes limits + style
    # ==========================================================================

    DEFAULT_XLIM = (-11, 11)
    DEFAULT_YLIM = (-4.8, TUBE_TOP_Y + 8.5)
    ZOOM_XLIM    = (-1.3, 1.3)
    ZOOM_YLIM    = (CORR_BOT_Y - 0.3, FOCAL_Y + 0.4)

    ax.set_xlim(*DEFAULT_XLIM)
    ax.set_ylim(*DEFAULT_YLIM)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


    # ==========================================================================
    # Controls
    # ==========================================================================

    def mkax(x, y, w, h, fc=PANEL_BG):
        a = fig.add_axes([x, y, w, h]); a.set_facecolor(fc)
        return a

    def mkbtn(x, y, w, h, label):
        a = mkax(x, y, w, h, '#2a2f4a')
        b = Button(a, label, color='#2a2f4a', hovercolor='#3d4470')
        b.label.set_color(TEXT); b.label.set_fontsize(9)
        return b

    btn_zoom   = mkbtn(0.06, 0.12, 0.14, 0.035, 'Zoom to focal plane')
    btn_reset  = mkbtn(0.21, 0.12, 0.10, 0.035, 'Reset view')
    btn_save   = mkbtn(0.32, 0.12, 0.10, 0.035, 'Save PNG')

    # Slider: telescope tilt (rotates tube around declination axis)
    tilt_ax = mkax(0.22, 0.17, 0.25, 0.022, '#10121f')
    tilt = Slider(tilt_ax, 'Telescope tilt (°)',
                  -60.0, 60.0, valinit=LATITUDE_DEG, valstep=1,
                  color=C_MIRROR, track_color='#20263e')
    tilt.label.set_color(TEXT); tilt.label.set_fontsize(9)
    tilt.valtext.set_color(TEXT); tilt.valtext.set_fontsize(9)

    # Status line
    status_ax = mkax(0.58, 0.12, 0.36, 0.035, '#10121f')
    status_ax.set_xticks([]); status_ax.set_yticks([])
    status = status_ax.text(0.02, 0.5, '', color=TEXT_DIM, fontsize=9,
                            va='center', transform=status_ax.transAxes)

    def info(msg):
        status.set_text(msg); fig.canvas.draw_idle()


    # ==========================================================================
    # Interactivity
    # ==========================================================================

    state = {'isolated': None, 'zoomed': False}

    def set_bundle_alpha(ang, dim=False):
        lines = ray_lines[ang]
        spot = focal_spots[ang]
        if dim:
            for ln in lines: ln.set_alpha(0.10)
            spot.set_alpha(0.25)
        else:
            # restore default alpha (remember types: incoming/mid/out)
            for i, ln in enumerate(lines):
                which = i % 3
                ln.set_alpha([0.55, 0.85, 0.75][which])
            spot.set_alpha(0.95)

    def isolate(ang):
        for a in FIELD_ANGLES:
            if np.isclose(a, ang):
                for i, ln in enumerate(ray_lines[a]):
                    which = i % 3
                    ln.set_alpha([0.85, 1.0, 0.95][which])
                    ln.set_linewidth([1.3, 1.8, 1.5][which])
                focal_spots[a].set_alpha(1.0)
                focal_spots[a].set_markersize(8)
            else:
                set_bundle_alpha(a, dim=True)

    def reset_isolation():
        for a in FIELD_ANGLES:
            set_bundle_alpha(a, dim=False)
            for i, ln in enumerate(ray_lines[a]):
                which = i % 3
                ln.set_linewidth([1.0, 1.4, 1.2][which])
            focal_spots[a].set_markersize(5)

    def _zoom_limits_for_tilt(deg):
        """Rotated bounding box around the corrector+focal-plane region so
           the zoom follows the telescope's current tilt."""
        x0, x1 = ZOOM_XLIM
        y0, y1 = ZOOM_YLIM
        corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        theta = np.deg2rad(deg)
        cx, cy = TUBE_PIVOT
        c, s = np.cos(theta), np.sin(theta)
        xs, ys = [], []
        for x, y in corners:
            dx, dy = x - cx, y - cy
            xs.append(cx + dx * c - dy * s)
            ys.append(cy + dx * s + dy * c)
        margin = 0.2
        return (min(xs) - margin, max(xs) + margin,
                min(ys) - margin, max(ys) + margin)


    def on_zoom(_):
        if state['zoomed']:
            ax.set_xlim(*DEFAULT_XLIM); ax.set_ylim(*DEFAULT_YLIM)
            state['zoomed'] = False
            btn_zoom.label.set_text('Zoom to focal plane')
            for a in incoming_rays + reflected_rays:
                a.set_visible(True)
            info('full telescope view')
        else:
            xl, xr, yl, yu = _zoom_limits_for_tilt(tilt.val)
            ax.set_xlim(xl, xr); ax.set_ylim(yl, yu)
            state['zoomed'] = True
            btn_zoom.label.set_text('Zoom out')
            # hide incoming + primary-reflected segments — they originate
            # outside the zoom and don't actually land on the FP
            for a in incoming_rays + reflected_rays:
                a.set_visible(False)
            info('zoomed to focal plane — see how each field angle lands at a '
                 'different x on the curved surface')
        fig.canvas.draw_idle()

    def on_reset(_):
        ax.set_xlim(*DEFAULT_XLIM); ax.set_ylim(*DEFAULT_YLIM)
        state['zoomed'] = False
        btn_zoom.label.set_text('Zoom to focal plane')
        for a in incoming_rays + reflected_rays:
            a.set_visible(True)
        reset_isolation()
        state['isolated'] = None
        info('view reset')
        fig.canvas.draw_idle()

    def on_save(_):
        out_dir = os.path.dirname(os.path.abspath(__file__)) \
                  if '__file__' in globals() else '.'
        out = os.path.join(
            out_dir,
            f'mayall_desi_raytrace_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png')
        fig.savefig(out, dpi=160, facecolor=BG, bbox_inches='tight')
        info(f'saved → {os.path.basename(out)}')
        print(f'saved → {out}')

    def _toggle_angle(ang):
        if state['isolated'] is not None and np.isclose(state['isolated'], ang):
            reset_isolation()
            state['isolated'] = None
            info('cleared isolation')
        else:
            isolate(ang)
            state['isolated'] = ang
            info(f'isolated {ang:+.1f}° field angle')
        fig.canvas.draw_idle()


    def on_pick(event):
        """Pick on legend swatches OR text labels — whole row clickable."""
        ang = legend_pick_to_angle.get(event.artist)
        if ang is not None:
            _toggle_angle(ang)

    btn_zoom.on_clicked(on_zoom)
    btn_reset.on_clicked(on_reset)
    btn_save.on_clicked(on_save)


    def on_tilt(v):
        apply_tube_tilt(float(v))
        if state['zoomed']:
            xl, xr, yl, yu = _zoom_limits_for_tilt(float(v))
            ax.set_xlim(xl, xr); ax.set_ylim(yl, yu)
        info(f'telescope tilt: {v:+.0f}° from vertical')
        # NOTE: info() already schedules a draw_idle; don't double-schedule
    tilt.on_changed(on_tilt)

    # Apply the default Kitt Peak latitude tilt once at startup
    apply_tube_tilt(LATITUDE_DEG)
    fig.canvas.mpl_connect('pick_event', on_pick)


    # Keep widget refs alive
    fig._keep_widgets = (btn_zoom, btn_reset, btn_save, tilt)

    info('click a legend entry to isolate · drag slider for a probe ray · '
         'zoom, toggle mount, save')

    plt.show()


if __name__ == "__main__":
    main()
