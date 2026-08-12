"""
BEM-Optimized Small Wind Turbine Blade Generator
=================================================
Computes the Betz/Schmitz-optimum chord and twist distribution for a
horizontal-axis wind turbine blade, builds a lofted 3D surface from
NACA 4-digit airfoil sections, and writes an ASCII STL mesh ready for
slicing/3D printing (or CNC/laser-cut ribs).

Design theory (why this maximizes energy conversion):
  For a given tip-speed ratio (TSR) and number of blades, there is a
  unique chord(r) and twist(r) that makes the local blade element
  operate at its airfoil's design angle of attack (max Cl/Cd) at every
  span station simultaneously. Deviating from this (e.g. constant
  chord, no twist) sharply drops the achievable power coefficient Cp.

  Flow angle:      phi(r) = (2/3) * atan(1 / lambda_r)
  Chord (Betz):     c(r)  = 8*pi*r*(1 - cos(phi)) / (B * Cl_design)
  Twist:          theta(r) = phi(r) - alpha_design
  where lambda_r = TSR * r / R
"""

import numpy as np

# ---------------------------------------------------------------
# 1. DESIGN PARAMETERS  (tune these for your build)
# ---------------------------------------------------------------
R = 0.25            # rotor radius, m (0.25 m -> 50 cm diameter rotor)
R_HUB = 0.03         # hub radius, m (blade root starts here)
B = 3                # number of blades (3 = good balance of Cp & smooth torque)
TSR_DESIGN = 6.0     # design tip-speed ratio (5-7 typical for 3-blade lift rotors)
CL_DESIGN = 0.9      # design lift coefficient of chosen airfoil section
ALPHA_DESIGN_DEG = 6.0   # angle of attack at max Cl/Cd for that airfoil (deg)
NACA = "4412"        # airfoil: cambered, forgiving at low Reynolds number
N_STATIONS = 14       # number of spanwise sections to loft
N_AIRFOIL_PTS = 40    # points per airfoil section (upper+lower)
TIP_CHORD_MIN = 0.008 # m, floor so tip doesn't taper to zero (manufacturability)

# ---------------------------------------------------------------
# 2. NACA 4-DIGIT AIRFOIL COORDINATE GENERATOR
# ---------------------------------------------------------------
def naca4_coords(code, n_pts=40, closed_te=True):
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:4]) / 100.0

    beta = np.linspace(0, np.pi, n_pts)
    x = (1 - np.cos(beta)) / 2  # cosine spacing: denser near LE/TE

    yt = 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                  + 0.2843*x**3 - (0.1036 if closed_te else 0.1015)*x**4)

    yc = np.where(x < p, m/p**2 * (2*p*x - x**2) if p > 0 else 0,
                  m/(1-p)**2 * ((1-2*p) + 2*p*x - x**2) if p > 0 else 0)
    dyc = np.where(x < p, 2*m/p**2 * (p - x) if p > 0 else 0,
                   2*m/(1-p)**2 * (p - x) if p > 0 else 0)
    theta = np.arctan(dyc)

    xu = x - yt*np.sin(theta); yu = yc + yt*np.cos(theta)
    xl = x + yt*np.sin(theta); yl = yc - yt*np.cos(theta)

    # single loop: TE(upper) -> LE -> TE(lower), i.e. closed polygon
    xs = np.concatenate([xu[::-1], xl[1:]])
    ys = np.concatenate([yu[::-1], yl[1:]])
    return xs, ys  # in chord fraction, LE at x=0, TE at x=1, quarter-chord ref x=0.25

# ---------------------------------------------------------------
# 3. BEM OPTIMUM CHORD & TWIST DISTRIBUTION
# ---------------------------------------------------------------
def bem_optimum(R, R_hub, B, tsr, cl, alpha_deg, n_stations):
    r = np.linspace(R_hub, R, n_stations)
    lam_r = tsr * r / R
    phi = (2.0/3.0) * np.arctan(1.0 / lam_r)          # flow angle, rad
    chord = 8 * np.pi * r * (1 - np.cos(phi)) / (B * cl)
    chord = np.maximum(chord, TIP_CHORD_MIN)
    twist = phi - np.radians(alpha_deg)                 # rad, blade setting angle
    return r, chord, twist, np.degrees(phi)

r_stations, chord, twist, phi_deg = bem_optimum(
    R, R_HUB, B, TSR_DESIGN, CL_DESIGN, ALPHA_DESIGN_DEG, N_STATIONS)

print(f"{'r (mm)':>8} {'chord (mm)':>11} {'twist (deg)':>12} {'phi (deg)':>10}")
for ri, ci, ti, pi_ in zip(r_stations, chord, twist, phi_deg):
    print(f"{ri*1000:8.1f} {ci*1000:11.2f} {np.degrees(ti):12.2f} {pi_:10.2f}")

# save table for reference / CNC rib cutting
np.savetxt("blade_schedule.csv",
           np.column_stack([r_stations, chord, np.degrees(twist), phi_deg]),
           header="r_m,chord_m,twist_deg,phi_deg", delimiter=",", comments="")

# ---------------------------------------------------------------
# 4. LOFT SECTIONS INTO A 3D BLADE SURFACE, EXPORT STL
# ---------------------------------------------------------------
def build_blade_mesh(r_stations, chord, twist, naca_code, n_pts):
    xs_af, ys_af = naca4_coords(naca_code, n_pts)  # chord-normalized 2D profile
    sections = []
    for ri, ci, ti in zip(r_stations, chord, twist):
        # place quarter-chord (x=0.25) on the blade pitch axis, scale by chord
        x = (xs_af - 0.25) * ci
        y = ys_af * ci
        ct, st = np.cos(ti), np.sin(ti)
        xr = x*ct - y*st      # rotate section by twist about pitch axis (in-plane)
        yr = x*st + y*ct
        pts3d = np.column_stack([np.full_like(xr, ri), xr, yr])  # (spanwise r, chordwise x, thickness y)
        sections.append(pts3d)
    return np.array(sections)  # shape (n_stations, n_pts_loop, 3)

sections = build_blade_mesh(r_stations, chord, twist, NACA, N_AIRFOIL_PTS)

def write_stl(sections, filename, root_cap=True, tip_cap=True):
    n_stations, n_loop, _ = sections.shape
    tris = []

    def tri(a, b, c):
        tris.append((a, b, c))

    # side faces: connect consecutive stations around the loop
    for s in range(n_stations - 1):
        for k in range(n_loop - 1):
            p1, p2 = sections[s, k], sections[s, k+1]
            p3, p4 = sections[s+1, k], sections[s+1, k+1]
            tri(p1, p2, p3)
            tri(p2, p4, p3)

    # root cap (fan triangulation from centroid)
    if root_cap:
        c = sections[0].mean(axis=0)
        for k in range(n_loop - 1):
            tri(c, sections[0, k+1], sections[0, k])

    # tip cap
    if tip_cap:
        c = sections[-1].mean(axis=0)
        for k in range(n_loop - 1):
            tri(c, sections[-1, k], sections[-1, k+1])

    with open(filename, "w") as f:
        f.write("solid blade\n")
        for a, b, c in tris:
            n = np.cross(b - a, c - a)
            norm = np.linalg.norm(n)
            n = n / norm if norm > 1e-12 else np.array([0, 0, 1])
            f.write(f" facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n")
            f.write("  outer loop\n")
            for v in (a, b, c):
                f.write(f"   vertex {v[0]:.6e} {v[1]:.6e} {v[2]:.6e}\n")
            f.write("  endloop\n endfacet\n")
        f.write("endsolid blade\n")

write_stl(sections, "blade.stl")
print("\nWrote blade.stl and blade_schedule.csv")
print(f"Rotor swept area: {np.pi*R**2:.4f} m^2  |  Rotor diameter: {2*R*100:.0f} cm")
