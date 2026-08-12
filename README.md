# BEM_Blade_Generator

A Python tool that computes aerodynamically optimized wind turbine blade geometry using Blade Element Momentum (BEM) theory and exports a ready-to-print 3D mesh.

Instead of guessing at a chord/twist profile, it derives the Betz-optimum shape for a given rotor radius, blade count, and design tip-speed ratio — then lofts that geometry through NACA 4-digit airfoil sections into a triangulated STL.

Computes the optimum chord and twist at each spanwise station using the Betz/Schmitz BEM equations, so every section of the blade operates near its airfoil's best lift-to-drag angle.

Generates NACA 4-digit airfoil coordinates for the blade cross-section (default: 4412).

Lofts the sections — scaled by chord, rotated by twist about the quarter-chord axis — into a 3D surface.

Exports the result as an ASCII STL (no external mesh library) and a CSV blade schedule.

Requirements:

Python 3.8+
NumPy
bash
pip install numpy
Usage
bash
python blade_generator.py

This prints the blade schedule to the console and writes two files to the working directory:

File	Description:

blade.stl	3D-printable blade mesh

blade_schedule.csv	Per-station radius, chord, twist, and flow angle — useful for CNC rib-cutting as a print-free alternative
Configuration

All design parameters live in the constants block at the top of blade_generator.py:

Parameter	Default	Description
R	0.25 m	Rotor radius
R_HUB	0.03 m	Hub radius (blade root start)
B	3	Number of blades
TSR_DESIGN	6.0	Design tip-speed ratio
CL_DESIGN	0.9	Design lift coefficient
ALPHA_DESIGN_DEG	6.0	Angle of attack at max Cl/Cd
NACA	"4412"	Airfoil (any valid NACA 4-digit code)
N_STATIONS	14	Number of spanwise sections to loft
N_AIRFOIL_PTS	40	Points per airfoil section
TIP_CHORD_MIN	0.008 m	Minimum tip chord for manufacturability

Edit these and rerun — the chord/twist schedule, mesh, and CSV all update accordingly.

Theory:

For a given tip-speed ratio and blade count, there's a unique chord(r) and twist(r) that puts every blade element at its airfoil's design angle of attack simultaneously:

lambda_r  = TSR * r / R

phi(r)    = (2/3) * atan(1 / lambda_r)          # flow angle

chord(r)  = 8*pi*r*(1 - cos(phi)) / (B * Cl)     # Betz-optimum chord

twist(r)  = phi(r) - alpha_design                # blade setting angle

Deviating from this (constant chord, no twist) sharply reduces the achievable power coefficient. This script implements that derivation directly, rather than approximating it.


Root: wide chord, steep twist. Tip: slim, nearly flat — the classic BEM-optimum profile.

Next steps:

Import blade.stl into a slicer for 3D printing, or use blade_schedule.csv to cut ribs by hand/CNC.
Combine with hub_nacelle.scad for a full parametric turbine assembly.

Swap in Reynolds-number-appropriate Cl/alpha values (e.g. from XFOIL) for a specific airfoil and blade size.

Limitations:

Assumes ideal inviscid BEM (no drag, no tip/hub losses) — real Cp will be somewhat lower than the theoretical optimum.

No Reynolds-number correction; CL_DESIGN/ALPHA_DESIGN_DEG should be adjusted for very small chords (low Re).

NACA 4-digit airfoils only; no support for custom .dat coordinate files (yet).
