import numpy as np
import matplotlib.pyplot as plt

# ============================================
# LOAD DATA
# ============================================
# u_visual is n1 x n1 (cell-centered)
u = np.loadtxt("u_visual.txt")

# Grid size
Ny, Nx = u.shape
assert Ny == Nx, "Grid must be square (n1 x n1)"

n1 = Nx

# ============================================
# DOMAIN
# ============================================
Lx = 1.0
Ly = 1.0

dx = Lx / n1
dy = Ly / n1

# Cell-centered y coordinates
y = np.linspace(dy/2, Ly - dy/2, n1)

# ============================================
# CENTERLINE INDEX (x = 0.5)
# ============================================
center_i = n1 // 2   # integer index of center column

# Extract u velocity along centerline
u_centerline = u[:, center_i]

# ============================================
# SAVE CENTERLINE DATA
# ============================================
np.savetxt(
    "u_centerline_vs_y.txt",
    np.column_stack((y, u_centerline)),
    header="y   u_centerline",
    comments=""
)

# ============================================
# PLOT: u vs y
# ============================================
plt.figure(figsize=(5, 7))
plt.plot(u_centerline, y, marker="o")
plt.xlabel("u velocity")
plt.ylabel("y")
plt.title("Centerline u velocity vs y")
plt.grid(True)
plt.tight_layout()
plt.show()
