import numpy as np
import matplotlib.pyplot as plt

# =================================================
# LOAD DATA
# =================================================
u = np.loadtxt("u_visual.txt")
v = np.loadtxt("v_visual.txt")
p = np.loadtxt("p_1.txt")

# Sanity check
assert u.shape == v.shape == p.shape, "u, v, p must have same shape"

Ny, Nx = u.shape
print(f"Grid detected: {Nx} x {Ny}")

# =================================================
# DOMAIN SIZE (unit square cavity)
# =================================================
Lx = 1.0
Ly = 1.0

dx = Lx / Nx
dy = Ly / Ny

# Cell-centered coordinates
x = np.linspace(dx/2, Lx - dx/2, Nx)
y = np.linspace(dy/2, Ly - dy/2, Ny)

X, Y = np.meshgrid(x, y)

# =================================================
# 1️⃣ U VELOCITY — COLORMAP
# =================================================
plt.figure(figsize=(7, 6))
plt.contourf(X, Y, u, levels=50, cmap="coolwarm")
plt.colorbar(label="u velocity")
plt.xlabel("x")
plt.ylabel("y")
plt.title("u velocity (colormap)")
plt.axis("equal")
plt.tight_layout()
plt.show()

# =================================================
# 2️⃣ U VELOCITY — ARROWS
# =================================================
plt.figure(figsize=(7, 6))
plt.quiver(X, Y, u, np.zeros_like(u), scale=50)
plt.xlabel("x")
plt.ylabel("y")
plt.title("u velocity (arrows)")
plt.axis("equal")
plt.tight_layout()
plt.show()

# =================================================
# 3️⃣ V VELOCITY — COLORMAP
# =================================================
plt.figure(figsize=(7, 6))
plt.contourf(X, Y, v, levels=50, cmap="coolwarm")
plt.colorbar(label="v velocity")
plt.xlabel("x")
plt.ylabel("y")
plt.title("v velocity (colormap)")
plt.axis("equal")
plt.tight_layout()
plt.show()

# =================================================
# 4️⃣ V VELOCITY — ARROWS
# =================================================
plt.figure(figsize=(7, 6))
plt.quiver(X, Y, np.zeros_like(v), v, scale=50)
plt.xlabel("x")
plt.ylabel("y")
plt.title("v velocity (arrows)")
plt.axis("equal")
plt.tight_layout()
plt.show()

# =================================================
# 5️⃣ PRESSURE — COLORMAP (scaled)
# =================================================
pmin, pmax = p.min(), p.max()
eps = 0.1 * (pmax - pmin)

plt.figure(figsize=(7, 6))
plt.contourf(
    X, Y, p,
    levels=50,
    cmap="viridis",
    vmin=pmin + eps,
    vmax=pmax - eps
)
plt.colorbar(label="pressure")
plt.xlabel("x")
plt.ylabel("y")
plt.title("pressure (enhanced contrast)")
plt.axis("equal")
plt.tight_layout()
plt.show()