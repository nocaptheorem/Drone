# Procedural Physics Lab: 6-DoF Flight Dynamics Engine

[![Godot Engine](https://img.shields.io/badge/Godot-v4.x--.NET-blue?logo=godotengine&logoColor=white)](https://godotengine.org)
[![.NET](https://img.shields.io/badge/.NET-v8.0-purple?logo=dotnet&logoColor=white)](https://dotnet.microsoft.com/)
[![C#](https://img.shields.io/badge/Language-C%23_12-green?logo=csharp&logoColor=white)](https://docs.microsoft.com/en-us/dotnet/csharp/)
[![License](https://img.shields.io/badge/License-Apache2-yellow.svg)](LICENSE)
[![Channel](https://img.shields.io/badge/YouTube-NoCapTheorem-red?logo=youtube&logoColor=white)](https://www.youtube.com/@NoCapTheorem)

> **A high-fidelity, first-principles 6-DoF drone simulation built in C# and Godot 4.**
> Modeled with explicit 4th-order Runge-Kutta tensor integration, dynamic center-of-mass shifts via the Parallel Axis Theorem, sub-step cable damping, fault-tolerant control allocation, and dynamic shader-driven environmental interactions.

---

## Quick Start & Installation

### Prerequisites
* **Godot Engine v4.x** (specifically the **.NET edition**)
* **.NET 8.0 SDK** (or higher)

### Build & Run
Clone the repository and compile/run the C# solution:

```bash
# Build C# solution and execute in Godot
dotnet build && godot --headless --build-solutions --verbose Main.tscn
```

---

## License & Channel

Distributed under the **Apache2 License**. See `LICENSE` for more information.

* **No-Cap Theorem Platform:** For deep-dive architectural breakdowns, formal specs, and video walkthroughs on systems engineering and simulation physics, check out **No-Cap Theorem**.

---

---

# Technical Specification & Documentation

## 1. What Is It Simulating?

This codebase simulates a **physics-driven, multi-rotor aerial vehicle (drone) operating within a dynamic, non-linear environment**. Rather than relying on simplified arcade movement or rigid-body presets, it models physical flight dynamics, aerodynamic forces, structural damage, dynamic center-of-mass shifts, and environment interactions.

### Key Systems Simulated

* **Full 6-DoF Rigid-Body Mechanics:** High-order integration of spatial position, linear velocity, orientation (quaternions), and angular velocity using full 3D inertia tensors.
* **Aerodynamics & Gyroscopic Dynamics:** Blade spin inertia, angular momentum ($\mathbf{L} = \mathbf{I}\boldsymbol{\omega}$), gyroscopic cross-coupling torques ($\boldsymbol{\omega} \times \mathbf{L}$), directional rotor drag, vortex ring states (VRS), and ground effect dynamics.
* **Structural Breakage & Real-time CoM Recalculation:** Dynamic recalculation of mass, center of mass (CoM), and principal axes of inertia when rotors detach during extreme force impacts.
* **Cascaded Flight Control System:** Multi-loop attitude stabilization, altitude PID control, feed-forward tension compensation, and closed-form dynamic thrust allocation for under-actuated (3-rotor or asymmetrical) configurations.
* **Dynamic Winch & Electromagnet Load System:** Critically damped spring-damper cable physics carrying external, variable-mass rigid bodies.
* **Environmental Interaction Fields:** Spatial turbulence volumes, rotor downwash drag cones acting on surrounding objects, dynamic terrain cratering via deformation shaders, seamless state-flipping spatial portals, and an orbiting gravity anomaly (pulsar/vortex).

---

## 2. Technical Breakdown

### Component 1: Flight State & 3D Tensor RK4 Integrator

Standard Euler integration introduces rapid numerical drift and energy explosions in high-speed rotational mechanics. To solve this, the integration loop uses an explicit **4th-Order Runge-Kutta (RK4)** numerical integrator evaluated across 4 sub-steps per frame ($120\text{ Hz} \times 4 = 480\text{ Hz}$ effective integration frequency).

#### State Integration Mechanics

* **Mass & Inertia Tensors:** Inertia is represented as a $3 \times 3$ matrix ($\mathbf{I}$) and its inverse ($\mathbf{I}^{-1}$).
* **Gyroscopic Precession:** At each evaluation step, the angular momentum $\mathbf{L} = \mathbf{I}\boldsymbol{\omega}$ is computed, yielding the gyroscopic cross-coupling torque:

$$\boldsymbol{\tau}_{\text{gyro}} = \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega})$$

* **Rotational Acceleration:** Angular acceleration is evaluated in local space as:

$$\boldsymbol{\alpha} = \mathbf{I}^{-1}(\boldsymbol{\tau}_{\text{local}} - \boldsymbol{\tau}_{\text{gyro}})$$

* **Quaternion Integration:** Orientation derivatives are computed via quaternion rates:

$$\dot{\mathbf{q}} = \frac{1}{2} \mathbf{q} \otimes \boldsymbol{\omega}$$

#### Dynamic CoM Shift & Parallel Axis Theorem

When a rotor breaks off at local position `$ \mathbf{r}_{\text{part}} $`, the mass drops to `$ M_{\text{new}} = M - m_{\text{part}} $`. The Center of Mass shifts by:

$$\Delta\mathbf{r}_{\text{CoM}} = \frac{-m_{\text{part}} \mathbf{r}_{\text{part}}}{M_{\text{new}}}$$

The code updates the primary inertia tensor by subtracting the removed component's inertia tensor via discrete cross-product matrix terms and applying the generalized Parallel Axis Theorem transformation:

$$\mathbf{I}_{\text{shift}} = m \left( (\mathbf{r} \cdot \mathbf{r})\mathbf{E} - \mathbf{r} \otimes \mathbf{r} \right)$$

Crucially, attachments like the winch hook offset ($\mathbf{r}_{\text{hook}}$) and remaining motor positions ($\mathbf{r}_{\text{rotor}, i}$) are updated relative to the shifted origin to preserve torque balance correctness ($\boldsymbol{\tau} = \mathbf{r} \times \mathbf{F}$).

---

### Component 2: Flight Controllers & Mixer

#### Cascaded Attitude & Altitude Controller

* **Outer Attitude Loop:** Computes orientation error quaternions $\mathbf{q}_{\text{err}} = \mathbf{q}_{\text{current}}^{-1} \otimes \mathbf{q}_{\text{target}}$, converts to axis-angle representation, and yields a target angular rate vector scaled by $K_{p,\text{outer}}$.
* **Inner Rate Loop:** A proportional-derivative rate controller maps target rate errors to commanded body torques:

$$\boldsymbol{\tau}_{\text{cmd}} = \mathbf{K}_{p,\text{rate}} (\boldsymbol{\omega}_{\text{target}} - \boldsymbol{\omega}) - \mathbf{K}_{d,\text{rate}} \boldsymbol{\omega}$$

* **Altitude PID & Anti-Windup:** Maintains vertical position. If motor thrust saturates ($T_i \ge T_{\text{max}}$) while an altitude deficit exists, the integral accumulator is frozen (`freezeIntegral = true`) to eliminate integrator windup.

#### Dynamic Allocation Matrix (3-Motor Fault-Tolerant Control)

When all 4 motors are active, a standard cross-mix matrix controls pitch, roll, and yaw. However, when a motor fails ($N=3$), standard mixing breaks down. The class construct solves a $3 \times 3$ linear allocation equation using basis matrix inversion:

$$\begin{bmatrix} T_{\text{total}} \\ \tau_x \\ \tau_z \end{bmatrix} = \begin{bmatrix} 1 & 1 & 1 \\ -z_1 & -z_2 & -z_3 \\ x_1 & x_2 & x_3 \end{bmatrix} \begin{bmatrix} T_1 \\ T_2 \\ T_3 \end{bmatrix}$$

By calculating $\mathbf{B}^{-1}$, the controller allocates asymmetric individual motor thrusts to maintain flat hover and directional roll/pitch authority, sacrificing yaw control to prioritize basic structural flight.

---

### Component 3: Cable Physics & Dynamic Damping

The winch system connects a secondary `RigidBody3D` (magnet body) to the integrated drone through sub-step cable constraints.

* **Dynamic Critical Damping Formulation:** Static damping factors cause severe numerical oscillations under variable mass loads. The system continuously calculates the effective reduced mass ($m_{\text{eff}}$) of the combined system:

$$m_{\text{eff}} = \frac{m_{\text{drone}} \cdot m_{\text{load}}}{m_{\text{drone}} + m_{\text{load}}}$$

* **Critical Damping Coefficient:** Computed continuously to adapt to attached loads:

$$c_{\text{crit}} = 2 \sqrt{k \cdot m_{\text{eff}}}, \quad c = c_{\text{crit}} \cdot \zeta$$

* **Feed-Forward Tension Injection:** The downward tension force $F_{\text{cable, Y}}$ generated by the payload is extracted within sub-steps and added directly to the motor command thrust $T_{\text{hover}}$ prior to feedback execution, eliminating altitude sag during heavy load pickups.

---

### Component 4: Environmental Interaction Dynamics

```
               [ Gravity Anomaly / Vortex ]
                            |
                     (Swirling Pull)
                            v
[ FastNoiseLite ] ---> [ Turbulence Volume ] ---> [ Rigid Body Integration ]
                            |                                 |
                            v                                 v
                     (Downwash Cone)                 (Terrain Deformation)
                            |                                 |
                            v                                 v
                   [ Dynamic Objects ]               [ Crater Shader Array ]

```

* **Ground Effect & VRS:** Thrust multiplication increases exponentially near ground surfaces ($h < 10\text{m}$) using $T_{\text{effective}} = T \cdot (1 + 0.5 e^{-2h})$. Sinking into downwash streams ($v_y < -2\text{ m/s}$) triggers Vortex Ring State (VRS), inducing random destabilizing torque spikes and thrust loss.
* **Rotor Downwash Fields:** Computes induced velocity $v_{\text{induced}} = \sqrt{\frac{T}{2 \rho A}}$ and applies a Gaussian radial drag distribution to downward physics objects.
* **Orbital Gravity Anomaly:** Generates an inward radial pull balanced with a $40\%$ orthogonal swirl force vector. It utilizes an $\epsilon$-softened denominator ($r^2 + \epsilon^2$) to flatten gravity spikes near the center, creating a stable vortex trajectory.

---

## 3. How to Control It and Push It to Its Limits

### Controls Reference Table

| Category | Input / Key | Action |
| --- | --- | --- |
| **Attitude** | `W` / `S` | Pitch Down / Pitch Up |
|  | `A` / `D` | Roll Left / Roll Right |
|  | `Q` / `E` | Yaw Rate Left / Yaw Rate Right |
| **Throttle** | `SPACE` / `SHIFT` | Increase Target Altitude / Decrease Target Altitude |
| **Motor Faults** | `1`, `2`, `3`, `4` | Toggle FL, FR, BL, BR Motor Power/Failure States |
| **Thrust Inversion** | `F` | Invert Thrust Vector (Inverted Flight / Downward Force) |
| **Winch Control** | `Z` / `C` | Reel Cable Out / Reel Cable In |
|  | `X` | Toggle Electromagnet On/Off |
| **System** | `TAB` | Toggle Camera Yaw Follow Mode |
|  | `ESC` | Release/Capture Mouse Cursor |

---

### Extreme Stress Tests & Edge Cases

#### 1. Asymmetric Motor Loss & Survival Flight

* **Action:** Press `1` to kill Motor 1 (Front-Left) during hover.
* **What Happens:** The system drops yaw control, recalculates the $3 \times 3$ allocation matrix using the active motor positions, and maintains stable flight.
* **Pushing the Limit:** Press `1` and `4` simultaneously. With two diagonal motors disabled, the allocation matrix determinant approaches zero ($\det(\mathbf{B}) \to 0$), causing loss of control and rotational tumbling.

#### 2. High-Speed Ground Impact & Structural Disassembly

* **Action:** Climb to altitude ($>30\text{m}$), press `F` to invert thrust downwards, and accelerate straight into the terrain.
* **What Happens:** Tip collision forces exceed `StructuralBreakForce` ($2500\text{ N}$), shearing off individual rotor assemblies. The RK4 integrator dynamically updates the remaining mass ($M$), recalculates the shifted CoM, recalibrates the $3 \times 3$ Inertia Tensor, and spawns physical spinning debris. The impact registers a deep crater in the ground mesh via the dynamic `craters` shader array.

#### 3. Heavy Payload Slingshot Dynamics

* **Action:** Lower the winch (`Z`), activate the magnet (`X`), attach a heavy $20\text{kg}+$ payload from the scattered objects, and enter the **Turbulence Volume** or orbit the **Gravity Anomaly**.
* **What Happens:** Exceeding standard payload thresholds tests the feed-forward compensation limits. High angular velocity turns the hanging mass into a double pendulum, transferring dynamic momentum back into the flight frame through local offset torques $\mathbf{r}_{\text{hook}} \times \mathbf{F}_{\text{cable}}$.

#### 4. Relativistic Portal Traversal

* **Action:** Fly at maximum velocity through one of the paired spatial portals.
* **What Happens:** The simulator extracts local linear and angular velocity vectors in the entry portal's reference frame, projects them through a $180^\circ$ coordinate inversion matrix, and reconstructs the body's spatial orientation at the destination portal without resetting integrator state matrices.
"""
