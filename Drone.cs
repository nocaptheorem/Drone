using Godot;
using System.Net.Sockets;
using System.Text;

namespace ProceduralPhysicsLab
{
    // =================================================================================
    // COMPONENT 1: FLIGHT STATE & 3D TENSOR RK4 INTEGRATOR
    // =================================================================================
    public struct RigidBodyState
    {
        public Vector3 Position;
        public Vector3 Velocity;
        public Quaternion Orientation;
        public Vector3 AngularVelocity;
    }

    public struct StateDerivative
    {
        public Vector3 DPosition;
        public Vector3 DVelocity;
        public Quaternion DOrientation;
        public Vector3 DAngularVelocity;
    }

    public class RK4Integrator
    {
        public float Mass = 2.0f;
        public Basis InertiaTensor = new Basis(
            new Vector3(0.05f, 0f, 0f),
            new Vector3(0f, 0.08f, 0f),
            new Vector3(0f, 0f, 0.05f)
        );
        public Basis InvInertiaTensor = new Basis(
            new Vector3(20.0f, 0f, 0f),
            new Vector3(0f, 12.5f, 0f),
            new Vector3(0f, 0f, 20.0f)
        );

        private const float MaxStableSpin = 180.0f;
        private const float MaxStableVelocity = 200.0f;

        public void RecalculateInertia(float brokenPartMass, Vector3 brokenPartLocalPos)
        {
            float newMass = Mass - brokenPartMass;

            Vector3 comShift = (-brokenPartMass * brokenPartLocalPos) / newMass;

            float bx = brokenPartLocalPos.X;
            float by = brokenPartLocalPos.Y;
            float bz = brokenPartLocalPos.Z;
            Basis brokenInertia = new Basis(
                new Vector3(by * by + bz * bz, -bx * by, -bx * bz) * brokenPartMass,
                new Vector3(-bx * by, bx * bx + bz * bz, -by * bz) * brokenPartMass,
                new Vector3(-bx * bz, -by * bz, bx * bx + by * by) * brokenPartMass
            );

            Basis remainingInertiaOldCom = new Basis(
                InertiaTensor.X - brokenInertia.X,
                InertiaTensor.Y - brokenInertia.Y,
                InertiaTensor.Z - brokenInertia.Z
            );

            float dx = comShift.X;
            float dy = comShift.Y;
            float dz = comShift.Z;
            Basis shiftInertia = new Basis(
                new Vector3(dy * dy + dz * dz, -dx * dy, -dx * dz) * newMass,
                new Vector3(-dx * dy, dx * dx + dz * dz, -dy * dz) * newMass,
                new Vector3(-dx * dz, -dy * dz, dx * dx + dy * dy) * newMass
            );

            InertiaTensor = new Basis(
                remainingInertiaOldCom.X - shiftInertia.X,
                remainingInertiaOldCom.Y - shiftInertia.Y,
                remainingInertiaOldCom.Z - shiftInertia.Z
            );

            Mass = newMass;
            InvInertiaTensor = InertiaTensor.Inverse();
        }

        public RigidBodyState Integrate(RigidBodyState state, float dt, Vector3 forceGlobal, Vector3 torqueLocal)
        {
            state.AngularVelocity = state.AngularVelocity.LimitLength(MaxStableSpin);

            StateDerivative k1 = Evaluate(state, 0.0f, new StateDerivative(), forceGlobal, torqueLocal);
            StateDerivative k2 = Evaluate(state, dt * 0.5f, k1, forceGlobal, torqueLocal);
            StateDerivative k3 = Evaluate(state, dt * 0.5f, k2, forceGlobal, torqueLocal);
            StateDerivative k4 = Evaluate(state, dt, k3, forceGlobal, torqueLocal);

            Vector3 dPos = (k1.DPosition + 2.0f * (k2.DPosition + k3.DPosition) + k4.DPosition) / 6.0f;
            Vector3 dVel = (k1.DVelocity + 2.0f * (k2.DVelocity + k3.DVelocity) + k4.DVelocity) / 6.0f;
            Vector3 dAngVel = (k1.DAngularVelocity + 2.0f * (k2.DAngularVelocity + k3.DAngularVelocity) + k4.DAngularVelocity) / 6.0f;

            Quaternion dOri = new Quaternion(
                (k1.DOrientation.X + 2.0f * (k2.DOrientation.X + k3.DOrientation.X) + k4.DOrientation.X) / 6.0f,
                (k1.DOrientation.Y + 2.0f * (k2.DOrientation.Y + k3.DOrientation.Y) + k4.DOrientation.Y) / 6.0f,
                (k1.DOrientation.Z + 2.0f * (k2.DOrientation.Z + k3.DOrientation.Z) + k4.DOrientation.Z) / 6.0f,
                (k1.DOrientation.W + 2.0f * (k2.DOrientation.W + k3.DOrientation.W) + k4.DOrientation.W) / 6.0f
            );

            Quaternion newOri = new Quaternion(
                state.Orientation.X + dOri.X * dt,
                state.Orientation.Y + dOri.Y * dt,
                state.Orientation.Z + dOri.Z * dt,
                state.Orientation.W + dOri.W * dt
            );

            if (newOri.LengthSquared() < 0.0001f) newOri = Quaternion.Identity;
            else newOri = newOri.Normalized();

            Vector3 finalVel = state.Velocity + dVel * dt;
            if (finalVel.LengthSquared() > MaxStableVelocity * MaxStableVelocity)
                finalVel = finalVel.Normalized() * MaxStableVelocity;

            return new RigidBodyState
            {
                Position = state.Position + dPos * dt,
                Velocity = finalVel,
                AngularVelocity = (state.AngularVelocity + dAngVel * dt).LimitLength(MaxStableSpin),
                Orientation = newOri
            };
        }

        private StateDerivative Evaluate(RigidBodyState initial, float dt, StateDerivative d, Vector3 forceGlobal, Vector3 torqueLocal)
        {
            Vector3 newVel = initial.Velocity + d.DVelocity * dt;
            if (newVel.LengthSquared() > MaxStableVelocity * MaxStableVelocity)
                newVel = newVel.Normalized() * MaxStableVelocity;

            Vector3 newAngVel = initial.AngularVelocity + d.DAngularVelocity * dt;
            if (newAngVel.LengthSquared() > MaxStableSpin * MaxStableSpin)
                newAngVel = newAngVel.Normalized() * MaxStableSpin;

            Quaternion newOri = new Quaternion(
                initial.Orientation.X + d.DOrientation.X * dt,
                initial.Orientation.Y + d.DOrientation.Y * dt,
                initial.Orientation.Z + d.DOrientation.Z * dt,
                initial.Orientation.W + d.DOrientation.W * dt
            );

            if (newOri.LengthSquared() < 0.0001f) newOri = Quaternion.Identity;
            else newOri = newOri.Normalized();

            RigidBodyState s = new RigidBodyState
            {
                Position = initial.Position + d.DPosition * dt,
                Velocity = newVel,
                AngularVelocity = newAngVel,
                Orientation = newOri
            };

            StateDerivative output = new StateDerivative
            {
                DPosition = s.Velocity,
                DVelocity = forceGlobal / Mass,
                DOrientation = new Quaternion(
                    0.5f * ( s.Orientation.W * s.AngularVelocity.X + s.Orientation.Y * s.AngularVelocity.Z - s.Orientation.Z * s.AngularVelocity.Y),
                    0.5f * ( s.Orientation.W * s.AngularVelocity.Y + s.Orientation.Z * s.AngularVelocity.X - s.Orientation.X * s.AngularVelocity.Z),
                    0.5f * ( s.Orientation.W * s.AngularVelocity.Z + s.Orientation.X * s.AngularVelocity.Y - s.Orientation.Y * s.AngularVelocity.X),
                    0.5f * (-s.Orientation.X * s.AngularVelocity.X - s.Orientation.Y * s.AngularVelocity.Y - s.Orientation.Z * s.AngularVelocity.Z)
                )
            };

            Vector3 angularMomentum = InertiaTensor * s.AngularVelocity;
            Vector3 gyroTorque = s.AngularVelocity.Cross(angularMomentum);

            output.DAngularVelocity = InvInertiaTensor * (torqueLocal - gyroTorque);

            return output;
        }
    }

    // =================================================================================
    // COMPONENT 2: CASCADED CONTROLLERS & ALTITUDE PID
    // =================================================================================
    public class CascadedAttitudeController
    {
        public float KpOuter = 10.0f;
        public Vector3 KpRate = new Vector3(12.0f, 12.0f, 10.0f);
        public Vector3 KdRate = new Vector3(2.0f, 2.0f, 1.5f);

        public Vector3 ComputeTorque(Quaternion currentRot, Quaternion targetRot, Vector3 currentAngVelLocal)
        {
            Quaternion qErr = currentRot.Inverse() * targetRot;
            if (qErr.W < 0) qErr = new Quaternion(-qErr.X, -qErr.Y, -qErr.Z, -qErr.W);
            qErr = qErr.Normalized();

            Vector3 axis = new Vector3(qErr.X, qErr.Y, qErr.Z);
            float angle = 2.0f * Mathf.Acos(Mathf.Clamp(qErr.W, -1.0f, 1.0f));

            if (axis.LengthSquared() > 0.0001f) axis = axis.Normalized();
            else axis = Vector3.Zero;

            Vector3 targetAngVelLocal = axis * angle * KpOuter;
            Vector3 rateError = targetAngVelLocal - currentAngVelLocal;

            return new Vector3(
                (rateError.X * KpRate.X) - (currentAngVelLocal.X * KdRate.X),
                (rateError.Y * KpRate.Y) - (currentAngVelLocal.Y * KdRate.Y),
                (rateError.Z * KpRate.Z) - (currentAngVelLocal.Z * KdRate.Z)
            );
        }
    }

    public class PIDController
    {
        public float Kp, Ki, Kd;
        private float _integralAccumulator, _integralLimit;
        public PIDController(float p, float i, float d, float limit = 20.0f) { Kp = p; Ki = i; Kd = d; _integralLimit = limit; }

        public float Update(float error, float currentVelocity, float dt, bool freezeIntegral = false)
        {
            if (dt <= 0.0001f) return 0f;
            if (!freezeIntegral)
            {
                _integralAccumulator += error * dt;
                _integralAccumulator = Mathf.Clamp(_integralAccumulator, -_integralLimit, _integralLimit);
            }
            return (Kp * error) + (Ki * _integralAccumulator) + (Kd * (0f - currentVelocity));
        }

        public void Reset() => _integralAccumulator = 0f;
        public void Bleed(float rate, float dt) => _integralAccumulator = Mathf.Lerp(_integralAccumulator, 0f, rate * dt);
    }

    // =================================================================================
    // COMPONENT 3: PROCEDURAL TURBULENCE VOLUME
    // =================================================================================
    public partial class TurbulenceVolume : Area3D
    {
        public Vector3 WindDirection = Vector3.Up;
        public float BaseWindForce = 15.0f;
        public float TurbulenceStrength = 8.0f;
        public Vector3 VolumeSize = new Vector3(10, 10, 10);
        public Color ZoneColor = new Color(0.8f, 0.2f, 0.1f, 0.2f);

        private FastNoiseLite _noise = null!;
        private MeshInstance3D _visualizer = null!;
        private StandardMaterial3D _material = null!;
        private double _timePassed;
        private HashSet<Node3D> _trackedBodies = new HashSet<Node3D>();

        public override void _Ready()
        {
            Monitoring = true;
            AddChild(new CollisionShape3D { Shape = new BoxShape3D { Size = VolumeSize } });

            _material = new StandardMaterial3D
            {
                AlbedoColor = ZoneColor, Transparency = BaseMaterial3D.TransparencyEnum.Alpha,
                DisableReceiveShadows = true, EmissionEnabled = true, Emission = ZoneColor, EmissionEnergyMultiplier = 0.5f
            };

            _visualizer = new MeshInstance3D { Mesh = new BoxMesh { Size = VolumeSize }, MaterialOverride = _material };
            AddChild(_visualizer);

            _noise = new FastNoiseLite { NoiseType = FastNoiseLite.NoiseTypeEnum.Simplex, Frequency = 0.5f };

            BodyEntered += (Node3D body) => { _trackedBodies.Add(body); };
            BodyExited += (Node3D body) => { _trackedBodies.Remove(body); };
        }

        public override void _PhysicsProcess(double delta)
        {
            _timePassed += delta;
            _material.AlbedoColor = new Color(ZoneColor.R, ZoneColor.G, ZoneColor.B, (Mathf.Sin((float)_timePassed * 2.0f) * 0.1f) + 0.2f);

            foreach (var node in _trackedBodies)
            {
                if (node.Name == "DroneBody" && node.GetParent() is Drone builder)
                {
                    Vector3 pos = node.GlobalPosition;
                    float nx = _noise.GetNoise3D(pos.X, pos.Y, pos.Z + (float)_timePassed * 5f);
                    float ny = _noise.GetNoise3D(pos.X + 100, pos.Y, pos.Z + (float)_timePassed * 5f);
                    float nz = _noise.GetNoise3D(pos.X, pos.Y + 100, pos.Z + (float)_timePassed * 5f);

                    Vector3 turbulenceForce = new Vector3(nx, ny, nz) * TurbulenceStrength;
                    Vector3 totalForce = (WindDirection * BaseWindForce) + turbulenceForce;
                    Vector3 turbulentOffset = new Vector3(nz, 0, nx) * 0.5f;

                    builder.ApplyExternalWind(totalForce, turbulentOffset);
                }
            }
        }
    }

    // =================================================================================
    // MAIN CLASS: PROCEDURAL PHYSICS LAB - FLIGHT CONTROLLER
    // =================================================================================
    [GlobalClass]
    public partial class Drone : Node3D
    {
        [ExportGroup("Telemetry & Debug")]
        [Export] public bool EnableTelemetry = true;
        [Export] public float TelemetryPrintRateHz = 4.0f;
        [Export] public bool EnableAnomalyDetector = true;
        private UdpClient _udpClient;
        private const string UDP_IP = "127.0.0.1";
        private const int UDP_PORT = 9870;
        private float _lastUdpErrorTime = -10.0f;
        private const float ERROR_LOG_INTERVAL_SEC = 2.0f; // Limit error prints to avoid spamming console

        [ExportGroup("Quadcopter Config")]
        [Export] public float ArmLength = 0.5f;
        [Export] public float MaxMotorThrust = 40.0f;
        [Export] public float YawDragFactor = 0.4f;
        [Export] public float MotorTimeConstant = 15.0f;
        [Export] public float InputSlewSpeed = 8.0f;

        [ExportGroup("Gyroscopic Physics")]
        [Export] public float RotorMass = 0.05f;
        [Export] public float RotorRadius = 0.4f;
        [Export] public float RpmScaleFactor = 150.0f;
        [Export] public float StructuralBreakForce = 1000.0f;
        [Export] public float RotorBreakageMultiplier = 1.0f;
        [Export] public bool EnableGravitationalAnomaly = false;


        [ExportGroup("Deformable Terrain")]
        [Export] public float MinimumCraterForce = 1500.0f;
        [Export] public float CraterForceSensitivity = 0.001f;
        [Export] public float MaxCraterRadius = 2.5f;

        [ExportGroup("Camera System")]
        [Export] public bool CameraFollowsDrone = true;
        [Export] public float MaxStabilizationAngle = 60.0f;
        [Export] public float BaseFOV = 75.0f;
        [Export] public float MaxFOV = 110.0f;
        [Export] public float FovSpeedThreshold = 30.0f;
        [Export] public float FovInterpolationRate = 5.0f;

        [ExportGroup("Winch System")]
        [Export] public float MaxCableLength = 20.0f;
        [Export] public float CableReelSpeed = 4.0f;
        [Export] public float CableSpringStiffness = 2500.0f;
        [Export] public float CableDampingRatio = 0.85f;
        [Export] public float HookCaptureRadius = 1.5f;

        private AnimatableBody3D _gravityAnomaly = null!;
        private float _anomalyOrbitAngle = 0.0f;
        private RigidBody3D _magnetBody = null!;
        private RigidBody3D _hookedPayload = null;

        private float _currentCableLength = 0.0f;
        private bool _electromagnetActive = false;
        private MeshInstance3D _cableVisual;
        private MeshInstance3D _magnetVisual;
        private Vector3 _hookLocalOffset = new Vector3(0, -0.1f, 0);
        private bool _anomalyActive = true;

        private RK4Integrator _rk4 = new RK4Integrator();
        private RigidBodyState _simState;

        private CascadedAttitudeController _attController = new CascadedAttitudeController();
        private PIDController _altPID = new PIDController(15.0f, 0.021777f, 9.0f, 200.0f);

        private float _targetAlt = 5.0f;
        private Quaternion _targetAttitude = Quaternion.Identity;
        private float _targetYawAccumulator = 0.0f;
        private float _smoothedPitch = 0.0f;
        private float _smoothedRoll = 0.0f;

        private AnimatableBody3D _drone = null!;
        private Camera3D _camera = null!;
        private Label _hud = null!;
        private Node3D[] _rotors = new Node3D[4];
        private Vector3[] _rotorPositions = new Vector3[4];

        private bool[] _motorActive = { true, true, true, true };
        private bool[] _rotorStructurallyIntact = { true, true, true, true };
        private float[] _actualMotorThrust = new float[4];
        private MeshInstance3D[] _thrustVectors = new MeshInstance3D[4];

        private GpuParticles3D[] _smokeEmitters = new GpuParticles3D[4];
        private GpuParticles3D _sparkEmitter = null!;
        private AudioStreamPlayer3D[] _motorAudio = new AudioStreamPlayer3D[4];

        private Area3D _downwashVolume = null!;
        private HashSet<RigidBody3D> _bodiesInDownwash = new HashSet<RigidBody3D>();

        private Godot.Collections.Array<Vector4> _craters = new Godot.Collections.Array<Vector4>();
        private int _craterIndex = 0;
        private ShaderMaterial _floorMaterial = null!;

        private float _thrustDirection = 1.0f;
        private float _telemetryLy = 0.0f;
        private Vector3 _telemetryGyroTorque = Vector3.Zero;
        private Vector3 _externalForceAccum = Vector3.Zero;
        private Vector3 _externalTorqueAccum = Vector3.Zero;

        private float _camYaw = 0.0f;
        private float _camPitch = 0.3f;
        private float _camDistance = 6.0f;
        private Shader _portalShader = null!;
        private Shader _floorShader = null!;
        private Shader _smokeShader = null!;
        private float _lastDroneYaw = 0.0f;

        private List<RigidBody3D> _activeDebris = new List<RigidBody3D>();
        private Shader _plumeShader = null!;

        private float _telemetryTimer = 0.0f;
        private float _lastWinchDist = 0.0f;
        private Vector3 _lastWinchForce = Vector3.Zero;
        private float _lastPidOutput = 0.0f;
        private Vector3 _lastExtForce = Vector3.Zero;
        private Vector3 _lastExtTorque = Vector3.Zero;
        private Vector3 _lastTorqueCmd = Vector3.Zero;

        private const string FLUID_PLUME_SHADER = @"
          shader_type spatial;
        render_mode blend_add, unshaded, cull_disabled, depth_draw_never;

        float hash(vec2 p) { return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453); }

        float noise(vec2 p) {
          vec2 i = floor(p); vec2 f = fract(p); f = f * f * (3.0 - 2.0 * f);
          return mix(mix(hash(i + vec2(0.0, 0.0)), hash(i + vec2(1.0, 0.0)), f.x),
              mix(hash(i + vec2(0.0, 1.0)), hash(i + vec2(1.0, 1.0)), f.x), f.y);
        }

        void vertex() {
          float expansion = pow(UV.y, 2.5);
          float flutter = noise(vec2(VERTEX.y * 15.0, TIME * 25.0));
          vec3 push = normalize(vec3(NORMAL.x, 0.0, NORMAL.z));
          VERTEX += push * (flutter * 0.35 * expansion);
        }

        void fragment() {
          float t1 = mod(TIME * 12.0, 1000.0);
          float t2 = mod(TIME * 22.0, 1000.0);

          vec2 uv1 = vec2(UV.x * 3.0, UV.y * 1.5 - t1);
          vec2 uv2 = vec2(UV.x * 5.0 + mod(TIME * 3.0, 100.0), UV.y * 2.5 - t2);

          float n = noise(uv1) * 0.6 + noise(uv2) * 0.4;
          float vertical_fade = pow(UV.y, 1.8);
          float horizontal_fade = sin(UV.x * 3.141592);
          horizontal_fade = pow(horizontal_fade, 2.0);

          float fire_mask = smoothstep(0.1, 0.8, n - (vertical_fade * 0.8));
          fire_mask *= horizontal_fade;

          vec3 core = vec3(1.0, 0.95, 0.75);
          vec3 edge = vec3(1.0, 0.25, 0.0);

          ALBEDO = mix(edge, core, fire_mask * 0.8) * 2.5;
          ALPHA = fire_mask * (1.0 - vertical_fade) * 0.75;
        }
        ";

        private const string PULSAR_CORE_SHADER = @"
          shader_type spatial;
        render_mode unshaded;
        uniform sampler2D noise_tex_a;
        uniform sampler2D noise_tex_b;
        void fragment() {
          float n1 = texture(noise_tex_a, UV + TIME * 0.03).r;
          float n2 = texture(noise_tex_b, UV - TIME * 0.02).r;
          float combined_noise = (n1 + n2) * 0.5;
          vec3 dark_base = vec3(0.05, 0.0, 0.15);
          vec3 mid_glow = vec3(0.3, 0.1, 0.8);
          vec3 bright_scales = vec3(0.5, 0.9, 1.0);
          vec3 color = mix(dark_base, mid_glow, smoothstep(0.2, 0.5, combined_noise));
          color = mix(color, bright_scales, smoothstep(0.6, 0.9, combined_noise));
          float fresnel = pow(1.0 - dot(NORMAL, VIEW), 2.5);
          ALBEDO = color + (bright_scales * fresnel * 0.5);
        }";

        private const string VIOLET_VORTEX_DISK_SHADER = @"
          shader_type spatial;
        render_mode unshaded, cull_disabled, blend_add;
        uniform sampler2D noise_tex;
        void fragment() {
          vec2 rel_uv = UV - vec2(0.5);
          float dist = length(rel_uv);
          float angle = atan(rel_uv.y, rel_uv.x);
          float noise = texture(noise_tex, vec2(dist - TIME * 0.1, angle * 0.2 + TIME * 0.05)).r;

          // Doppler effect creates the illusion of a spinning vortex
          float doppler = dot(normalize(rel_uv), vec2(1.0, 0.2)) * 0.5 + 0.5;

          vec3 deep_purple = vec3(0.5, 0.1, 1.0) * 8.0;
          vec3 cyan_shift = vec3(0.1, 0.8, 1.0) * 4.0;
          vec3 final_color = mix(cyan_shift, deep_purple, doppler);

          // Fades out at the edges and hollows out the center for the core
          float alpha_mask = smoothstep(0.5, 0.2, dist) * smoothstep(0.05, 0.15, dist);

          ALBEDO = final_color * noise;
          ALPHA = alpha_mask * noise * (doppler + 0.3);
        }";

        private void UpdateAnomalyPhysics(float dt)
        {
          if (_gravityAnomaly == null) return;


          // 1. Unstoppable "On-Rails" Orbit
          float orbitRadius = 80.0f;
          float orbitSpeed = 0.08f; // Radians per second - adjust for faster/slower orbits

          _anomalyOrbitAngle += dt * orbitSpeed;

          _gravityAnomaly.GlobalPosition = new Vector3(
              Mathf.Cos(_anomalyOrbitAngle) * orbitRadius,
              20.0f,
              Mathf.Sin(_anomalyOrbitAngle) * orbitRadius
              );

          // Add a visual spin
          _gravityAnomaly.RotateY(dt * 2.5f);

          // 2. Pull the Drone (Tuned for Soft Orbits)
          Vector3 toAnomaly = _gravityAnomaly.GlobalPosition - _simState.Position;
          float dist = toAnomaly.Length();

          if (dist < 250.0f)
          {
            // --- TWEAK 1: Lower overall strength and flatten the core ---
            float gravityStrength = 80000.0f; // Lowered from 80000

            // A massive epsilon flattens the gravity well. Instead of spiking to infinity
            // at the center, gravity peaks at a distance, then smoothly drops to zero.
            float epsilon = 20.0f;

            float gravityAccel = gravityStrength / (dist * dist + epsilon * epsilon);

            // Note: The 'dist < 50.0f' proximity multiplier block has been DELETED entirely.

            // --- TWEAK 2: The "Vortex Assist" ---
            // Pure gravity pulls straight in. If you don't hit it with perfect tangential
            // velocity, you fall in. This adds a gentle horizontal swirl to naturally
            // guide your drone into a circular path.
            Vector3 pullDir = toAnomaly.Normalized();
            Vector3 orbitTangent = pullDir.Cross(Vector3.Up).Normalized();

            // Mix 100% inward pull with 40% sideways swirl
            Vector3 forceDir = (pullDir + (orbitTangent * 0.4f)).Normalized();

            float totalForceMag = gravityAccel * _rk4.Mass;

            // Cap the force lower so you can easily throttle out of it
            float maxDroneMotorForce = MaxMotorThrust * 1.8f;
            float cappedForceMag = Mathf.Min(totalForceMag, maxDroneMotorForce);

            // Inject into the RK4 Integrator using the new swirling direction
            ApplyExternalWind(forceDir * cappedForceMag, Vector3.Zero);

            // 3. Pull Debris (also using the soft parameters so they orbit smoothly too)
            foreach (var debris in _activeDebris)
            {
              if (IsInstanceValid(debris))
              {
                Vector3 debrisToAnomaly = _gravityAnomaly.GlobalPosition - debris.GlobalPosition;
                float dDist = debrisToAnomaly.Length();
                if (dDist < 100.0f) {
                  Vector3 dPull = debrisToAnomaly.Normalized();
                  Vector3 dTangent = dPull.Cross(Vector3.Up).Normalized();
                  Vector3 dForceDir = (dPull + (dTangent * 0.6f)).Normalized();

                  float dForceMag = (gravityStrength * debris.Mass) / (dDist * dDist + epsilon * epsilon);
                  debris.ApplyCentralForce(dForceDir * dForceMag);
                }
              }
            }
          }
        }

        private void BuildGravityAnomaly()
        {
          _gravityAnomaly = new AnimatableBody3D {
            Position = new Vector3(180, 20, 0),
            SyncToPhysics = false
          };

          // 1. Generate Noise Textures dynamically for the shaders
          var noiseA = new NoiseTexture2D {
            Noise = new FastNoiseLite { Frequency = 0.02f, Seed = 100 }, Seamless = true
          };
          var noiseB = new NoiseTexture2D {
            Noise = new FastNoiseLite { Frequency = 0.03f, Seed = 200 }, Seamless = true
          };

          // 2. Build the Pulsar Core (The physical center)
          var coreMat = new ShaderMaterial { Shader = new Shader { Code = PULSAR_CORE_SHADER } };
          coreMat.SetShaderParameter("noise_tex_a", noiseA);
          coreMat.SetShaderParameter("noise_tex_b", noiseB);

          var coreMesh = new MeshInstance3D {
            Mesh = new SphereMesh { Radius = 4.0f },
            MaterialOverride = coreMat
          };
          _gravityAnomaly.AddChild(coreMesh);

          // 3. Build the Swirling Accretion Disk (Visualizing the vortex pull)
          var diskMat = new ShaderMaterial { Shader = new Shader { Code = VIOLET_VORTEX_DISK_SHADER } };
          diskMat.SetShaderParameter("noise_tex", noiseA);

          var diskMesh = new MeshInstance3D {
            // A flat plane is perfect for UV mapping the accretion disk shader
            Mesh = new PlaneMesh { Size = new Vector2(30.0f, 30.0f) },
            MaterialOverride = diskMat
          };
          _gravityAnomaly.AddChild(diskMesh);

          // 4. Physical Collision
          _gravityAnomaly.AddChild(new CollisionShape3D { Shape = new SphereShape3D { Radius = 4.0f } });

          AddChild(_gravityAnomaly);
        }

        public override void _Ready()
        {
            Engine.PhysicsTicksPerSecond = 120;
            _simState = new RigidBodyState { Position = new Vector3(0, 5, 0), Orientation = Quaternion.Identity };

            for (int i = 0; i < 16; i++) _craters.Add(Vector4.Zero);

            CompileShaders();
            SetupEnvironment();
            BuildPortals();
            BuildDrone();
            BuildWinchSystem();
            if (EnableGravitationalAnomaly) {
              BuildGravityAnomaly();
            }
            ScatterPhysicsObjects();
            SetupDownwashVolume();
            SetupHUD();
            Input.MouseMode = Input.MouseModeEnum.Captured;
            try
            {
              _udpClient = new UdpClient();
            }
            catch (Exception ex)
            {
              GD.PrintErr($"[TELEMETRY INIT ERROR] Failed to instantiate UdpClient: {ex.Message}");
            }
        }

        private void BuildWinchSystem()
        {
            _cableVisual = new MeshInstance3D
            {
                Mesh = new CylinderMesh { TopRadius = 0.015f, BottomRadius = 0.015f, Height = 1.0f },
                MaterialOverride = new StandardMaterial3D { AlbedoColor = Colors.DarkSlateGray, Metallic = 0.8f }
            };
            AddChild(_cableVisual);

            _magnetBody = new RigidBody3D
            {
                Mass = 3.5f,
                Position = _simState.Position + _hookLocalOffset,
                LinearDamp = 0.2f,
                AngularDamp = 0.5f
            };

            var magnetShape = new CollisionShape3D { Shape = new CylinderShape3D { Radius = 0.15f, Height = 0.08f } };
            _magnetBody.AddChild(magnetShape);

            _magnetVisual = new MeshInstance3D
            {
                Mesh = new CylinderMesh { TopRadius = 0.15f, BottomRadius = 0.15f, Height = 0.08f },
                MaterialOverride = new StandardMaterial3D
                {
                    AlbedoColor = Colors.Black,
                    EmissionEnabled = true,
                    Emission = Colors.Red,
                    EmissionEnergyMultiplier = 0.0f
                }
            };
            _magnetBody.AddChild(_magnetVisual);

            AddChild(_magnetBody);

            if (_drone != null)
            {
                _magnetBody.AddCollisionExceptionWith(_drone);
                _drone.AddCollisionExceptionWith(_magnetBody);
            }
        }

        private void ScatterPhysicsObjects()
        {
          RandomNumberGenerator rng = new RandomNumberGenerator();
          rng.Randomize();

          Color[] palette = { Colors.Coral, Colors.Teal, Colors.Gold, Colors.RebeccaPurple, Colors.DodgerBlue };

          for (int i = 0; i < 0; i++)
          {
            // 1. Generate the physical size first (from small 0.5m chunks to massive 3.0m boulders)
            float size = rng.RandfRange(0.5f, 3.0f);

            // 2. Calculate mass based on volume.
            // Volume scales with the cube of the size. We multiply by a base density.
            float density = 2.5f;
            float calculatedMass = Mathf.Pow(size, 3.0f) * density;

            RigidBody3D rb = new RigidBody3D
            {
              Mass = calculatedMass,
              Position = new Vector3(rng.RandfRange(-40, 40), 10.0f + (i * 2.0f), rng.RandfRange(-40, 40)),
              LinearDamp = 0.5f,
              AngularDamp = 0.5f
            };

            Mesh mesh = null;
            Shape3D shape = null;
            int type = rng.RandiRange(0, 2);

            if (type == 0) // Cube
            {
              mesh = new BoxMesh { Size = Vector3.One * size };
              shape = new BoxShape3D { Size = Vector3.One * size };
            }
            else if (type == 1) // Sphere
            {
              mesh = new SphereMesh { Radius = size / 2.0f, Height = size };
              shape = new SphereShape3D { Radius = size / 2.0f };
            }
            else // Cylinder
            {
              mesh = new CylinderMesh { TopRadius = size / 2.0f, BottomRadius = size / 2.0f, Height = size };
              shape = new CylinderShape3D { Radius = size / 2.0f, Height = size };
            }

            var mat = new StandardMaterial3D
            {
              AlbedoColor = palette[rng.RandiRange(0, palette.Length - 1)],
              Roughness = rng.RandfRange(0.1f, 0.9f),
              Metallic = rng.RandfRange(0.0f, 0.8f)
            };

            rb.AddChild(new MeshInstance3D { Mesh = mesh, MaterialOverride = mat });
            rb.AddChild(new CollisionShape3D { Shape = shape });

            AddChild(rb);
          }
        }

        // Phase 1: Cleaned up magnet catch/release logic (Cable integration moved to ApplyFlightDynamics)
        private void UpdateElectromagnet()
        {
            Vector3 magnetPos = _magnetBody.GlobalPosition;

            if (_electromagnetActive)
            {
                ((StandardMaterial3D)_magnetVisual.MaterialOverride).EmissionEnergyMultiplier = _hookedPayload != null ? 8.0f : 4.0f;
                ((StandardMaterial3D)_magnetVisual.MaterialOverride).Emission = _hookedPayload != null ? Colors.LimeGreen : Colors.Red;

                if (_hookedPayload == null)
                {
                    var spaceState = GetWorld3D().DirectSpaceState;
                    var shapeParams = new PhysicsShapeQueryParameters3D
                    {
                        Shape = new SphereShape3D { Radius = HookCaptureRadius },
                        Transform = new Transform3D(Basis.Identity, magnetPos)
                    };

                    var results = spaceState.IntersectShape(shapeParams);
                    foreach (var res in results)
                    {
                        if (res["collider"].AsGodotObject() is RigidBody3D rb && rb.Name != "DroneBody" && rb != _magnetBody)
                        {
                            _hookedPayload = rb;
                            _hookedPayload.AddCollisionExceptionWith(_drone);
                            _hookedPayload.AddCollisionExceptionWith(_magnetBody);
                            break;
                        }
                    }
                }
                else if (IsInstanceValid(_hookedPayload))
                {
                    Vector3 payloadPos = _hookedPayload.GlobalPosition;
                    Vector3 diff = magnetPos - payloadPos;
                    float payloadDist = diff.Length();

                    if (payloadDist > 0.05f)
                    {
                        Vector3 dir = diff / payloadDist;
                        Vector3 relVel = _magnetBody.LinearVelocity - _hookedPayload.LinearVelocity;

                        // Critical damping logic for the payload capture
                        float m_eff = (_magnetBody.Mass * _hookedPayload.Mass) / (_magnetBody.Mass + _hookedPayload.Mass);
                        float k = 5000.0f;
                        float c = 2.0f * Mathf.Sqrt(k * m_eff);

                        float pullForce = (payloadDist * k) + (relVel.Dot(dir) * c);
                        pullForce = Mathf.Clamp(pullForce, 0.0f, 8000.0f);

                        Vector3 forceVec = dir * pullForce;
                        _hookedPayload.ApplyCentralForce(forceVec);
                        _magnetBody.ApplyCentralForce(-forceVec);
                    }
                }
            }
            else
            {
                ((StandardMaterial3D)_magnetVisual.MaterialOverride).EmissionEnergyMultiplier = 0.0f;
                ((StandardMaterial3D)_magnetVisual.MaterialOverride).Emission = Colors.Red;

                if (_hookedPayload != null && IsInstanceValid(_hookedPayload))
                {
                    _hookedPayload.RemoveCollisionExceptionWith(_drone);
                    _hookedPayload.RemoveCollisionExceptionWith(_magnetBody);
                    _hookedPayload = null;
                }
                else if (_hookedPayload != null)
                {
                    _hookedPayload = null;
                }
            }
        }

        private void UpdateDebrisPhysics(float dt)
        {
          const float debrisHalfHeight = 0.01f;
          const float groundClearance  = 0.05f;

          for (int i = _activeDebris.Count - 1; i >= 0; i--)
          {
            var d = _activeDebris[i];
            if (!IsInstanceValid(d) || d.Sleeping)
            {
              _activeDebris.RemoveAt(i);
              continue;
            }

            float groundY = GetTerrainHeight(new Vector2(d.GlobalPosition.X, d.GlobalPosition.Z));
            float debrisBottomY = d.GlobalPosition.Y - debrisHalfHeight;

            if (debrisBottomY < groundY)
            {
              var trans = d.GlobalTransform;
              trans.Origin.Y = groundY + debrisHalfHeight + groundClearance;
              d.GlobalTransform = trans;

              d.LinearVelocity = new Vector3(
                  d.LinearVelocity.X * 0.7f,
                  Mathf.Max(0, d.LinearVelocity.Y * -0.4f),
                  d.LinearVelocity.Z * 0.7f);
              d.AngularVelocity *= 0.6f;

              if (d.LinearVelocity.LengthSquared() < 0.1f &&
                  d.AngularVelocity.LengthSquared() < 0.5f)
              {
                d.LinearVelocity = Vector3.Zero;
                d.AngularVelocity = Vector3.Zero;
                d.Sleeping = true;
                _activeDebris.RemoveAt(i);
              }
            }
          }
        }

        public override void _PhysicsProcess(double delta)
        {
            float dt = (float)delta;

            HandleKeyboardInput(dt);
            UpdateElectromagnet();
            if (EnableGravitationalAnomaly) {
              UpdateAnomalyPhysics(dt);
            }

            ApplyFlightDynamics(dt);

            ApplyRotorDownwash(dt);
            UpdateVisualsAndAudio(dt);
            UpdateDebrisPhysics(dt);
            UpdateHUD();

            if (EnableAnomalyDetector) DetectAnomalies();

            if (EnableTelemetry)
            {
                _telemetryTimer += dt;
                if (TelemetryPrintRateHz > 0 && _telemetryTimer >= (1.0f / TelemetryPrintRateHz))
                {
                    LogTelemetry();
                    _telemetryTimer = 0f;
                }
            }

            _externalForceAccum = Vector3.Zero;
            _externalTorqueAccum = Vector3.Zero;
        }

        private void DetectAnomalies()
        {
            if (float.IsNaN(_simState.Position.X) || float.IsInfinity(_simState.Position.X))
                GD.PrintErr("[ANOMALY] NaN or Infinity detected in Drone Position! ", _simState.Position);

            if (float.IsNaN(_simState.Velocity.X) || float.IsInfinity(_simState.Velocity.X))
                GD.PrintErr("[ANOMALY] NaN or Infinity detected in Drone Velocity! ", _simState.Velocity);

            if (_lastWinchForce.LengthSquared() > 20000f * 20000f)
                GD.PrintErr($"[ANOMALY] Extreme Winch Force detected: {_lastWinchForce.Length():F1} N. Threshold breached.");

            if (float.IsNaN(_actualMotorThrust[0]) || float.IsInfinity(_actualMotorThrust[0]))
                GD.PrintErr($"[ANOMALY] Invalid Motor Thrust detected. Value: {_actualMotorThrust[0]}");

            if (Mathf.Abs(_simState.Orientation.LengthSquared() - 1.0f) > 0.05f)
                GD.PrintErr($"[ANOMALY] Quaternion denormalization detected! Length: {_simState.Orientation.Length():F4}");
        }

        private void LogTelemetry()
        {
          if (_udpClient == null)
          {
            ReportTelemetryError("UdpClient is uninitialized or null.");
            return;
          }

          var metrics = new
          {
            timestamp = Time.GetTicksMsec() / 1000.0f,
            alt_target = _targetAlt,
            alt_actual = _simState.Position.Y,
            alt_error = _targetAlt - _simState.Position.Y,
            pid_out_n = _lastPidOutput,
            vel_y = _simState.Velocity.Y,
            cmd_torque_pitch = _lastTorqueCmd.X,
            cmd_torque_yaw = _lastTorqueCmd.Y,
            cmd_torque_roll = _lastTorqueCmd.Z,
            thrust_fl = _actualMotorThrust[0],
            thrust_fr = _actualMotorThrust[1],
            thrust_bl = _actualMotorThrust[2],
            thrust_br = _actualMotorThrust[3],
            winch_force_n = _lastWinchForce.Length()
          };

          try
          {
            string jsonString = System.Text.Json.JsonSerializer.Serialize(metrics);
            byte[] payload = Encoding.UTF8.GetBytes(jsonString);

            _udpClient.Send(payload, payload.Length, UDP_IP, UDP_PORT);
          }
          catch (SocketException ex)
          {
            ReportTelemetryError($"SocketException on port {UDP_PORT}: {ex.Message} (Code: {ex.SocketErrorCode})");
          }
          catch (Exception ex)
          {
            ReportTelemetryError($"Unexpected telemetry serialization/transmission error: {ex.Message}");
          }
        }

        private void ReportTelemetryError(string message)
        {
          float currentTime = Time.GetTicksMsec() / 1000.0f;

          // Rate-limit the error reporting so it doesn't saturate the stdout/stderr stream or drop frames
          if (currentTime - _lastUdpErrorTime >= ERROR_LOG_INTERVAL_SEC)
          {
            GD.PrintErr($"[TELEMETRY ERROR] {message}");
            _lastUdpErrorTime = currentTime;
          }
        }

        public override void _ExitTree()
        {
          _udpClient?.Close();
          _udpClient?.Dispose();
          _udpClient = null;
        }

        public void ApplyExternalWind(Vector3 force, Vector3 localOffset)
        {
            _externalForceAccum += force;
            Vector3 globalOffset = new Basis(_simState.Orientation) * localOffset;
            _externalTorqueAccum += globalOffset.Cross(force);
        }

        public override void _Input(InputEvent @event)
        {
            if (@event is InputEventMouseMotion mouseMotion && Input.MouseMode == Input.MouseModeEnum.Captured)
            {
                _camYaw += mouseMotion.Relative.X * 0.005f;
                _camPitch = Mathf.Clamp(_camPitch + mouseMotion.Relative.Y * 0.005f, -Mathf.Pi / 2.2f, Mathf.Pi / 2.2f);
            }
            if (@event is InputEventKey key && key.Pressed && !key.Echo)
            {
                if (key.Keycode == Key.Key1) _motorActive[0] = !_motorActive[0];
                if (key.Keycode == Key.Key2) _motorActive[1] = !_motorActive[1];
                if (key.Keycode == Key.Key3) _motorActive[2] = !_motorActive[2];
                if (key.Keycode == Key.Key4) _motorActive[3] = !_motorActive[3];
                if (key.Keycode == Key.F) _thrustDirection *= -1.0f;
                if (key.Keycode == Key.X) _electromagnetActive = !_electromagnetActive;
                if (key.Keycode == Key.Tab) CameraFollowsDrone = !CameraFollowsDrone;
                if (key.Keycode == Key.Escape) Input.MouseMode = Input.MouseMode == Input.MouseModeEnum.Captured ? Input.MouseModeEnum.Visible : Input.MouseModeEnum.Captured;
                if (key.Keycode == Key.G) _anomalyActive = !_anomalyActive;
            }
        }

        private void HandleKeyboardInput(float dt)
        {
            int intactCount = 0;
            for (int i = 0; i < 4; i++) if (_rotorStructurallyIntact[i] && _motorActive[i]) intactCount++;

            if (Input.IsKeyPressed(Key.Space)) _targetAlt += 5.0f * dt;
            if (Input.IsKeyPressed(Key.Shift)) _targetAlt -= 5.0f * dt;
            _targetAlt = Mathf.Max(_targetAlt, 0.5f);

            if (Input.IsKeyPressed(Key.Z)) _currentCableLength = Mathf.Min(_currentCableLength + CableReelSpeed * dt, MaxCableLength);
            if (Input.IsKeyPressed(Key.C)) _currentCableLength = Mathf.Max(_currentCableLength - CableReelSpeed * dt, 0.0f);

            float maxTilt = 1.2f;
            float targetPitch = 0, targetRoll = 0, yawRate = 0;

            if (Input.IsKeyPressed(Key.W)) targetPitch = -maxTilt;
            if (Input.IsKeyPressed(Key.S)) targetPitch = maxTilt;
            if (Input.IsKeyPressed(Key.A)) targetRoll = -maxTilt;
            if (Input.IsKeyPressed(Key.D)) targetRoll = maxTilt;
            if (Input.IsKeyPressed(Key.Q)) yawRate = -2.0f;
            if (Input.IsKeyPressed(Key.E)) yawRate = 2.0f;

            _smoothedPitch = Mathf.MoveToward(_smoothedPitch, targetPitch, InputSlewSpeed * dt);
            _smoothedRoll = Mathf.MoveToward(_smoothedRoll, targetRoll, InputSlewSpeed * dt);

            Vector3 forward = new Basis(_simState.Orientation).Z;
            float currentYaw = Mathf.Atan2(forward.X, forward.Z);

            if (intactCount < 4)
            {
                _targetYawAccumulator = currentYaw;
            }
            else
            {
                _targetYawAccumulator += yawRate * dt;
                float yawError = Mathf.AngleDifference(currentYaw, _targetYawAccumulator);
                _targetYawAccumulator = currentYaw + Mathf.Clamp(yawError, -0.5f, 0.5f);
            }

            _targetAttitude = Basis.FromEuler(new Vector3(_smoothedPitch, _targetYawAccumulator, _smoothedRoll)).GetRotationQuaternion();
        }

        private void RegisterCrater(Vector3 pos, float force)
        {
            if (force < MinimumCraterForce) return;
            float radius = Mathf.Min(force * CraterForceSensitivity, MaxCraterRadius);
            _craters[_craterIndex] = new Vector4(pos.X, pos.Y, pos.Z, radius);
            _craterIndex = (_craterIndex + 1) % 16;
            _floorMaterial.SetShaderParameter("craters", _craters);
        }

        private float GetTerrainHeight(Vector2 posXZ)
        {
            float displacement = 0.0f;
            for (int i = 0; i < 16; i++)
            {
                Vector4 c = _craters[i];
                if (c.W <= 0.0001f) continue;
                float dist = posXZ.DistanceTo(new Vector2(c.X, c.Z));
                float radius = c.W;
                if (dist < radius * 2.0f)
                {
                    float influence = 1.0f - Mathf.SmoothStep(0.0f, radius, dist);
                    displacement -= influence * radius * 0.4f;
                    float rimInfluence = 1.0f - Mathf.SmoothStep(radius * 0.7f, radius * 1.5f, dist);
                    displacement += rimInfluence * influence * radius * 0.15f;
                }
            }
            return displacement;
        }

        private void ApplyFlightDynamics(float frameDt)
        {
            int subSteps = 4;
            float dt = frameDt / subSteps;

            Vector3 frameWindForce = _externalForceAccum;
            Vector3 frameWindTorque = _externalTorqueAccum;
            _lastExtForce = frameWindForce;
            _lastExtTorque = frameWindTorque;

            int intactCount = 0;
            List<int> activeIndices = new List<int>();
            for (int i = 0; i < 4; i++)
            {
                if (_rotorStructurallyIntact[i] && _motorActive[i])
                {
                    intactCount++;
                    activeIndices.Add(i);
                }
            }

            Vector3 frameMagnetForceAccum = Vector3.Zero;

            for (int step = 0; step < subSteps; step++)
            {
                Basis currentBasis = new Basis(_simState.Orientation);
                Vector3 stepForceAccum = frameWindForce;
                Vector3 stepTorqueAccum = frameWindTorque;

                // =========================================================
                // PHASE 1 & 4: Sub-Step Cable Integration & Critical Damping
                // =========================================================
                // Extrapolate magnet state into the integration sub-step
                Vector3 magnetPos = _magnetBody.GlobalPosition + _magnetBody.LinearVelocity * (step * dt);
                Vector3 magnetVel = _magnetBody.LinearVelocity;

                Vector3 globalHookRoot = _simState.Position + (currentBasis * _hookLocalOffset);
                Vector3 hookVel = _simState.Velocity + (currentBasis * _simState.AngularVelocity).Cross(currentBasis * _hookLocalOffset);

                Vector3 cableVec = magnetPos - globalHookRoot;
                float dist = cableVec.Length();
                Vector3 springForce = Vector3.Zero;

                if (dist > _currentCableLength)
                {
                    Vector3 cableDir = cableVec / dist;
                    float extension = dist - _currentCableLength;
                    Vector3 relVel = magnetVel - hookVel;

                    // Dynamic Critical Damping Formulation
                    float m_drone = _rk4.Mass;
                    float m_load = _magnetBody.Mass + (_hookedPayload != null && IsInstanceValid(_hookedPayload) ? _hookedPayload.Mass : 0f);
                    float m_eff = (m_drone * m_load) / (m_drone + m_load);
                    float c_crit = 2.0f * Mathf.Sqrt(CableSpringStiffness * m_eff);
                    float dampingCoef = c_crit * CableDampingRatio;

                    float dampingForce = relVel.Dot(cableDir) * dampingCoef;
                    float springForceMag = (extension * CableSpringStiffness) + dampingForce;
                    springForceMag = Mathf.Clamp(springForceMag, 0.0f, 15000.0f); // Prevents impulse blowout

                    springForce = cableDir * springForceMag;

                    stepForceAccum += springForce;
                    stepTorqueAccum += (currentBasis * _hookLocalOffset).Cross(springForce);

                    // Accumulate equal/opposite reaction to apply to Godot's physics solver at the end
                    frameMagnetForceAccum -= springForce;
                }

                if (step == subSteps - 1)
                {
                    _lastWinchDist = dist;
                    _lastWinchForce = springForce;
                }

                // =========================================================
                // PHASE 2: Feed-Forward Control Augmentation
                // =========================================================
                float feedForwardDownwardThrust = Mathf.Max(0f, -springForce.Y);

                Vector3 localAngVel = _simState.AngularVelocity;
                Vector3 globalGravity = new Vector3(0, -9.81f * _rk4.Mass, 0);

                Vector3 torqueCmd = _attController.ComputeTorque(_simState.Orientation, _targetAttitude, localAngVel);
                if (step == subSteps - 1) _lastTorqueCmd = torqueCmd;

                float analyticalGroundY = GetTerrainHeight(new Vector2(_simState.Position.X, _simState.Position.Z));
                float distToGround = _simState.Position.Y - analyticalGroundY;

                var spaceState = GetWorld3D().DirectSpaceState;
                var query = PhysicsRayQueryParameters3D.Create(_simState.Position, _simState.Position + Vector3.Down * 10.0f);
                var result = spaceState.IntersectRay(query);
                if (result.Count > 0)
                {
                    float physicalHitY = (float)result["position"].AsVector3().Y;
                    distToGround = Mathf.Min(distToGround, _simState.Position.Y - physicalHitY);
                }

                float groundEffect = 1.0f;
                if (distToGround < 10.0f && distToGround > 0f) groundEffect += 0.5f * Mathf.Exp(-distToGround * 2.0f);

                float tiltFactor = currentBasis.Y.Y;
                float clampedTilt = Mathf.Max(0.5f, tiltFactor);

                bool motorsSaturated = false;
                for (int i = 0; i < 4; i++) if (_actualMotorThrust[i] >= MaxMotorThrust * 0.95f) motorsSaturated = true;

                float altError = _targetAlt - _simState.Position.Y;
                bool freezeIntegral = motorsSaturated && altError > 0;

                float altCmd = 0f;
                if (tiltFactor > 0.4f) {
                    altCmd = _altPID.Update(altError, _simState.Velocity.Y, dt, freezeIntegral);
                } else {
                    _altPID.Bleed(5.0f, dt);
                }
                if (step == subSteps - 1) _lastPidOutput = altCmd;

                float hoverThrust = (_rk4.Mass * 9.81f) / (intactCount > 0 ? (intactCount * groundEffect * clampedTilt) : 1f);

                // Inject feed forward tension compensation immediately to fight sag
                float ffThrust = feedForwardDownwardThrust / (intactCount > 0 ? (intactCount * clampedTilt) : 1f);
                float baseThrust = hoverThrust + ffThrust + altCmd;

                float totalThrustCmd = baseThrust * intactCount;

                float[] cmdThrust = new float[4];

                if (intactCount == 4)
                {
                    float pitchMix = torqueCmd.X / (4.0f * ArmLength);
                    float rollMix = torqueCmd.Z / (4.0f * ArmLength);
                    float yawMix = torqueCmd.Y / (4.0f * YawDragFactor);

                    float yawAuthorityLimit = MaxMotorThrust * 0.35f;
                    yawMix = Mathf.Clamp(yawMix, -yawAuthorityLimit, yawAuthorityLimit);

                    cmdThrust[0] = baseThrust + pitchMix - rollMix + yawMix; // FL
                    cmdThrust[1] = baseThrust + pitchMix + rollMix - yawMix; // FR
                    cmdThrust[2] = baseThrust - pitchMix - rollMix - yawMix; // BL
                    cmdThrust[3] = baseThrust - pitchMix + rollMix + yawMix; // BR
                }
                else if (intactCount == 3)
                {
                    int i1 = activeIndices[0], i2 = activeIndices[1], i3 = activeIndices[2];

                    Vector3 col1 = new Vector3(1.0f, -_rotorPositions[i1].Z, _rotorPositions[i1].X);
                    Vector3 col2 = new Vector3(1.0f, -_rotorPositions[i2].Z, _rotorPositions[i2].X);
                    Vector3 col3 = new Vector3(1.0f, -_rotorPositions[i3].Z, _rotorPositions[i3].X);

                    Basis allocator = new Basis(col1, col2, col3);

                    if (Mathf.Abs(allocator.Determinant()) > 0.001f)
                    {
                        Basis invAllocator = allocator.Inverse();
                        Vector3 command = new Vector3(totalThrustCmd, torqueCmd.X, torqueCmd.Z);

                        Vector3 solved = invAllocator * command;
                        cmdThrust[i1] = solved.X;
                        cmdThrust[i2] = solved.Y;
                        cmdThrust[i3] = solved.Z;
                    }
                    else
                    {
                        for (int i = 0; i < 4; i++) cmdThrust[i] = activeIndices.Contains(i) ? baseThrust : 0f;
                    }
                }
                else
                {
                    for (int i = 0; i < 4; i++) cmdThrust[i] = activeIndices.Contains(i) ? baseThrust : 0f;
                }

                float maxRequested = -float.MaxValue;
                float minRequested = float.MaxValue;
                for (int i = 0; i < 4; i++)
                {
                    if (activeIndices.Contains(i)) {
                        maxRequested = Mathf.Max(maxRequested, cmdThrust[i]);
                        minRequested = Mathf.Min(minRequested, cmdThrust[i]);
                    }
                }

                if (maxRequested > MaxMotorThrust)
                {
                    float shiftDown = maxRequested - MaxMotorThrust;
                    for (int i = 0; i < 4; i++) if (activeIndices.Contains(i)) cmdThrust[i] -= shiftDown;
                }
                if (minRequested < 0)
                {
                    float shiftUp = -minRequested;
                    for (int i = 0; i < 4; i++) if (activeIndices.Contains(i)) cmdThrust[i] += shiftUp;
                }

                for (int i = 0; i < 4; i++)
                {
                    float targetThrust = (activeIndices.Contains(i)) ? Mathf.Clamp(cmdThrust[i], 0, MaxMotorThrust) : 0f;
                    float dThrust = (targetThrust - _actualMotorThrust[i]) * MotorTimeConstant;
                    _actualMotorThrust[i] = Mathf.Clamp(_actualMotorThrust[i] + dThrust * dt, 0f, MaxMotorThrust);
                }

                float totalThrust = (_actualMotorThrust[0] + _actualMotorThrust[1] + _actualMotorThrust[2] + _actualMotorThrust[3]) * _thrustDirection;
                Vector3 thrustForceLocal = new Vector3(0, totalThrust, 0);

                if (_simState.Velocity.Y < -2.0f && Mathf.Abs(_simState.Velocity.X) < 1.5f && Mathf.Abs(_simState.Velocity.Z) < 1.5f)
                {
                    float verticalSinkRate = Mathf.Abs(_simState.Velocity.Y);
                    float vrsFactor = Mathf.Clamp((verticalSinkRate - 2.0f) / 6.0f, 0.0f, 1.0f);
                    if (vrsFactor > 0.1f)
                    {
                        thrustForceLocal.Y *= (1.0f - (0.3f * vrsFactor));
                        stepTorqueAccum += currentBasis * new Vector3((GD.Randf() - 0.5f) * 4.0f * vrsFactor, 0f, (GD.Randf() - 0.5f) * 4.0f * vrsFactor);
                    }
                }

                Vector3 aeroTorqueLocal = new Vector3(
                    (_actualMotorThrust[0] + _actualMotorThrust[1] - _actualMotorThrust[2] - _actualMotorThrust[3]) * ArmLength,
                    (_actualMotorThrust[0] + _actualMotorThrust[3] - _actualMotorThrust[1] - _actualMotorThrust[2]) * YawDragFactor * _thrustDirection,
                    (_actualMotorThrust[1] + _actualMotorThrust[3] - _actualMotorThrust[0] - _actualMotorThrust[2]) * ArmLength
                );

                Vector3 totalRotorDragGlobal = Vector3.Zero;
                Vector3 totalRotorDragTorqueLocal = Vector3.Zero;
                float totalLocalLy = 0f;
                float rotorInertia = 0.5f * RotorMass * (RotorRadius * RotorRadius);

                float groundStiffness = 5000.0f;
                float groundDamping = 200.0f;

                for (int i = 0; i < 4; i++)
                {
                    if (!_rotorStructurallyIntact[i]) continue;

                    float spinDirection = (i == 0 || i == 3) ? -1.0f : 1.0f;
                    totalLocalLy += rotorInertia * spinDirection * Mathf.Sqrt(_actualMotorThrust[i]) * RpmScaleFactor;

                    Vector3 rotorVelGlobal = _simState.Velocity + (currentBasis * localAngVel).Cross(currentBasis * _rotorPositions[i]);
                    Vector3 localRotorVel = currentBasis.Inverse() * rotorVelGlobal;

                    Vector3 rotorDragLocal = new Vector3(
                        -localRotorVel.X * Mathf.Abs(localRotorVel.X) * 0.15f,
                        -localRotorVel.Y * Mathf.Abs(localRotorVel.Y) * 0.05f,
                        -localRotorVel.Z * Mathf.Abs(localRotorVel.Z) * 0.15f
                    );

                    totalRotorDragGlobal += currentBasis * rotorDragLocal;
                    totalRotorDragTorqueLocal += _rotorPositions[i].Cross(rotorDragLocal);

                    Vector3 globalTip = _simState.Position + (currentBasis * _rotorPositions[i]);
                    float floorY = GetTerrainHeight(new Vector2(globalTip.X, globalTip.Z));

                    if (globalTip.Y < floorY)
                    {
                        float penetration = floorY - globalTip.Y;
                        float dampingForce = -rotorVelGlobal.Y * groundDamping;
                        float upwardForce = Mathf.Max(0, (penetration * groundStiffness) + dampingForce);

                        Vector3 lateralVel = new Vector3(rotorVelGlobal.X, 0, rotorVelGlobal.Z);
                        Vector3 frictionForce = -lateralVel * (upwardForce * 0.8f);

                        float maxFriction = (lateralVel.Length() * _rk4.Mass) / (4.0f * dt);
                        if (frictionForce.LengthSquared() > maxFriction * maxFriction) {
                            frictionForce = frictionForce.Normalized() * maxFriction;
                        }

                        Vector3 totalTipForce = new Vector3(0, upwardForce, 0) + frictionForce;
                        stepForceAccum += totalTipForce;
                        stepTorqueAccum += (currentBasis * _rotorPositions[i]).Cross(totalTipForce);

                        _simState.AngularVelocity = _simState.AngularVelocity.Lerp(Vector3.Zero, 4.0f * dt);

                        if (step == subSteps - 1)
                        {
                            if (upwardForce > MinimumCraterForce) RegisterCrater(globalTip, upwardForce);
                            if (rotorVelGlobal.LengthSquared() > 4.0f && !_sparkEmitter.Emitting)
                            {
                              _sparkEmitter.GlobalPosition = globalTip;
                              _sparkEmitter.Restart();
                              _sparkEmitter.Emitting = true;
                            }

                            if (RotorBreakageMultiplier > 0.0f)
                            {
                                float scaledThreshold = StructuralBreakForce * RotorBreakageMultiplier;
                                float structuralForceFactor = upwardForce * Mathf.Max(0.1f, Mathf.Abs(rotorVelGlobal.Normalized().Y));

                                if (structuralForceFactor > scaledThreshold)
                                {
                                    HandleComponentFailure(i);
                                }
                            }
                        }
                    }
                }

                float halfChassisHeight = 0.1f;
                if (_simState.Position.Y - halfChassisHeight < analyticalGroundY)
                {
                    float penetration = analyticalGroundY - (_simState.Position.Y - halfChassisHeight);
                    float dampingForce = -_simState.Velocity.Y * (groundDamping * 2.0f);
                    float upwardForce = Mathf.Max(0, (penetration * groundStiffness * 2.0f) + dampingForce);

                    stepForceAccum += new Vector3(0, upwardForce, 0);

                    if (step == subSteps - 1 && upwardForce > MinimumCraterForce)
                        RegisterCrater(_simState.Position - new Vector3(0, halfChassisHeight, 0), upwardForce);
                }

                Vector3 angularDragLocal = -localAngVel * localAngVel.Length() * 0.08f;
                Vector3 totalGlobalForce = globalGravity + totalRotorDragGlobal + (currentBasis * thrustForceLocal) + stepForceAccum;
                Vector3 totalLocalTorque = aeroTorqueLocal + totalRotorDragTorqueLocal + angularDragLocal + (currentBasis.Inverse() * stepTorqueAccum);

                _simState = _rk4.Integrate(_simState, dt, totalGlobalForce, totalLocalTorque);

                if (step == subSteps - 1)
                {
                    _telemetryLy = totalLocalLy;
                    _telemetryGyroTorque = _rk4.InertiaTensor * _simState.AngularVelocity;
                }
            }

            // Apply the averaged sub-step tension back to Godot's Jolt/PhysicsServer magnet body
            _magnetBody.ApplyCentralForce(frameMagnetForceAccum / subSteps);

            _drone.GlobalTransform = new Transform3D(new Basis(_simState.Orientation), _simState.Position);
        }

        private void HandleComponentFailure(int index)
        {
          _rotorStructurallyIntact[index] = false;
          _motorActive[index] = false;

          float massLost = RotorMass;
          Vector3 relativePos = _rotorPositions[index];

          Vector3 comShiftLocal = (-massLost * relativePos) / (_rk4.Mass - massLost);
          Vector3 comShiftGlobal = new Basis(_simState.Orientation) * comShiftLocal;

          _simState.Position += comShiftGlobal;

          for(int i = 0; i < 4; i++)
          {
            if (_rotorStructurallyIntact[i]) _rotorPositions[i] -= comShiftLocal;
          }

          // PHASE 3: Make torque mechanically CoM-Aware
          _hookLocalOffset -= comShiftLocal;

          _rk4.RecalculateInertia(massLost, relativePos);
          _rotors[index].Visible = false;

          RigidBody3D debris = new RigidBody3D { };
          debris.CollisionLayer = 0;
          debris.CollisionMask = 0;

          MeshInstance3D debrisMesh = new MeshInstance3D { Mesh = new CylinderMesh { TopRadius = RotorRadius, BottomRadius = RotorRadius, Height = 0.02f } };
          debris.AddChild(debrisMesh);

          var debrisSmoke = new GpuParticles3D
          {
            ProcessMaterial = _smokeEmitters[index].ProcessMaterial,
            DrawPass1 = _smokeEmitters[index].DrawPass1,
            Amount = 16,
            Emitting = true,
            Lifetime = 1.0f
          };
          debris.AddChild(debrisSmoke);

          GetParent().AddChild(debris);

          debris.LinearVelocity = _simState.Velocity + (new Basis(_simState.Orientation) * _simState.AngularVelocity).Cross(new Basis(_simState.Orientation) * relativePos);
          debris.AngularVelocity = new Basis(_simState.Orientation) * _simState.AngularVelocity + new Vector3(GD.Randf() - 0.5f, GD.Randf() - 0.5f, GD.Randf() - 0.5f) * 15.0f;

          _activeDebris.Add(debris);
          debris.GlobalPosition = _drone.GlobalPosition + comShiftGlobal + relativePos;
        }

        private void SetupDownwashVolume()
        {
            _downwashVolume = new Area3D();
            CollisionShape3D coneCol = new CollisionShape3D { Shape = new CylinderShape3D { Radius = 5.0f, Height = 10.0f }, Position = new Vector3(0, -5.0f, 0) };
            _downwashVolume.AddChild(coneCol);
            _drone.AddChild(_downwashVolume);

            _downwashVolume.BodyEntered += (Node3D body) => { if (body is RigidBody3D rb) _bodiesInDownwash.Add(rb); };
            _downwashVolume.BodyExited += (Node3D body) => { if (body is RigidBody3D rb) _bodiesInDownwash.Remove(rb); };
        }

        private void ApplyRotorDownwash(float dt)
        {
            float totalThrust = _actualMotorThrust[0] + _actualMotorThrust[1] + _actualMotorThrust[2] + _actualMotorThrust[3];
            if (totalThrust < 0.1f) return;

            float diskArea = 4.0f * Mathf.Pi * (RotorRadius * RotorRadius);
            float inducedVelocity = Mathf.Sqrt(totalThrust / (2.0f * 1.225f * diskArea));

            Basis droneBasis = _drone.GlobalTransform.Basis;
            Vector3 downVector = -droneBasis.Y;

            foreach (var rb in _bodiesInDownwash)
            {
                Vector3 relativePos = rb.GlobalPosition - _drone.GlobalPosition;
                float verticalDist = relativePos.Dot(downVector);
                if (verticalDist <= 0) continue;

                Vector3 lateralPos = relativePos - (downVector * verticalDist);
                float radialDist = lateralPos.Length();

                float wakeRadius = RotorRadius * 2.0f + (0.1f * verticalDist);
                if (radialDist > wakeRadius * 2.0f) continue;

                float localVelocity = inducedVelocity * (RotorRadius / wakeRadius) * Mathf.Exp(-Mathf.Pow(radialDist / wakeRadius, 2.0f));

                float dragForceMag = 0.5f * 1.225f * 1.0f * 0.25f * (localVelocity * localVelocity);
                Vector3 forceVector = downVector + (lateralPos.Normalized() * 0.2f);

                rb.ApplyCentralForce(forceVector.Normalized() * dragForceMag);
            }
        }

        private void CompileShaders()
        {
            _plumeShader = new Shader { Code = FLUID_PLUME_SHADER };
            _portalShader = new Shader { Code = @"
                shader_type spatial; render_mode unshaded, cull_disabled;
                uniform vec4 portal_color : source_color = vec4(0.0, 0.5, 1.0, 1.0);
                void fragment() {
                    vec2 uv = UV * 2.0 - 1.0; float dist = length(uv); if (dist > 1.0) discard;
                    float ripple = sin(dist * 30.0 - TIME * 15.0) * 0.5 + 0.5;
                    ALBEDO = portal_color.rgb * (ripple * smoothstep(0.0, 0.8, dist) + smoothstep(0.8, 1.0, dist) * 2.0); ALPHA = 1.0;
                }" };

            _floorShader = new Shader { Code = @"
                shader_type spatial;
                uniform vec4 craters[16];
                varying vec3 v_world_pos;

                void vertex() {
                    v_world_pos = (MODEL_MATRIX * vec4(VERTEX, 1.0)).xyz;
                    float displacement = 0.0;

                    for(int i = 0; i < 16; i++) {
                        float dist = length(v_world_pos.xz - craters[i].xz);
                        float radius = craters[i].w;

                        if (radius > 0.0 && dist < radius * 2.0) {
                            float influence = smoothstep(radius, 0.0, dist);
                            displacement -= influence * radius * 0.4;
                            displacement += smoothstep(radius * 1.5, radius * 0.7, dist) * influence * radius * 0.15;
                        }
                    }
                    VERTEX.y += displacement;
                }

                void fragment() {
                    vec2 grid = floor(v_world_pos.xz * 2.0);
                    ALBEDO = mix(vec3(0.15), vec3(0.25), mod(grid.x + grid.y, 2.0));
                }" };

                _smokeShader = new Shader { Code = @"
                  shader_type spatial;
                  render_mode blend_mix, depth_draw_never, cull_disabled;

                  uniform sampler2D depth_texture : hint_depth_texture, repeat_disable, filter_nearest;

                  float hash(vec2 p) { return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453); }
                  float noise(vec2 p) {
                    vec2 i = floor(p); vec2 f = fract(p); f = f * f * (3.0 - 2.0 * f);
                    return mix(mix(hash(i), hash(i + vec2(1.0, 0.0)), f.x),
                        mix(hash(i + vec2(0.0, 1.0)), hash(i + vec2(1.0, 1.0)), f.x), f.y);
                  }

                  void vertex() {
                    MODELVIEW_MATRIX = VIEW_MATRIX * mat4(INV_VIEW_MATRIX[0], INV_VIEW_MATRIX[1], INV_VIEW_MATRIX[2], MODEL_MATRIX[3]);
                    MODELVIEW_NORMAL_MATRIX = mat3(MODELVIEW_MATRIX);
                  }

                  void fragment() {
                    vec2 uv = UV * 2.0 - 1.0;
                    float dist = length(uv);
                    if (dist > 1.0) discard;

                    float n = noise(uv * 2.5 + vec2(COLOR.r * 15.0));
                    float alpha_mask = smoothstep(1.0, 0.2, dist) * n;

                    float depth = texture(depth_texture, SCREEN_UV).x;
                    vec3 ndc = vec3(SCREEN_UV * 2.0 - 1.0, depth);
                    vec4 view_pos = INV_PROJECTION_MATRIX * vec4(ndc, 1.0);
                    float depth_z = view_pos.z / view_pos.w;
                    float proximity_fade = smoothstep(0.0, 1.5, VERTEX.z - depth_z);

                    vec3 smoke_color = mix(vec3(0.05), vec3(0.35), n);
                    ALBEDO = smoke_color * COLOR.rgb;
                    ALPHA = alpha_mask * COLOR.a * proximity_fade * 0.8;
                  }"
                };
        }

        private void BuildPortals()
        {
            BuildPortalPair(new Vector3(15, 20, 0), Vector3.Down, Colors.DeepSkyBlue, new Vector3(0, 10, -30), Vector3.Back, Colors.Orange);
            BuildPortalPair(new Vector3(-20, 0.1f, -10), Vector3.Up, Colors.LimeGreen, new Vector3(25, 35, -20), Vector3.Right, Colors.MediumPurple);
            BuildPortalPair(new Vector3(-35, 15, 20), Vector3.Right, Colors.Cyan, new Vector3(10, 40, 20), Vector3.Down, Colors.Crimson);
        }

        private void BuildPortalPair(Vector3 posA, Vector3 fwdA, Color colA, Vector3 posB, Vector3 fwdB, Color colB)
        {
            var pA = CreatePortal(posA, fwdA, colA); var pB = CreatePortal(posB, fwdB, colB);
            ConnectPortals(pA, pB); ConnectPortals(pB, pA);
        }

        private Area3D CreatePortal(Vector3 pos, Vector3 forward, Color color)
        {
            var p = new Area3D();
            AddChild(p);
            p.GlobalPosition = pos;
            p.GlobalBasis = Basis.LookingAt(forward, Math.Abs(forward.Dot(Vector3.Up)) > 0.99f ? Vector3.Right : Vector3.Up);
            var mat = new ShaderMaterial { Shader = _portalShader }; mat.SetShaderParameter("portal_color", color);
            p.AddChild(new MeshInstance3D { Mesh = new PlaneMesh { Size = new Vector2(6, 6) }, MaterialOverride = mat, RotationDegrees = new Vector3(90, 0, 0) });
            p.AddChild(new CollisionShape3D { Shape = new CylinderShape3D { Radius = 3.0f, Height = 0.5f }, RotationDegrees = new Vector3(90, 0, 0) });
            return p;
        }

        private void ConnectPortals(Area3D entry, Area3D exit)
        {
            entry.BodyEntered += (Node3D body) => {
                if (body.Name == "DroneBody") {
                    ulong now = Time.GetTicksMsec();
                    ulong last = body.HasMeta("portal_time") ? (ulong)body.GetMeta("portal_time") : 0ul;
                    if (now - last < 300) return;
                    body.SetMeta("portal_time", now);
                    TeleportEntity(entry, exit);
                }
            };
        }

        private void TeleportEntity(Area3D inPortal, Area3D outPortal)
        {
            Transform3D relativeTrans = inPortal.GlobalTransform.AffineInverse() * new Transform3D(new Basis(_simState.Orientation), _simState.Position);
            Transform3D flipY = new Transform3D(new Basis(Vector3.Up, Mathf.Pi), Vector3.Zero);
            Transform3D newTrans = outPortal.GlobalTransform * flipY * relativeTrans;
            newTrans.Origin += -outPortal.GlobalBasis.Z * 1.5f;

            Vector3 localLinVel = inPortal.GlobalBasis.Inverse() * _simState.Velocity;
            Vector3 localAngVel = inPortal.GlobalBasis.Inverse() * (new Basis(_simState.Orientation) * _simState.AngularVelocity);

            _simState.Position = newTrans.Origin;
            _simState.Orientation = newTrans.Basis.GetRotationQuaternion();
            _simState.Velocity = outPortal.GlobalBasis * (flipY.Basis * localLinVel);

            Vector3 newGlobalAngVel = outPortal.GlobalBasis * (flipY.Basis * localAngVel);
            _simState.AngularVelocity = newTrans.Basis.Inverse() * newGlobalAngVel;

            _altPID.Reset();
            _targetAlt = _simState.Position.Y;

            Vector3 forward = newTrans.Basis.Z;
            _camYaw = Mathf.Atan2(forward.X, forward.Z);

            _targetYawAccumulator = _camYaw;
            _targetAttitude = _simState.Orientation;

            Vector3 newForward = new Basis(_simState.Orientation).Z;
            _lastDroneYaw = Mathf.Atan2(newForward.X, newForward.Z);

            // 1. Teleport the Winch Magnet
            Transform3D relMagnetTrans = inPortal.GlobalTransform.AffineInverse() * _magnetBody.GlobalTransform;
            Transform3D newMagnetTrans = outPortal.GlobalTransform * flipY * relMagnetTrans;
            newMagnetTrans.Origin += -outPortal.GlobalBasis.Z * 1.5f; // Apply the same exit offset

            Vector3 localMagnetVel = inPortal.GlobalBasis.Inverse() * _magnetBody.LinearVelocity;
            Vector3 localMagnetAngVel = inPortal.GlobalBasis.Inverse() * _magnetBody.AngularVelocity;

            _magnetBody.GlobalTransform = newMagnetTrans;
            _magnetBody.LinearVelocity = outPortal.GlobalBasis * (flipY.Basis * localMagnetVel);
            _magnetBody.AngularVelocity = outPortal.GlobalBasis * (flipY.Basis * localMagnetAngVel);

            // 2. Teleport the Payload (if one is currently attached)
            if (_hookedPayload != null && IsInstanceValid(_hookedPayload))
            {
              Transform3D relPayloadTrans = inPortal.GlobalTransform.AffineInverse() * _hookedPayload.GlobalTransform;
              Transform3D newPayloadTrans = outPortal.GlobalTransform * flipY * relPayloadTrans;
              newPayloadTrans.Origin += -outPortal.GlobalBasis.Z * 1.5f;

              Vector3 localPayloadVel = inPortal.GlobalBasis.Inverse() * _hookedPayload.LinearVelocity;
              Vector3 localPayloadAngVel = inPortal.GlobalBasis.Inverse() * _hookedPayload.AngularVelocity;

              _hookedPayload.GlobalTransform = newPayloadTrans;
              _hookedPayload.LinearVelocity = outPortal.GlobalBasis * (flipY.Basis * localPayloadVel);
              _hookedPayload.AngularVelocity = outPortal.GlobalBasis * (flipY.Basis * localPayloadAngVel);
            }
        }

        private void BuildDrone()
        {
          _drone = new AnimatableBody3D { Name = "DroneBody", Position = _simState.Position, SyncToPhysics = false };
          AddChild(_drone);

          var shinySilverMat = new StandardMaterial3D
          {
            AlbedoColor = new Color(0.85f, 0.85f, 0.88f),
            Metallic = 1.0f,
            Roughness = 0.05f
          };

          var payloadMesh = new MeshInstance3D
          {
            Mesh = new SphereMesh { Radius = 0.35f, Height = 0.7f },
            MaterialOverride = shinySilverMat,
            Scale = new Vector3(1.0f, 0.5f, 1.0f)
          };
          _drone.AddChild(payloadMesh);

          _drone.AddChild(new CollisionShape3D {
              Shape = new CylinderShape3D { Radius = 0.35f, Height = 0.35f }
              });

          _rotorPositions[0] = new Vector3(-ArmLength, 0, -ArmLength);
          _rotorPositions[1] = new Vector3(ArmLength, 0, -ArmLength);
          _rotorPositions[2] = new Vector3(-ArmLength, 0, ArmLength);
          _rotorPositions[3] = new Vector3(ArmLength, 0, ArmLength);

          var smokeScaleCurve = new Curve();
          smokeScaleCurve.AddPoint(new Vector2(0.0f, 0.2f));
          smokeScaleCurve.AddPoint(new Vector2(1.0f, 2.5f));

          var smokeGradient = new Gradient();
          smokeGradient.SetColors(new[] { new Color(1,1,1, 0.8f), new Color(0.3f,0.3f,0.3f, 0.5f), new Color(0,0,0, 0.0f) });
          smokeGradient.SetOffsets(new[] { 0.0f, 0.3f, 1.0f });

          var smokeProcessMat = new ParticleProcessMaterial {
            EmissionShape = ParticleProcessMaterial.EmissionShapeEnum.Sphere,
            EmissionSphereRadius = 0.05f,
            Direction = new Vector3(0, 1, 0),
            Spread = 35.0f,
            InitialVelocityMin = 1.0f,
            InitialVelocityMax = 2.5f,
            Gravity = new Vector3(0, 3.0f, 0),
            AngleMin = 0.0f,
            AngleMax = 360.0f,
            ScaleCurve = new CurveTexture { Curve = smokeScaleCurve },
            ColorRamp = new GradientTexture1D { Gradient = smokeGradient }
          };

          AudioStream motorStreamAsset = null;
          if (ResourceLoader.Exists("res://audio/motor_loop.wav")) {
            motorStreamAsset = GD.Load<AudioStream>("res://audio/motor_loop.wav");
          }
          var smokeMaterial = new ShaderMaterial { Shader = _smokeShader };
          var smokeDrawPass = new QuadMesh { Material = smokeMaterial, Size = new Vector2(0.2f, 0.2f) };

          for (int i = 0; i < 4; i++)
          {
            float exactYawAngle = Mathf.RadToDeg(Mathf.Atan2(_rotorPositions[i].X, _rotorPositions[i].Z));

            var arm = new MeshInstance3D
            {
              Mesh = new CylinderMesh { TopRadius = 0.04f, BottomRadius = 0.04f, Height = ArmLength * Mathf.Sqrt(2.0f),
              },
              MaterialOverride = shinySilverMat,
              RotationDegrees = new Vector3(90, exactYawAngle, 0),
              Position = _rotorPositions[i] / 2.0f
            };
            _drone.AddChild(arm);

            var rotor = new MeshInstance3D { Mesh = new CylinderMesh { TopRadius = RotorRadius, BottomRadius = RotorRadius, Height = 0.02f }, MaterialOverride = new StandardMaterial3D { AlbedoColor = i < 2 ? new Color(1, 0, 0, 0.5f) : new Color(0, 0, 0, 0.5f), Transparency = BaseMaterial3D.TransparencyEnum.Alpha }, Position = _rotorPositions[i] + new Vector3(0, 0.1f, 0) };
            _drone.AddChild(rotor); _rotors[i] = rotor;

                var plumeMat = new ShaderMaterial { Shader = _plumeShader };

                var vectorMesh = new MeshInstance3D
                {
                  Mesh = new CylinderMesh
                  {
                    TopRadius = 0.03f,
                    BottomRadius = 0.04f,
                    Height = 0.6f,
                    CapTop = false,
                    CapBottom = false
                  },
                  MaterialOverride = plumeMat,
                  Position = new Vector3(0, -0.5f, 0)
                };

                var vectorPivot = new Node3D { Position = _rotorPositions[i] };
                vectorPivot.AddChild(vectorMesh);
                _drone.AddChild(vectorPivot);
                _thrustVectors[i] = vectorMesh;

                _smokeEmitters[i] = new GpuParticles3D {
                 ProcessMaterial = smokeProcessMat,
                 DrawPass1 = smokeDrawPass,
                 Position = _rotorPositions[i], Emitting = false, Amount = 16
                };
                _drone.AddChild(_smokeEmitters[i]);

                if (motorStreamAsset != null)
                {
                    _motorAudio[i] = new AudioStreamPlayer3D { Stream = motorStreamAsset, Autoplay = true, VolumeDb = -80f };
                    _drone.AddChild(_motorAudio[i]);
                }
            }

          var sparkScaleCurve = new Curve();
          sparkScaleCurve.AddPoint(new Vector2(0.0f, 1.0f));
          sparkScaleCurve.AddPoint(new Vector2(1.0f, 0.0f));

          var sparkGradient = new Gradient();
          sparkGradient.SetColors(new[] { Colors.White, Colors.Yellow, Colors.OrangeRed, new Color(0, 0, 0, 0) });
          sparkGradient.SetOffsets(new[] { 0.0f, 0.1f, 0.6f, 1.0f });

          var sparkProcessMat = new ParticleProcessMaterial {
            Direction = new Vector3(0, 1, 0),
            Spread = 60,
            InitialVelocityMin = 5f,
            InitialVelocityMax = 12f,
            Gravity = new Vector3(0, -15.0f, 0),
            ScaleCurve = new CurveTexture { Curve = sparkScaleCurve },
            ColorRamp = new GradientTexture1D { Gradient = sparkGradient },
          };

          var sparkVisualMat = new StandardMaterial3D {
            Transparency = BaseMaterial3D.TransparencyEnum.Alpha,
            VertexColorUseAsAlbedo = true,
            EmissionEnabled = true,
            Emission = Colors.White,
            EmissionEnergyMultiplier = 4.0f,
            BillboardMode = BaseMaterial3D.BillboardModeEnum.Particles
          };

          _sparkEmitter = new GpuParticles3D {
            ProcessMaterial = sparkProcessMat,
            DrawPass1 = new QuadMesh { Material = sparkVisualMat, Size = new Vector2(0.04f, 0.04f) },
            Emitting = true,
            OneShot = true,
            Explosiveness = 0.95f,
            Amount = 40,
            Lifetime = 0.6f,
            LocalCoords = false
          };
          AddChild(_sparkEmitter);
        }

        private void UpdateVisualsAndAudio(float dt)
        {
            // --- Update Winch Visuals here using latest integration states ---
            Vector3 globalHookRoot = _drone.GlobalTransform * _hookLocalOffset;
            Vector3 magnetPos = _magnetBody.GlobalPosition;
            Vector3 cableVec = magnetPos - globalHookRoot;
            float dist = cableVec.Length();

            if (dist > 0.01f)
            {
                _cableVisual.Visible = true;
                _cableVisual.GlobalPosition = globalHookRoot + (cableVec * 0.5f);
                Vector3 cableDir = cableVec / dist;
                Vector3 upReference = Mathf.Abs(cableDir.Dot(Vector3.Up)) > 0.999f ? Vector3.Right : Vector3.Up;
                _cableVisual.LookAt(magnetPos, upReference);
                _cableVisual.RotateObjectLocal(Vector3.Right, Mathf.Pi / 2.0f);
                _cableVisual.Scale = new Vector3(1, dist, 1);
            }
            else
            {
                _cableVisual.Visible = false;
            }

            _magnetVisual.GlobalPosition = magnetPos;
            _magnetVisual.GlobalRotation = _magnetBody.GlobalRotation;
            // -------------------------------------------------------------

            for (int m = 0; m < 4; m++)
            {
                if (_rotors[m] != null && _rotorStructurallyIntact[m])
                    _rotors[m].RotateY(Mathf.Sqrt(_actualMotorThrust[m]) * 10.0f * ((m == 0 || m == 3) ? -1.0f : 1.0f) * _thrustDirection * dt);

                if (_thrustVectors[m] != null)
                {
                    float ts = _actualMotorThrust[m] / MaxMotorThrust;
                    _thrustVectors[m].Visible = ts > 0.01f && _rotorStructurallyIntact[m];
                    if (ts > 0.01f) { _thrustVectors[m].Scale = new Vector3(1.0f + ts, ts * 4.0f, 1.0f + ts); _thrustVectors[m].Position = new Vector3(0, -(ts * 4.0f) / 2.0f * _thrustDirection, 0); _thrustVectors[m].RotationDegrees = new Vector3(_thrustDirection < 0 ? 180 : 0, 0, 0); }
                }

                _smokeEmitters[m].Emitting = !_motorActive[m];

                if (_motorAudio[m] != null && _motorAudio[m].Playing)
                {
                    float normalizedThrust = _actualMotorThrust[m] / MaxMotorThrust;
                    _motorAudio[m].PitchScale = Mathf.Lerp(0.5f, 2.0f, normalizedThrust);
                    _motorAudio[m].VolumeDb = _motorActive[m] && _rotorStructurallyIntact[m] && normalizedThrust > 0.01f ? Mathf.LinearToDb(normalizedThrust * 0.5f) : -80f;
                }
            }

            if (_camera != null && _drone != null)
            {
              Vector3 droneUp = _drone.GlobalTransform.Basis.Y;
              float tiltAngle = Mathf.RadToDeg(Mathf.Acos(Mathf.Clamp(droneUp.Dot(Vector3.Up), -1.0f, 1.0f)));

              if (CameraFollowsDrone && tiltAngle < MaxStabilizationAngle)
              {
                Vector3 droneForward = new Basis(_simState.Orientation).Z;
                float currentDroneYaw = Mathf.Atan2(droneForward.X, droneForward.Z);

                float yawDelta = Mathf.AngleDifference(_lastDroneYaw, currentDroneYaw);
                _camYaw += yawDelta;

                _lastDroneYaw = currentDroneYaw;
              }
              else
              {
                Vector3 droneForward = new Basis(_simState.Orientation).Z;
                _lastDroneYaw = Mathf.Atan2(droneForward.X, droneForward.Z);
              }

              float upDot = _drone.GlobalTransform.Basis.Y.Dot(Vector3.Up);
              float stress = Mathf.Clamp(1.0f - upDot, 0.0f, 2.0f);

              Vector3 offset = new Vector3(Mathf.Cos(_camPitch) * Mathf.Sin(_camYaw), Mathf.Sin(_camPitch), Mathf.Cos(_camPitch) * Mathf.Cos(_camYaw)) * _camDistance;
              Vector3 targetCamPos = _drone.GlobalPosition + offset;

              if (stress > 0.15f)
              {
                float shakeAmt = stress * 0.15f;
                targetCamPos += new Vector3((GD.Randf() - 0.5f) * shakeAmt, (GD.Randf() - 0.5f) * shakeAmt, (GD.Randf() - 0.5f) * shakeAmt);
              }

              float lerpFactor = 1.0f - Mathf.Exp(-15.0f * dt);
              if (_camera.GlobalPosition.DistanceSquaredTo(targetCamPos) > 200.0f)
                _camera.GlobalPosition = targetCamPos;
              else
                _camera.GlobalPosition = _camera.GlobalPosition.Lerp(targetCamPos, lerpFactor);

              Vector3 camForward = (_drone.GlobalPosition - _camera.GlobalPosition).Normalized();
              Vector3 referenceUp = Vector3.Up;

              Vector3 camRight = referenceUp.Cross(camForward).Normalized();
              Vector3 camUp = camForward.Cross(camRight).Normalized();

              _camera.GlobalTransform = new Transform3D(new Basis(camRight, camUp, -camForward), _camera.GlobalPosition);
            }

            if (_camera != null)
            {
                float speed = _simState.Velocity.Length();
                float fovAlpha = Mathf.Clamp(speed / FovSpeedThreshold, 0f, 1f);
                fovAlpha *= fovAlpha;

                float targetFov = Mathf.Lerp(BaseFOV, MaxFOV, fovAlpha);
                _camera.Fov = Mathf.Lerp(_camera.Fov, targetFov, dt * FovInterpolationRate);
            }
        }

        private void SetupEnvironment()
        {
          _camera = new Camera3D { Current = true, Fov = BaseFOV }; AddChild(_camera);
          AddChild(new DirectionalLight3D { Position = new Vector3(5, 50, 10), ShadowEnabled = true, RotationDegrees = new Vector3(-45, 45, 0) });

          var environment = new Godot.Environment
          {
            BackgroundMode = Godot.Environment.BGMode.Sky,
            Sky = new Sky { SkyMaterial = new ProceduralSkyMaterial { SkyTopColor = new Color(0.35f, 0.55f, 0.85f), SkyHorizonColor = new Color(0.65f, 0.75f, 0.85f) } },
            AmbientLightSource = Godot.Environment.AmbientSource.Sky,
            AmbientLightColor = new Color(0.2f, 0.2f, 0.2f),
            AmbientLightSkyContribution = 0.5f
          };

          var worldEnv = new WorldEnvironment { Environment = environment };
          AddChild(worldEnv);

          float floorSize = 999.0f;
          int subdivisions = 800;

          var floor = new StaticBody3D();
          _floorMaterial = new ShaderMaterial { Shader = _floorShader };
          _floorMaterial.SetShaderParameter("craters", _craters);

          var terrainMesh = new PlaneMesh
          {
            Size = new Vector2(floorSize, floorSize),
            SubdivideWidth = subdivisions,
            SubdivideDepth = subdivisions
          };

          floor.AddChild(new MeshInstance3D { Mesh = terrainMesh, MaterialOverride = _floorMaterial });

          floor.AddChild(new CollisionShape3D
              {
              Shape = new BoxShape3D { Size = new Vector3(floorSize, 1.0f, floorSize) },
              Position = new Vector3(0, -0.5f, 0)
              });

          AddChild(floor);

          AddChild(new TurbulenceVolume { Position = new Vector3(15, 5, 0), VolumeSize = new Vector3(8, 20, 8), WindDirection = Vector3.Up, BaseWindForce = 35.0f, TurbulenceStrength = 15.0f, ZoneColor = new Color(1.0f, 0.5f, 0.0f, 0.2f) });
          AddChild(new TurbulenceVolume { Position = new Vector3(-15, 5, 0), VolumeSize = new Vector3(10, 10, 10), WindDirection = Vector3.Right, BaseWindForce = 25.0f, TurbulenceStrength = 20.0f, ZoneColor = new Color(0.0f, 0.5f, 1.0f, 0.2f) });
        }

        private void SetupHUD()
        {
          var canvas = new CanvasLayer();
          AddChild(canvas);

          _hud = new Label
          {
            LabelSettings = new LabelSettings
            {
              FontSize = 18,
              FontColor = Colors.Cyan,
              OutlineSize = 3,
              OutlineColor = Colors.Black
            },
            HorizontalAlignment = HorizontalAlignment.Right,
            VerticalAlignment = VerticalAlignment.Bottom
          };

          // Anchor to bottom-right of screen
          _hud.SetAnchorsAndOffsetsPreset(Control.LayoutPreset.BottomRight, Control.LayoutPresetMode.KeepSize, 20);
          // Shift margin inward from the screen edge
          _hud.Position = new Vector2(_hud.Position.X - 230, _hud.Position.Y - 250);

          canvas.AddChild(_hud);
        }

        private void UpdateHUD()
        {
          float verticalVelocity = _simState.Velocity.Y;

          string payloadInfo = _hookedPayload == null ? "EMPTY" : $"{_hookedPayload.Mass:F1}kg Attached";

          _hud.Text = $"TARGET ALT: {_targetAlt:F1}m\n" +
            $"ACTUAL ALT: {_simState.Position.Y:F1}m\n" +
            $"CAM FOLLOW (TAB): {(CameraFollowsDrone ? "ON" : "OFF")}\n\n" +
            $"[1] FL: {(_rotorStructurallyIntact[0] ? (_motorActive[0] ? "ON" : "FAIL") : "MISSING")} | Pwr: {_actualMotorThrust[0] * _thrustDirection:F1}\n" +
            $"[2] FR: {(_rotorStructurallyIntact[1] ? (_motorActive[1] ? "ON" : "FAIL") : "MISSING")} | Pwr: {_actualMotorThrust[1] * _thrustDirection:F1}\n" +
            $"[3] BL: {(_rotorStructurallyIntact[2] ? (_motorActive[2] ? "ON" : "FAIL") : "MISSING")} | Pwr: {_actualMotorThrust[2] * _thrustDirection:F1}\n" +
            $"[4] BR: {(_rotorStructurallyIntact[3] ? (_motorActive[3] ? "ON" : "FAIL") : "MISSING")} | Pwr: {_actualMotorThrust[3] * _thrustDirection:F1}\n\n";

          _hud.LabelSettings.FontColor = _thrustDirection > 0 ? Colors.Cyan : Colors.Red;
        }
    }
}
