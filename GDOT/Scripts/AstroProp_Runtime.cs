// Jonathan C. Ross 
// 09 01 2023 C


// Fuck you, Unity.

// Hello, welcome to the math 
// If you're coming from Unity, you will have a problem assigning references to nodes in the code
// If you need to set a reference at runtime, simply leave the global variable empty and assign it in the 'Ready()' function

// Happy hardcoding o7





using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Xml.Serialization;
using Godot;
using static AstroProp_Runtime;


public partial class AstroProp_Runtime : Node3D
{
    public void SY4(ref Godot.Vector3 PosCartesian, ref Godot.Vector3 VelCartesian, Godot.Vector3 InstantaneousAccel, double MET)
    {
        double curt2 = 1.25992104989;
        double w0 = -(curt2 / (2 - curt2));
        double w1 = (1 / (2 - curt2));

        double c1 = w1 / 2;
        double c4 = w1 / 2;
        double c2 = (w0 + w1) / 2;
        double c3 = (w0 + w1) / 2;

        double d1 = w1;
        double d3 = w1;
        double d2 = w0;

        // PosCartesian = PosCartesian * (float)ScaleConversion("ToRealUnits"); do this when calling le function

       // GD.Print(PosCartesian);
        Godot.Vector3 x1 = PosCartesian + VelCartesian * (float)(c1 * Reference.Dynamics.TimeStep);
        Godot.Vector3 a1 = new Godot.Vector3();
        GravityMain_SOI(x1, MET, ref a1);
        foreach (var CelestialRender in KeplerContainers)
        {
            Godot.Vector3 TempAccel = new Godot.Vector3();
            GravityGradient(CelestialRender, x1, MET, ref TempAccel);
            //GD.Print(TempAccel);
            a1 += TempAccel;
            //Debug.Log(CelestialRender.Name.ToString());
        }
       // GD.Print(a1);
        a1 += InstantaneousAccel;

        Godot.Vector3 v1 = VelCartesian + (float)(d1) * a1 * (float)Reference.Dynamics.TimeStep;
        Godot.Vector3 x2 = x1 + (float)c2 * v1 * (float)Reference.Dynamics.TimeStep;
        
        Godot.Vector3 a2 = new Godot.Vector3();
        GravityMain_SOI(x2, MET, ref a2);
        foreach (var CelestialRender in KeplerContainers)
        {
            Godot.Vector3 TempAccel = new Godot.Vector3();
            GravityGradient(CelestialRender, x2, MET, ref TempAccel);
            a2 += TempAccel;
            //Debug.Log(CelestialRender.Name.ToString());
        }
        a2 += InstantaneousAccel;

        Godot.Vector3 v2 = v1 + (float)(d2) * a2 * (float)Reference.Dynamics.TimeStep;
        Godot.Vector3 x3 = x2 + (float)c3 * v2 * (float)Reference.Dynamics.TimeStep;

        Godot.Vector3 a3 = new Godot.Vector3();
        GravityMain_SOI(x3, MET, ref a3);
        foreach (var CelestialRender in KeplerContainers)
        {
            Godot.Vector3 TempAccel = new Godot.Vector3();
            GravityGradient(CelestialRender, x3, MET, ref TempAccel);
            a3 += TempAccel;
            //Debug.Log(CelestialRender.Name.ToString());
        }
        a3 += InstantaneousAccel;

        Godot.Vector3 v3 = v2 + (float)(d3) * a3 * (float)Reference.Dynamics.TimeStep;
        Godot.Vector3 x4 = x3 + (float)c4 * v3 * (float)Reference.Dynamics.TimeStep;

        Godot.Vector3 v4 = v3;
        //GD.Print(v4 );//* (float)ScaleConversion("ToUnityUnits")
        PosCartesian = x4;
        VelCartesian = v4;
    }
    
    public class StateVectors
    {
        public Godot.Vector3 PosCartesian = new Godot.Vector3();
        public Godot.Vector3 VelCartesian = new Godot.Vector3();

        public Godot.Vector3 PrevPosLerp = new Godot.Vector3();

        //noballs variable lmao:
        public Godot.Vector3 InstantaneousAccel = new Godot.Vector3(); // Instantaneous propulsion by engines to be considered by ship
    }
    public class SegmentStepFrame
    {
        public double MET = 0;

        public Node3D ObjectRef;
        
        
        public Godot.Vector3 PosCartesian = new Godot.Vector3();
        public Godot.Vector3 VelCartesian = new Godot.Vector3();

        public Godot.Vector3 PrevPosLerp = new Godot.Vector3();

            //noballs variable lmao:
        public Godot.Vector3 InstantaneousAccel = new Godot.Vector3(); // Instantaneous propulsion by engines to be considered by ship
        

        public double DeltaV;

        public string Event1; // misc register name
        public string Event2;
        public string Event3;

        public double Data1; // register data
        public double Data2;
        public double Data3;

       
        
        public SegmentStepFrame(
            double MET,
            Godot.Vector3 Prev_Pos,
            Godot.Vector3 Prev_Vel,
            Godot.Vector3 Des_Accel,

            bool PredictPrevPosLerp
           )
        {
            this.MET = MET;
            //  this.PosCartesian;
            this.PosCartesian = Prev_Pos;
            this.VelCartesian = Prev_Vel;
            this.InstantaneousAccel = Des_Accel;

            
        }
    }
    
    public class ProjectOry
    {
        public string Name;
        public string Description;

       

        public int StartMET;
        public int InterruptMET; // public ObjectRef = GetNode<Spatial>("%MyUniqueNodeName");

        public Godot.ImmediateMesh TrackStripMesh;

        public Node3D ParentRef;// = GetNode<Node3D>("Global"); // default is just the global scene (main soi relative)
        public MeshInstance3D ObjectRef; //this is the node3d that the primitive line strip is parented to

        public StateVectors StateVectors = new StateVectors(); // state at the start of the trajectory
        //LineMat.sha = Godot.BaseMaterial3D.ShadingModeEnum = 1;


        public List<SegmentStepFrame> Trajectory;
        public ProjectOry( // rendered trajectory
           Node3D ParentRef,
           string Name,
           string Description,
           Godot.Vector3 PosCartesian,
           Godot.Vector3 VelCartesian,
           Godot.Vector3 InstantaneousAccel,
           int StartMET,
           int InterruptMET
            
           )
        {
            this.ParentRef = ParentRef;
            this.Name = Name;
            this.Description = Description;
            this.StateVectors.PosCartesian = PosCartesian;
            this.StateVectors.VelCartesian = VelCartesian;
            this.StateVectors.InstantaneousAccel = InstantaneousAccel;

            this.StartMET = StartMET;
            this.InterruptMET = InterruptMET;

        }

       

        // create another list of the projected segments, then tell them to remove themselves as their MET is surpassed
    }
    public void SetUpProjectOry(
           ref ProjectOry ProjectOry,
           Node3D ParentRef
          


    )
    {
        GD.Print(
            ProjectOry,
            ParentRef
            );
        Godot.OrmMaterial3D LineMat = new Godot.OrmMaterial3D();
        LineMat.ShadingMode = BaseMaterial3D.ShadingModeEnum.PerPixel; //LineMat.ShadingMode = BaseMaterial3D.ShadingModeEnum.Unshaded;
        LineMat.EmissionEnabled = true;
        LineMat.DisableAmbientLight = true;
        LineMat.DisableReceiveShadows = true;
        LineMat.Emission = Color.FromHtml("#FF14AF");
        LineMat.EmissionIntensity = 10;
        
        LineMat.AlbedoColor = Color.FromHtml("#FF14AF");

        ProjectOry.ParentRef = ParentRef;
        ProjectOry.Name = Name;
        

        ProjectOry.ObjectRef = new Godot.MeshInstance3D(); // the track itself
        ProjectOry.TrackStripMesh = new Godot.ImmediateMesh();
        ProjectOry.TrackStripMesh.SurfaceBegin(Godot.Mesh.PrimitiveType.LineStrip, LineMat);

        
        ProjectOry.ObjectRef.Mesh = ProjectOry.TrackStripMesh;
        // this.StateVectors.InstantaneousAccel = InstantaneousAccel;
        

        Godot.Vector3 LS_P = ProjectOry.StateVectors.PosCartesian;
        Godot.Vector3 LS_V = ProjectOry.StateVectors.VelCartesian;
        ProjectOry.TrackStripMesh.SurfaceAddVertex(LS_P * (float)ScaleConversion("ToUnityUnits"));
        GD.Print(ProjectOry.StartMET, ProjectOry.InterruptMET);
        GD.Print(LS_P * (float)ScaleConversion("ToUnityUnits"));
        
        ProjectOry.Trajectory = new List<SegmentStepFrame>(ProjectOry.InterruptMET - ProjectOry.StartMET);
        for (int i = ProjectOry.StartMET; i < ProjectOry.InterruptMET; i++)
        {
            SY4(ref LS_P, ref LS_V, ProjectOry.StateVectors.InstantaneousAccel, i);
            SegmentStepFrame Iter_Frame = new SegmentStepFrame(i, LS_P, LS_V, ProjectOry.StateVectors.InstantaneousAccel, false);
            //GD.Print(LS_V);
            //GD.Print(i);
            ProjectOry.Trajectory.Add(Iter_Frame); // store instantaneous trajectory data here 
            ProjectOry.TrackStripMesh.SurfaceAddVertex(LS_P * (float)ScaleConversion("ToUnityUnits"));
            //ProjectOry.TrackStripMesh.SurfaceAddVertex(LS_P*(float)ScaleConversion("ToUnityUnits"));
            // GD.Print(LS_P * (float)ScaleConversion("ToUnityUnits"));
        };
        //GD.Print(LS_P * (float)ScaleConversion("ToUnityUnits"));
        // GD.Print("Done rendering");
        //GD.Print(ProjectOry.Trajectory);
        //GD.Print(ProjectOry.TrackStripMesh.SurfaceGetArrays();
        ProjectOry.TrackStripMesh.SurfaceEnd();
        ProjectOry.ObjectRef.TopLevel = true; // disable if relative to another body besides the main soi
        ParentRef.AddChild(ProjectOry.ObjectRef);
        ProjectOry.ObjectRef.Position = new Vector3();
        ProjectOry.ObjectRef.CastShadow = 0;
        // this.NBodyRef = Instantiate(NBodyRef, new Godot.Vector3(0, 0, 0), Godot.Quaternion.identity);
    }
    public void RunBallisticTrack(

           NBodyAffected Vessel
    )
    {
        double OrbitalPeriod = CalculateOrbital(
            new Godot.Vector3(),
            Vessel.StateVectors.PosCartesian,
            Vessel.StateVectors.VelCartesian,
            new Godot.Vector3(),
            Reference.SOI.MainReference.GravitationalParameter
            );
        GD.Print((OrbitalPeriod));
        Vessel.Trajectory = new ProjectOry(
            Vessel.ObjectRef,
            "Ballistic Trajectory",
            "",
            Vessel.StateVectors.PosCartesian,
            Vessel.StateVectors.VelCartesian,
            new Godot.Vector3(),
            (int)Reference.Dynamics.MET,
            259200*3
            //(int)(OrbitalPeriod * 1.5)
            );
        SetUpProjectOry(ref Vessel.Trajectory, Vessel.ObjectRef);
    }
    public class NBodyAffected // ballistic and nonballistic object
    {
        public string Name;
        public string Description;

        public Node3D ObjectRef;

        public StateVectors StateVectors = new StateVectors();

        public ProjectOry Trajectory;

        public NBodyAffected(
            Node3D ObjectRef,
            string Name,
            string Description,
            Godot.Vector3 PosCartesian,
            Godot.Vector3 VelCartesian,
            Godot.Vector3 InstantaneousAccel


            )
        {
            this.ObjectRef = ObjectRef;
            this.Name = Name;
            this.Description = Description;
            this.StateVectors.PosCartesian = PosCartesian;
            this.StateVectors.VelCartesian = VelCartesian;
           // this.StateVectors.InstantaneousAccel = InstantaneousAccel;


            // this.NBodyRef = Instantiate(NBodyRef, new Godot.Vector3(0, 0, 0), Godot.Quaternion.identity);

        }
    }
    public class CelestialRender
    {
        public string Name;
        public string Description;

        public Node3D ObjectRef;//   = GetNode<Node>("Global/Earth");



        // default is moon

        public double Mass = 7.35 * Mathf.Pow(10, 22); //kg
        public double GravitationalParameter = 4.904 * Mathf.Pow(10, 12); // m^3 s^-2
        public double SOI = 64300;


        public double SMA = 384400 * 1000; //meters lmao
        public double Inclination = 5.145 * AstroProp_Runtime.Reference.Dynamics.DegToRads; //degrees to the ecliptic, in rads
        public double Eccentricity = .0549; // just eccentricity
        public double ArgPeri = 0 * AstroProp_Runtime.Reference.Dynamics.DegToRads;
        public double LongAscen = 0 * AstroProp_Runtime.Reference.Dynamics.DegToRads;
        public double MeanAnom = 0;




        public CelestialRender(
            string Name,
            string Description,
            Node3D ObjectRef,
            double Mass,
            double GravitationalParameter,
            double SOI,
            double SMA,
            double Inclination,
            double Eccentricity,
            double ArgPeri,
            double LongAscen,
            double MeanAnom
            )
        {
            this.Name = Name;
            this.Description = Description;
            this.ObjectRef = ObjectRef;
            this.Mass = Mass;
            this.GravitationalParameter = GravitationalParameter;
            this.SOI = SOI;
            this.SMA = SMA;
            this.Inclination = -Inclination;
            this.Eccentricity = Eccentricity;
            this.ArgPeri = ArgPeri;
            this.LongAscen = LongAscen;
            this.MeanAnom = MeanAnom;
        }
    }



    public class Reference
    {
        
        public class Dynamics
        {
            public static double StandardGravParam = (6.67 * Mathf.Pow(10, -11)); // Unmodifiable

            public static double TimeStep = 1 / 1; // seconds in between
            public static double MET = 0; //mean elapsed time

            public static double RandomAssConstant = 1;//9, 8.4 for some odd reason
            public static double TimeCompression = (10000); //100000; // default is 1, 2360448 is 1 lunar month per second

            public static double DegToRads = Math.PI / 180;
        };

        public class SOI
        {
            // first class is always the origin soi
            public class MainReference
            {
                public Node3D ObjectRef; //= GetNode<Mesh>("Global/Earth");
                    //GetNode<Mesh>("Global/Earth");

                public static double Mass = 5.97 * Mathf.Pow(10, 24); //kg
                public static double GravitationalParameter = 3.98 * Mathf.Pow(10, 14); // m^3 s^-2
                public static double Radius = 6378.1 * 1000; //Meters
            };





            //public object Moon = new CelestialRender();
        };
    };
    List<CelestialRender> KeplerContainers = new List<CelestialRender>(100); //List<CelestialRender> KeplerContainers = new List<CelestialRender>(1);
    List<NBodyAffected> NByContainers = new List<NBodyAffected>(30);



    //List<CelestialRender> KeplerContainers = new List<CelestialRender>();

    public void ReturnEccentricAnomaly(double M, double e, ref double E)
    {
        //M = E - e*MathF.Sin(E);

        double Tolerance = Mathf.Pow(10, -3);

        // float e = (float) e; //.ToString(); // gives 2.9257


        for (int i = 0; i <= 100; i = i + 2)
        {
            //if (Mathf.Floor(i/20) == i/20)
            double Residual = E - (e * System.Math.Sin(E)) - M;

            double Derivative = 1 - (e * System.Math.Sin(E));

            E = E - Residual / Derivative;
            if (System.Math.Abs(Residual) < Tolerance)
            {
                break;
            }
        }
    }

    public void ModStateVector_Kep(CelestialRender SOI, double MET_Frame, ref Godot.Vector3 PosCartesian, ref Godot.Vector3 VelCartesian)
    {
        double E = 0;

        double dtT = MET_Frame;

        double M = dtT*(Reference.Dynamics.RandomAssConstant * System.Math.Sqrt(SOI.GravitationalParameter / System.Math.Pow(SOI.SMA, 3)));
        // mean anom
        ReturnEccentricAnomaly(M, SOI.Eccentricity, ref E);
        // Class Keplerian = SOI.GetType()

        double TrueAnom = 2 * System.Math.Atan2(System.Math.Sqrt(1 + SOI.Eccentricity) * System.Math.Sin(E / 2), System.Math.Sqrt(1 - SOI.Eccentricity) * System.Math.Cos(E / 2));

        double NewRadius = SOI.SMA * (1 - SOI.Eccentricity * Math.Cos(E));

        // SOI.Inclination -= 90 * Reference.Dynamics.DegToRads;

        // Debug.Log(TrueAnom.ToString());
        // Debug.Log(NewRadius.ToString());

        Godot.Vector3 O = new Godot.Vector3(
            (float)(Math.Cos(TrueAnom) * NewRadius),
            (float)(Math.Sin(TrueAnom) * NewRadius),
            0 // UpVector in perifocal, so no value
            );

        double O_Dot_OP = Math.Sqrt(SOI.GravitationalParameter * SOI.SMA) / NewRadius;

        Godot.Vector3 O_Dot = new Godot.Vector3(
            (float)(-Math.Sin(E) * O_Dot_OP),
            (float)(Math.Sqrt(1 - Math.Pow(((float)(SOI.Eccentricity)), ((float)(2)))) * Math.Cos(E) * O_Dot_OP),
            0 // UpVector in perifocal, so no value
            );
        // Debug.Log(O_Dot.ToString());
        PosCartesian = new Godot.Vector3(
           (float)(O.X * (Math.Cos(SOI.ArgPeri) * Math.Cos(SOI.LongAscen) - Math.Sin(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Sin(SOI.LongAscen)))
           -
           (float)(O.Y * (Math.Sin(SOI.ArgPeri) * Math.Cos(SOI.LongAscen) + Math.Cos(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Sin(SOI.LongAscen)))
           ,
           (float)(O.X * (Math.Cos(SOI.ArgPeri) * Math.Sin(SOI.LongAscen) + Math.Sin(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Cos(SOI.LongAscen)))
           +
           (float)(O.Y * (Math.Cos(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Cos(SOI.LongAscen) - Math.Sin(SOI.ArgPeri) * Math.Sin(SOI.LongAscen)))
            ,
           (float)(O.X * Math.Sin(SOI.ArgPeri) * Math.Sin(SOI.Inclination)) + (float)(O.Y * Math.Cos(SOI.ArgPeri) * Math.Sin(SOI.Inclination))
           );

        VelCartesian = new Godot.Vector3(
           (float)(O_Dot.X * (Math.Cos(SOI.ArgPeri) * Math.Cos(SOI.LongAscen) - Math.Sin(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Sin(SOI.LongAscen)))
           -
           (float)(O_Dot.Y * (Math.Sin(SOI.ArgPeri) * Math.Cos(SOI.LongAscen) + Math.Cos(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Sin(SOI.LongAscen)))
           ,
           (float)(O_Dot.X * (Math.Cos(SOI.ArgPeri) * Math.Sin(SOI.LongAscen) + Math.Sin(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Cos(SOI.LongAscen)))
           +
           (float)(O_Dot.Y * (Math.Cos(SOI.ArgPeri) * Math.Cos(SOI.Inclination) * Math.Cos(SOI.LongAscen) - Math.Sin(SOI.ArgPeri) * Math.Sin(SOI.LongAscen)))
            ,
           (float)(O_Dot.X * Math.Sin(SOI.ArgPeri) * Math.Sin(SOI.Inclination)) + (float)(O_Dot.Y * Math.Cos(SOI.ArgPeri) * Math.Sin(SOI.Inclination))
           );


    }
    public void GravityMain_SOI(Godot.Vector3 PosCartesian, double MET_Frame, ref Godot.Vector3 Acceleration)
    {



        Godot.Vector3 SOI_PosCartesian = new Godot.Vector3(0, 0, 0);







        Godot.Vector3 FromProjToSOI = SOI_PosCartesian - PosCartesian;

        double DistanceExponent = Math.Pow(FromProjToSOI.Length(), 2);

        Godot.Vector3 GravityUnit = FromProjToSOI / FromProjToSOI.Length();

        double GravityAccelerant = (Reference.SOI.MainReference.GravitationalParameter / DistanceExponent);

        Acceleration = GravityUnit * (float)GravityAccelerant;

        // return FromProjToSOI;
    }
    public void GravityGradient(CelestialRender SOI, Godot.Vector3 PosCartesian, double MET_Frame, ref Godot.Vector3 Acceleration)
    {



        Godot.Vector3 SOI_PosCartesian = new Godot.Vector3();
        Godot.Vector3 SOI_VelCartesian = new Godot.Vector3();

        ModStateVector_Kep(SOI, MET_Frame, ref SOI_PosCartesian, ref SOI_VelCartesian);





        Godot.Vector3 FromProjToSOI = (SOI_PosCartesian - PosCartesian);
        

        double DistanceExponent = Math.Pow(FromProjToSOI.Length() ,2);

        Godot.Vector3 GravityUnit = FromProjToSOI / FromProjToSOI.Length();

        double GravityAccelerant = (Reference.SOI.MainReference.GravitationalParameter / DistanceExponent);

        Acceleration = GravityUnit * (float)GravityAccelerant;

        // return FromProjToSOI;
    }
    

    public static double ScaleConversion(string Units)
    {
        // Earth Radius is 6378.1 * 1000 meters
        // earth diameter in the astrogation viewport is 10 meters, so radius is 5
        double Coefficient = 1;
        double UnityDefaulUnit_To_RealMeters = (6378.1 * 1000) / 5;
        if (Units == "ToUnityUnits")
        {
            Coefficient = 1 / UnityDefaulUnit_To_RealMeters;
        }
        else
        {
            Coefficient = UnityDefaulUnit_To_RealMeters;
        }
        return Coefficient;
    }

    public static double CalculateOrbital(Godot.Vector3 SOI_Pos, Godot.Vector3 Subject_Pos, Godot.Vector3 SOI_Vel, Godot.Vector3 Subject_Vel, double SOI_Mu)
    {
        double Relative_Vel = (Subject_Vel - SOI_Vel).Length();
        double Distance = (Subject_Pos - SOI_Pos).Length();

        double SMA = -(SOI_Mu*Distance)/(Distance*Relative_Vel-2*SOI_Mu);
        double Period = 2 * System.Math.PI * System.Math.Sqrt(Math.Pow(SMA, 3)/SOI_Mu);
        GD.Print(Period);
        return Period;
    }
    public void MoveCelestial(CelestialRender SOI, double MET)
    {
        Godot.Vector3 PosCartesian = new Godot.Vector3();
        Godot.Vector3 VelCartesian = new Godot.Vector3();
        double MET_Frame = MET;

        ModStateVector_Kep(SOI, MET_Frame, ref PosCartesian, ref VelCartesian);

        double LocalScale = ScaleConversion("ToUnityUnits");
        // LocalScale = (float)(LocalScale);

        PosCartesian = new Godot.Vector3((float)(PosCartesian.X * LocalScale), (float)(PosCartesian.Y * LocalScale), (float)(PosCartesian.Z * LocalScale));

       // SOI.ObjectRef.Translate(PosCartesian);
        SOI.ObjectRef.Position = (PosCartesian);
        // GD.Print(PosCartesian);
        //.transform.position = PosCartesian;
    }

    public void MoveNBy(NBodyAffected Object, double MET)
    {
        float LocalScale = (float)ScaleConversion("ToUnityUnits");

        Godot.Vector3 PosCartesian = Object.StateVectors.PosCartesian;
        Godot.Vector3 VelCartesian = Object.StateVectors.VelCartesian;
        //double MET_Frame = MET;
        Object.StateVectors.PrevPosLerp = PosCartesian;
        SY4(ref PosCartesian, ref VelCartesian, Object.StateVectors.InstantaneousAccel, MET);

       
        // LocalScale = (float)(LocalScale);

        
        // SOI.ObjectRef.Translate(PosCartesian);
        //Object.ObjectRef.Position = (new Godot.Vector3((float)(PosCartesian.X * LocalScale), (float)(PosCartesian.Y * LocalScale), (float)(PosCartesian.Z * LocalScale)));


        Object.StateVectors.PosCartesian = PosCartesian;
        Object.StateVectors.VelCartesian = VelCartesian;
        // GD.Print(PosCartesian);
        //.transform.position = PosCartesian;
    }

    void BeginStepOps()
    {
        if (Reference.Dynamics.MET == (2360448))
        {
            GD.Print("Lunar orbit debug");
            Reference.Dynamics.TimeCompression = 1;
        }

        Reference.Dynamics.MET += Reference.Dynamics.TimeStep;
        // Begin to do the Celestials 
        // naaa fuck that

        // Perform N-Body Ballistics

        foreach (var NBodyAffected in NByContainers)
        {
            MoveNBy(NBodyAffected, Reference.Dynamics.MET);
            //Debug.Log(CelestialRender.Name.ToString());
        }

    }

    // Start is called before the first frame update
    public override void _Ready()
    {
        Reference.Dynamics.TimeCompression = 1000;
        GD.Print("Boostrapper Startup");
        // Line ends after this, wtf???
        GD.Print(GetNode<Node3D>("Global/Moon").Position);
        
        GetNode<Node3D>("Global/Moon").Position = (new Vector3(0, 0, 0));
        GD.Print(GetNode<Node3D>("Global/Moon").Position);
        //=(new Godot.Vector3(10, 0, 0));
        GD.Print("Can Translate");                                           //(new Godot.Vector3(10, 0, 0));
                                                                        // Reference.SOI.MainReference. = GetNode<Node>("Global/Earth");
                                                                        //Print("AstroPropV2");
                                                                        //MotionContainer.
                                                                        //KeplerContainers.Clear();
        KeplerContainers.Add(new CelestialRender(
            "Moon",
            "Our closest rock",
            GetNode<Node3D>("Global/Moon"),
            7.35 * System.Math.Pow(10, 22),
            4.904 * System.Math.Pow(10, 12),
            64300,
            384400 * 1000,
            (5.145 + 90) * AstroProp_Runtime.Reference.Dynamics.DegToRads,
            .0549, //.0549
            0 * AstroProp_Runtime.Reference.Dynamics.DegToRads,
            0 * AstroProp_Runtime.Reference.Dynamics.DegToRads,
            0

        ));
        NByContainers.Add(new NBodyAffected(
            GetNode<Node3D>("Global/Sagitta"),
            "Sagitta",
            "a spaceship",
            new Godot.Vector3(410789000, 0, 0), //6789000 iss altitude meters
            new Godot.Vector3(0, 00, -1600), //-6576 iss velocity m/s
            new Godot.Vector3(0, 0, 0) // zero propulsion


        ));
        foreach (var NBodyAffected in NByContainers)
        {
            RunBallisticTrack(NBodyAffected);
        }
        // Debug.Log(Reference.SOI.KeplerContainer);
        // Reference.SOI.PrintProperties()

    }
    public static double LastStep = 0;
    public static double NextStep = Time.GetUnixTimeFromSystem() + Reference.Dynamics.TimeStep / (Reference.Dynamics.TimeCompression);
    
    public override void _Process(double Delta)
    {
       // Reference.Dynamics.TimeCompression = GetNode<Control>("Control").GetMeta("TimeCompression", new float = 1.0)); just change from witin another script
        double RealTimeInterpolate = Reference.Dynamics.MET;
        bool debug = false;

        
       // GD.Print((int)Math.Ceiling((Time.GetUnixTimeFromSystem() - NextStep) / (1 / (Reference.Dynamics.TimeCompression))));
        if (Time.GetUnixTimeFromSystem() >= NextStep & (! debug)) 
        {
            int OverfillFrame = (int)Math.Ceiling((Time.GetUnixTimeFromSystem() - NextStep) / (1 / (Reference.Dynamics.TimeCompression))); //(int)Math.Ceiling(Time.time - NextStep);
          // overfillframe just gives UnixTime for some reason?? fixed --10052023 jcr
            LastStep = Time.GetUnixTimeFromSystem();
            NextStep = (Time.GetUnixTimeFromSystem()) + Reference.Dynamics.TimeStep / (Reference.Dynamics.TimeCompression);
            for (int i = 0; i < OverfillFrame; i++)
            {
                
                BeginStepOps();
            
            };

            RealTimeInterpolate = Reference.Dynamics.MET;
        }
        else
        {
            RealTimeInterpolate = Reference.Dynamics.MET + (Time.GetUnixTimeFromSystem() - NextStep);

        }
       // Reference.Dynamics.TimeCompression = 2;
        foreach (var NBodyAffected in NByContainers)
        {
            //GD.Print(RealTimeInterpolate);
            float LerpFloat = (float)(1-(Reference.Dynamics.MET-RealTimeInterpolate));
            Godot.Vector3 LerpV3 = NBodyAffected.StateVectors.PosCartesian - NBodyAffected.StateVectors.PrevPosLerp;
           // GD.Print(LerpFloat,NBodyAffected.StateVectors.PrevPosLerp);
            NBodyAffected.ObjectRef.Position = (float)ScaleConversion("ToUnityUnits")*(NBodyAffected.StateVectors.PrevPosLerp + LerpV3 * LerpFloat);



            //NBodyAffected.ObjectRef.Position = NBodyAffected.StateVectors.PosCartesian*(float)ScaleConversion("ToUnityUnits") + LerpV3 * LerpFloat;
            //GD.Print(NBodyAffected.StateVectors.PosCartesian);
            //Debug.Log(CelestialRender.Name.ToString());
        }

        foreach (var CelestialRender in KeplerContainers)
        {
            MoveCelestial(CelestialRender, RealTimeInterpolate);
            //Debug.Log(CelestialRender.Name.ToString());
        }

        // Debug.Log((RealTimeInterpolate));

        // AstroProp_Runtime.Reference.Dynamics.StandardGravParam += 1;
        // Debug.Log(AstroProp_Runtime.Reference.Dynamics.StandardGravParam.ToString()); ;
    }
}

