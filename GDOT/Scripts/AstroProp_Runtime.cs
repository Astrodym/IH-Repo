// Jonathan C. Ross 
// 09 01 2023 
// All Rights Reserved

// Fuck you, Unity.







using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Xml.Serialization;
using Godot;


public partial class AstroProp_Runtime : Node3D
{
    public class SegmentStepFrame
    {
        public Double MET = 0;

        public class StateVectors
        {
            public Godot.Vector3 PosCartesian = new Godot.Vector3();
            public Godot.Vector3 VelCartesian = new Godot.Vector3();

            public Godot.Vector3 NextPosLerp = new Godot.Vector3();

            //noballs variable lmao:
            public Godot.Vector3 InstantaneousAccel = new Godot.Vector3(); // Instantaneous propulsion by engines to be considered by ship
        }


        public string Event1; // misc register name
        public string Event2;
        public string Event3;

        public double Data1; // register data
        public double Data2;
        public double Data3;
    }
    public class ProjectOry
    {
        public string Name;
        public string Description;

        // public ObjectRef = GetNode<Spatial>("%MyUniqueNodeName");

        public Node3D ObjectRef;

        List<SegmentStepFrame> Trajectory = new List<SegmentStepFrame>();

        // create another list of the projected segments, then tell them to remove themselves as their MET is surpassed
    }

    public class NBodyAffected // ballistic and nonballistic object
    {
        public string Name;
        public string Description;

        public Node3D ObjectRef;

        public SegmentStepFrame.StateVectors StateVectors = new SegmentStepFrame.StateVectors();

        public ProjectOry Trajectory = new ProjectOry();

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
            this.Inclination = Inclination;
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

            public static double RandomAssConstant = 8.4;
            public static double TimeCompression = 1; //100000; // default is 1, 2548800 is 1 lunar month per second

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

        double M = dtT * System.Math.Sqrt(SOI.GravitationalParameter / System.Math.Pow(SOI.SMA, 3));
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
    public void GravityMain_SOI(Godot.Vector3 PosCartesian, double MET_Frame, Godot.Vector3 Acceleration)
    {



        Godot.Vector3 SOI_PosCartesian = new Godot.Vector3(0, 0, 0);







        Godot.Vector3 FromProjToSOI = SOI_PosCartesian - PosCartesian;

        double DistanceExponent = FromProjToSOI.Length() * Math.Exp(2);

        Godot.Vector3 GravityUnit = new Godot.Vector3(
            FromProjToSOI.X / FromProjToSOI.Length(),
            FromProjToSOI.Y / FromProjToSOI.Length(),
            FromProjToSOI.Z / FromProjToSOI.Length()
        );

        double GravityAccelerant = (Reference.SOI.MainReference.GravitationalParameter / DistanceExponent);

        Acceleration = new Godot.Vector3(
            (float)(FromProjToSOI.X * GravityAccelerant),
            (float)(FromProjToSOI.Y * GravityAccelerant),
            (float)(FromProjToSOI.Z * GravityAccelerant)
        );

        // return FromProjToSOI;
    }
    public void GravityGradient(CelestialRender SOI, Godot.Vector3 PosCartesian, double MET_Frame, Godot.Vector3 Acceleration)
    {



        Godot.Vector3 SOI_PosCartesian = new Godot.Vector3();
        Godot.Vector3 SOI_VelCartesian = new Godot.Vector3();

        ModStateVector_Kep(SOI, MET_Frame, ref SOI_PosCartesian, ref SOI_VelCartesian);





        Godot.Vector3 FromProjToSOI = SOI_PosCartesian - PosCartesian;

        double DistanceExponent = FromProjToSOI.Length() * Math.Exp(2);

        Godot.Vector3 GravityUnit = new Godot.Vector3(
            FromProjToSOI.X / FromProjToSOI.Length(),
            FromProjToSOI.Y / FromProjToSOI.Length(),
            FromProjToSOI.Z / FromProjToSOI.Length()
        );

        double GravityAccelerant = (SOI.GravitationalParameter / DistanceExponent);

        Acceleration = new Godot.Vector3(
            (float)(FromProjToSOI.X * GravityAccelerant),
            (float)(FromProjToSOI.Y * GravityAccelerant),
            (float)(FromProjToSOI.Z * GravityAccelerant)
        );

        // return FromProjToSOI;
    }

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


        Godot.Vector3 x1 = PosCartesian + VelCartesian * (float)(c1 * Reference.Dynamics.TimeStep);
        Godot.Vector3 a1 = new Godot.Vector3();
        GravityMain_SOI(x1, MET, a1);
        foreach (var CelestialRender in KeplerContainers)
        {
            Godot.Vector3 TempAccel = new Godot.Vector3();
            GravityGradient(CelestialRender, x1, MET, TempAccel);
            a1 += TempAccel;
            //Debug.Log(CelestialRender.Name.ToString());
        }
        a1 += InstantaneousAccel;

        Godot.Vector3 v1 = VelCartesian + (float)(d1) * a1 * (float)Reference.Dynamics.TimeStep;
        Godot.Vector3 x2 = x1 + (float)c2 * v1 * (float)Reference.Dynamics.TimeStep;

        Godot.Vector3 a2 = new Godot.Vector3();
        GravityMain_SOI(x2, MET, a2);
        foreach (var CelestialRender in KeplerContainers)
        {
            Godot.Vector3 TempAccel = new Godot.Vector3();
            GravityGradient(CelestialRender, x2, MET, TempAccel);
            a2 += TempAccel;
            //Debug.Log(CelestialRender.Name.ToString());
        }
        a2 += InstantaneousAccel;

        Godot.Vector3 v2 = v1 + (float)(d2) * a2 * (float)Reference.Dynamics.TimeStep;
        Godot.Vector3 x3 = x2 + (float)c3 * v2 * (float)Reference.Dynamics.TimeStep;

        Godot.Vector3 a3 = new Godot.Vector3();
        GravityMain_SOI(x3, MET, a3);
        foreach (var CelestialRender in KeplerContainers)
        {
            Godot.Vector3 TempAccel = new Godot.Vector3();
            GravityGradient(CelestialRender, x3, MET, TempAccel);
            a3 += TempAccel;
            //Debug.Log(CelestialRender.Name.ToString());
        }
        a3 += InstantaneousAccel;

        Godot.Vector3 v3 = v2 + (float)(d3) * a3 * (float)Reference.Dynamics.TimeStep;
        Godot.Vector3 x4 = x3 + (float)c4 * v3 * (float)Reference.Dynamics.TimeStep;

        Godot.Vector3 v4 = v3;

        PosCartesian = x4;
        VelCartesian = v4;
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
        return SMA;
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
        Godot.Vector3 PosCartesian = Object.StateVectors.PosCartesian;
        Godot.Vector3 VelCartesian = Object.StateVectors.VelCartesian;
        double MET_Frame = MET;

        SY4(ref PosCartesian, ref VelCartesian, Object.StateVectors.InstantaneousAccel, MET);

        double LocalScale = ScaleConversion("ToUnityUnits");
        // LocalScale = (float)(LocalScale);

        PosCartesian = new Godot.Vector3((float)(PosCartesian.X * LocalScale), (float)(PosCartesian.Y * LocalScale), (float)(PosCartesian.Z * LocalScale));

        // SOI.ObjectRef.Translate(PosCartesian);
        Object.ObjectRef.Position = (PosCartesian);


        Object.StateVectors.PosCartesian = PosCartesian;
        Object.StateVectors.VelCartesian = VelCartesian;
        // GD.Print(PosCartesian);
        //.transform.position = PosCartesian;
    }

    void BeginStepOps()
    {
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
            "Our second closest rock",
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
            "a dumbfuck",
            new Godot.Vector3(1000000,0,0),
            new Godot.Vector3(9000, 0, 0),
            new Godot.Vector3(000, 0, 0)


        ));

        // Debug.Log(Reference.SOI.KeplerContainer);
        // Reference.SOI.PrintProperties()

    }
    public static double LastStep = 0;
    public static double NextStep = Reference.Dynamics.TimeStep;
    // Update is called once per frame
    public override void _Process(double Delta)
    {

        double RealTimeInterpolate = Reference.Dynamics.MET;

        if (Time.GetUnixTimeFromSystem() >= NextStep)
        {
            int OverfillFrame = (int)Math.Ceiling((Time.GetUnixTimeFromSystem() - NextStep) / (1 / (Reference.Dynamics.TimeCompression * Reference.Dynamics.RandomAssConstant))); //(int)Math.Ceiling(Time.time - NextStep);
            // Debug.Log((Time.time - NextStep)/ (1 / Reference.Dynamics.TimeCompression));
            //Debug.Log(Time.time + ">=" + NextStep);
            // Change the next update (current second+1)
            LastStep = Time.GetUnixTimeFromSystem();
            NextStep = (Time.GetUnixTimeFromSystem()) + Reference.Dynamics.TimeStep / (Reference.Dynamics.TimeCompression * Reference.Dynamics.RandomAssConstant);
            // Call your fonction

            // loop begins to empty the overfill of the timecompression (actions occurring more often than the frames can accommodate, so loop to catch up)
            for (int i = 0; i < OverfillFrame; i++)
            {
                // Debug.Log("Calc");
                BeginStepOps();
                //Debug.Log(OverfillFrame);
            };

            RealTimeInterpolate = Reference.Dynamics.MET;
        }
        else
        {
            RealTimeInterpolate = Reference.Dynamics.MET + (Time.GetUnixTimeFromSystem() - NextStep);

        }

        foreach (var NBodyAffected in NByContainers)
        {
            float LerpFloat = (float)RealTimeInterpolate;
            Godot.Vector3 LerpV3 = NBodyAffected.StateVectors.NextPosLerp - NBodyAffected.StateVectors.PosCartesian;

            NBodyAffected.ObjectRef.Position = NBodyAffected.StateVectors.PosCartesian + LerpV3 * LerpFloat;
            
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

