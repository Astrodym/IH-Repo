shader_type spatial;
render_mode cull_disabled;



uniform vec3 StartPos = vec3(0,10,3);
uniform vec3 EndPos = vec3(1,10,0);



uniform float LineThickness = .0001;

uniform vec3 Color = vec3(0,0,1);
//const StartPos = get_meta("StartPos");

vec3 YieldCross(vec3 CameraPos,vec3 PointPos,vec3 TerminalPos){
	vec3 CameraToVector = CameraPos-PointPos;
	vec3 InlineVector = TerminalPos-PointPos;
	
	float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));
	float IVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));
	CameraToVector = CameraToVector/CTVMag;
	InlineVector = InlineVector/IVMag;
	//CameraToVector = CameraToVector/CTVMag;
	
	//CameraToVector.x = CameraToVector.x/CTVMag;
	//CameraToVector.y = CameraToVector.y/CTVMag;
	//CameraToVector.z = CameraToVector.z/CTVMag;
	
	vec3 a = CameraToVector;
	vec3 b = InlineVector;
	vec3 Cross;
	
	
	Cross.x = (a.y + b.z) - (a.z + b.y);
	Cross.y = (a.z + b.x) - (a.x + b.z);
	Cross.z = (a.x + b.y) - (a.y + b.x)	;
	
	float CrossMag = sqrt(pow(Cross.x,2) + pow(Cross.y,2) + pow(Cross.z,2));
	//sqrt((Cross.x^2 + Cross.y^2 + Cross.z^2));
	Cross = Cross/CrossMag;
	vec3 ReturnCross = Cross;
	
	return ReturnCross;
}

//vec3 SetVertex(Vertex, StartPos, EndPos, Top){

void vertex() {
	// Sets: 0,2 and 1,3
	// 2 is bottom left, 0 is top left
	if (VERTEX_ID == 0){
		vec3 Cross = YieldCross(CAMERA_POSITION_WORLD,StartPos,EndPos);
		vec3 CameraToVector = StartPos-CAMERA_POSITION_WORLD;
		float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));

		float Coefficient = (LineThickness*CTVMag);
		vec3 NewPos = StartPos+(Cross*Coefficient);
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
	}
	//here comes the debug
	if (VERTEX_ID == 2){
		vec3 Cross = YieldCross(CAMERA_POSITION_WORLD,StartPos,EndPos);
		vec3 CameraToVector = StartPos-CAMERA_POSITION_WORLD;
		float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));

		float Coefficient = (LineThickness*CTVMag);
		vec3 NewPos = StartPos-(Cross*Coefficient);
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
	}
	if (VERTEX_ID == 1){
		vec3 Cross = YieldCross(CAMERA_POSITION_WORLD,EndPos,StartPos);
		vec3 CameraToVector = EndPos-CAMERA_POSITION_WORLD;
		float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));

		float Coefficient = (LineThickness*CTVMag);
		vec3 NewPos = EndPos+(Cross*Coefficient);
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
	}
	//here comes the debug
	if (VERTEX_ID == 3){
		vec3 Cross = YieldCross(CAMERA_POSITION_WORLD,EndPos,StartPos);
		vec3 CameraToVector = EndPos-CAMERA_POSITION_WORLD;
		float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));

		float Coefficient = (LineThickness*CTVMag);
		vec3 NewPos = EndPos-(Cross*Coefficient);
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
	}
}

void fragment() {
	// Place fragment code here.
	ALBEDO = Color;
	EMISSION = Color;
	RIM = 1.0;
	RIM_TINT = 1.0;
	
	
	// set_vertex( int indx, Vector3 vertex )
	
	//vec3 StartPos =  
}
