shader_type spatial;

uniform vec3 StartPos = vec3(30,0,5);
uniform vec3 EndPos = vec3(-30,0,5);

uniform float LineThickness = 1;
//const StartPos = get_meta("StartPos");

vec3 YieldCross(vec3 CameraPos,vec3 PointPos,vec3 TerminalPos){
	vec3 CameraToVector = PointPos-CameraPos;
	vec3 InlineVector = TerminalPos-PointPos;
	
	float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));
	
	//CameraToVector = CameraToVector/CTVMag;
	
	//CameraToVector.x = CameraToVector.x/CTVMag;
	//CameraToVector.y = CameraToVector.y/CTVMag;
	//CameraToVector.z = CameraToVector.z/CTVMag;
	
	vec3 a = CameraToVector;
	vec3 b = InlineVector;
	vec3 Cross;
	
	
	Cross.x = (a.y * b.z) - (a.z * b.y);
	Cross.y = (a.z * b.x) - (a.x * b.z);
	Cross.z = (a.x * b.y) - (a.y * b.x);
	
	float CrossMag = sqrt(pow(Cross.x,2) + pow(Cross.y,2) + pow(Cross.z,2));
	//sqrt((Cross.x^2 + Cross.y^2 + Cross.z^2));
	Cross = Cross/CrossMag;
	vec3 ReturnCross;
	
	ReturnCross.x = Cross.x;
	ReturnCross.y = Cross.y;
	ReturnCross.z = Cross.z;
	//Cross.x = Cross.x/CrossMag;
	//Cross.y = Cross.y/CrossMag;
	//Cross.z = Cross.z/CrossMag;
	
	// what the FUCK godot?
	
	//vec3 Cross = CameraToVector.cross(InlineVector);
	//CameraToVector.cross(InlineVector);
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
		//sqrt((Cross.x^2 + Cross.y^2 + Cross.z^2));
		//Cross = Cross/CrossMag;
		CTVMag = 1;
		vec3 NewPos = StartPos+(Cross*(LineThickness*CTVMag));
		
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
		//VERTEX.x += StartPos.x;
	}
	if (VERTEX_ID == 1){
		vec3 Cross = YieldCross(CAMERA_POSITION_WORLD,EndPos,StartPos);
		vec3 CameraToVector = EndPos-CAMERA_POSITION_WORLD;
		float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));
		//sqrt((Cross.x^2 + Cross.y^2 + Cross.z^2));
		//Cross = Cross/CrossMag;
		CTVMag = 1;
		vec3 NewPos = EndPos-(Cross*(LineThickness*CTVMag));
		
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
		//VERTEX.x += StartPos.x;
	}
	if (VERTEX_ID == 2){
		vec3 Cross = YieldCross(CAMERA_POSITION_WORLD,StartPos,EndPos);
		vec3 CameraToVector = StartPos-CAMERA_POSITION_WORLD;
		float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));
		//sqrt((Cross.x^2 + Cross.y^2 + Cross.z^2));
		//Cross = Cross/CrossMag;
		CTVMag = 1;
		vec3 NewPos = StartPos-(Cross*(LineThickness*CTVMag));
		
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
		//VERTEX.x += StartPos.x;
	}
	if (VERTEX_ID == 3){
		vec3 Cross = YieldCross(CAMERA_POSITION_WORLD,EndPos,StartPos);
		vec3 CameraToVector = EndPos-CAMERA_POSITION_WORLD;
		float CTVMag = sqrt(pow(CameraToVector.x,2) + pow(CameraToVector.y,2) + pow(CameraToVector.z,2));
		//sqrt((Cross.x^2 + Cross.y^2 + Cross.z^2));
		//Cross = Cross/CrossMag;
		CTVMag = 1;
		vec3 NewPos = EndPos+(Cross*(LineThickness*CTVMag));
		
		VERTEX.x = NewPos.x;
		VERTEX.y = NewPos.y;
		VERTEX.z = NewPos.z;
		//VERTEX.x += StartPos.x;
	}
	



}

void fragment() {
	// Place fragment code here.
	
	// set_vertex( int indx, Vector3 vertex )
	
	//vec3 StartPos =  
}
