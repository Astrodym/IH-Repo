shader_type spatial;

//This shader takes in a mask and an image texture, then adjusts colors based on the mask and a temperature
//Perhaps a normal map could be used for the intensity of the radiated energy.
//For now temperature can be a value, possibly alongside the thermal energy being emitted by the radiator
//Perhaps temperature could be stored in some kind of map (data texture) which the cpu is doing thermal physics things with

uniform float temperature; //Temperature of the radiator
varying float wavelength; //Wavelength of the light being emitted

//A function that uses Planck's Law to calculate the color of light emitted by a black body at a given temperature (https://www.photo-mark.com/notes/plancks-despair)
float planck(in float t, inout float w) {
	float h = 6.62607004e-34; //Planck's constant
	float c = 299792458.; //Speed of light
	float k = 1.38064852e-23; //Boltzmann constant
	float e = 2.7182818284590452353602874713527; //Euler's number
	float a = 2. * h * c * c;
	float b = h * c / (w * k * t);
	return a / (pow(w, 5.) * (pow(e, b) - 1.));
}
void fragment() {
	wavelength = planck(temperature, wavelength);
	//ALBEDO.r = wavelength;
}
void light() {

}