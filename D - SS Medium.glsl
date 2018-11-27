/*===============================================================================*\
|########################        [Ishiiruka FX 0.8]        ######################||
|| Credist to:                                                                   ||
|| Asmodean (DolphinFX)                                                          ||
|| Matso (MATSODOF)                                                              ||
|| Gilcher Pascal aka Marty McFly (MATSODOF original port to MCFX)               ||
|| Daniel Rákos (Efficient Gaussian blur with linear sampling)                   ||
|| mudlord (FXAA)                                                                ||
|| Cozimo (Default shader settings for Twilight Princess Mod presets)            ||
|| Drakonas (Integrated DolphinFX effects for Twilight Princess Mod Team)        ||
|| Drakonas (With Tino's help, also integrated ReShade's Gaussian AnamFlare)     ||
|############################        By Tino          ############################|
\*===============================================================================*/
/*
[configuration]
[OptionBool]
GUIName = MATSO DOF
OptionName = B_MATSODOF
DefaultValue = True

[OptionBool]
GUIName = Use depth range focus
OptionName = DOF_A_FOCUSPOINT_RANGE
DefaultValue = True
DependentOption = B_MATSODOF

[OptionRangeFloat]
GUIName = Focus Point
OptionName = DOF_B_FOCUSPOINT
MinValue = 0.0, 0.0
MaxValue = 1.0, 1.0
DefaultValue = 0.18, 0.72
StepAmount = 0.01, 0.01
DependentOption = B_MATSODOF

[OptionRangeFloat]
GUIName = Near Blue Curve
OptionName = DOF_NEARBLURCURVE
MinValue = 0.4
MaxValue = 2.0
DefaultValue = 2
StepAmount = 0.01
DependentOption = B_MATSODOF

[OptionRangeFloat]
GUIName = FarBlue Curve
OptionName = DOF_FARBLURCURVE
MinValue = 0.4
MaxValue = 2.0
DefaultValue = 0.87
StepAmount = 0.01
DependentOption = B_MATSODOF

[OptionRangeFloat]
GUIName = Blur Radius
OptionName = DOF_BLURRADIUS
MinValue = 5.0
MaxValue = 50.0
DefaultValue = 12
StepAmount = 1.0
DependentOption = B_MATSODOF

[OptionRangeInteger]
GUIName = Vignette
OptionName = DOF_VIGNETTE
MinValue = 0
MaxValue = 1000
DefaultValue = 0
StepAmount = 10
DependentOption = B_MATSODOF

[OptionBool]
GUIName = CROMA Enable
OptionName = bMatsoDOFChromaEnable
DefaultValue = True
DependentOption = B_MATSODOF

[OptionRangeFloat]
GUIName = Chroma Pow
OptionName = fMatsoDOFChromaPow
MinValue = 0.2
MaxValue = 3.0
DefaultValue = 0.65
StepAmount = 0.01
DependentOption = B_MATSODOF

[OptionRangeFloat]
GUIName = Bokeh Curve
OptionName = fMatsoDOFBokehCurve
MinValue = 0.5
MaxValue = 20.0
DefaultValue = 8.16
StepAmount = 0.01
DependentOption = B_MATSODOF

[OptionRangeFloat]
GUIName = Bokeh Light
OptionName = fMatsoDOFBokehLight
MinValue = 0.0
MaxValue = 2.0
DefaultValue = 0.037
StepAmount = 0.001
DependentOption = B_MATSODOF

[OptionRangeInteger]
GUIName = Bokeh Quality
OptionName = iMatsoDOFBokehQuality
MinValue = 1
MaxValue = 10
DefaultValue = 2
StepAmount = 1
DependentOption = B_MATSODOF
ResolveAtCompilation = True

[OptionRangeInteger]
GUIName = Bokeh Angle
OptionName = fMatsoDOFBokehAngle
MinValue = 0
MaxValue = 360
DefaultValue = 42
StepAmount = 1
DependentOption = B_MATSODOF

[OptionBool]
GUIName = Bloom
OptionName = C_BLOOM
DefaultValue = true

[OptionBool]
GUIName = Bloom Only
OptionName = A_BLOOMONLY
DefaultValue = false
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Bloom Width
OptionName = A_BLOOMWIDTH
MinValue = 0.5
MaxValue = 3.0
StepAmount = 0.01
DefaultValue = 1.05
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Bloom Power
OptionName = B_BLOOMPOWER
MinValue = 1.0
MaxValue = 8.0
StepAmount = 0.1
DefaultValue = 8
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Bloom Intensity
OptionName = C_BLOOMINTENSITY
MinValue = 0.1
MaxValue = 1.0
StepAmount = 0.01
DefaultValue = 0.12
DependentOption = C_BLOOM

[OptionBool]
GUIName = Light Scattering
OptionName = D_SCATTERRING
DefaultValue = True
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Density
OptionName = E_SDENSITY
MinValue = 0.0
MaxValue = 1.0
StepAmount = 0.01
DefaultValue = 0.12
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Start
OptionName = F_SSTART
MinValue = 0.0
MaxValue = 0.5
StepAmount = 0.01
DefaultValue = 0.09
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = End
OptionName = G_SEND
MinValue = 0.5
MaxValue = 2.0
StepAmount = 0.01
DefaultValue = 2.0
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Default Color
OptionName = H_SCOLOR
MinValue = 0.1, 0.1, 0.1
MaxValue = 1.0, 1.0, 1.0
StepAmount = 0.01, 0.01, 0.01
DefaultValue = 0.53, 0.89, 0.92
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Scattering Intensity
OptionName = I_SINTENSITY
MinValue = 0.0
MaxValue = 1.0
StepAmount = 0.01
DefaultValue = 0.58
DependentOption = C_BLOOM

[OptionBool]
GUIName = Scene Tonemapping
OptionName = C_TONEMAP_PASS
DefaultValue = True

[OptionRangeInteger]
GUIName = TonemapType
OptionName = A_TONEMAP_TYPE
MinValue = 0
MaxValue = 3
StepAmount = 1
DefaultValue = 0
DependentOption = C_TONEMAP_PASS

[OptionRangeInteger]
GUIName = FilmOperator
OptionName = A_TONEMAP_FILM
MinValue = 0
MaxValue = 1
StepAmount = 1
DefaultValue = 1
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = ToneAmount
OptionName = B_TONE_AMOUNT
MinValue = 0.05
MaxValue = 2.00
StepAmount = 0.01
DefaultValue = 0.37
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = FilmStrength
OptionName = B_TONE_FAMOUNT
MinValue = 0.00
MaxValue = 1.00
StepAmount = 0.01
DefaultValue = 0.36
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = BlackLevels
OptionName = C_BLACK_LEVELS
MinValue = 0.00
MaxValue = 1.00
StepAmount = 0.01
DefaultValue = 0
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = Exposure
OptionName = D_EXPOSURE
MinValue = 0.50
MaxValue = 1.50
StepAmount = 0.01
DefaultValue = 1.04
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = Luminance
OptionName = E_LUMINANCE
MinValue = 0.50
MaxValue = 1.50
StepAmount = 0.01
DefaultValue = 0.89
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = WhitePoint
OptionName = F_WHITEPOINT
MinValue = 0.50
MaxValue = 1.50
StepAmount = 0.01
DefaultValue = 0.86
DependentOption = C_TONEMAP_PASS

[OptionBool]
GUIName = Colour Correction
OptionName = D_COLOR_CORRECTION
DefaultValue = true

[OptionRangeInteger]
GUIName = CorrectionPalette
OptionName = A_PALETTE
MinValue = 0
MaxValue = 4
StepAmount = 1
DefaultValue = 3
DependentOption = D_COLOR_CORRECTION

[OptionRangeFloat]
GUIName = Channels R|Y|X|H|Y
OptionName = B_RED_CORRECTION
MinValue = 0.10
MaxValue = 8.00
StepAmount = 0.10
DefaultValue = 0.60
DependentOption = D_COLOR_CORRECTION

[OptionRangeFloat]
GUIName = Channels G|X|Y|S|U
OptionName = C_GREEN_CORRECTION
MinValue = 0.10
MaxValue = 8.00
StepAmount = 0.10
DefaultValue = 3.40
DependentOption = D_COLOR_CORRECTION

[OptionRangeFloat]
GUIName = Channels B|Y|Z|V|V
OptionName = D_BLUE_CORRECTION
MinValue = 0.10
MaxValue = 8.00
StepAmount = 0.10
DefaultValue = 0.10
DependentOption = D_COLOR_CORRECTION

[OptionRangeFloat]
GUIName = PaletteStrength
OptionName = E_CORRECT_STR
MinValue = 0.00
MaxValue = 2.00
StepAmount = 0.01
DefaultValue = 1.18
DependentOption = D_COLOR_CORRECTION

[OptionBool]
GUIName = Paint Shading 1
OptionName = E_PAINT_SHADING
DefaultValue = True

[OptionRangeInteger]
GUIName = Paint method
GUIDescription = The algorithm used for paint effect. 1: water painting, 0: oil painting. You may want to readjust the radius between the two.
OptionName = PaintMethod
MinValue = 0
MaxValue = 1
StepAmount = 1
DefaultValue = 0
DependentOption = E_PAINT_SHADING

[OptionRangeInteger]
GUIName = Paint radius
GUIDescription = Radius of the painted effect, increasing the radius also requires longer loops, so higher values require more performance.
OptionName = PaintRadius
MinValue = 2
MaxValue = 8
StepAmount = 1
DefaultValue = 3
ResolveAtCompilation = True
DependentOption = E_PAINT_SHADING

[OptionRangeFloat]
GUIName = Paint radius
GUIDescription = The overall interpolated strength of the paint effect. Where 1.0 equates to 100% strength.
OptionName = PaintStrength
MinValue = 0.0
MaxValue = 1.0
StepAmount = 0.01
DefaultValue = 0.38
DependentOption = E_PAINT_SHADING

[OptionBool]
GUIName = Paint Shading 2
OptionName = E_PAINT_SHADING_2
DefaultValue = True

[OptionRangeInteger]
GUIName = Paint method
GUIDescription = The algorithm used for paint effect. 1: water painting, 0: oil painting. You may want to readjust the radius between the two.
OptionName = PaintMethod2
MinValue = 0
MaxValue = 1
StepAmount = 1
DefaultValue = 0
DependentOption = E_PAINT_SHADING_2

[OptionRangeInteger]
GUIName = Paint radius
GUIDescription = Radius of the painted effect, increasing the radius also requires longer loops, so higher values require more performance.
OptionName = PaintRadius2
MinValue = 2
MaxValue = 8
StepAmount = 1
DefaultValue = 4
ResolveAtCompilation = True
DependentOption = E_PAINT_SHADING_2

[OptionRangeFloat]
GUIName = Paint radius
GUIDescription = The overall interpolated strength of the paint effect. Where 1.0 equates to 100% strength.
OptionName = PaintStrength2
MinValue = 0.0
MaxValue = 1.0
StepAmount = 0.01
DefaultValue = 0.35
DependentOption = E_PAINT_SHADING_2

[OptionBool]
GUIName = Pixel Vibrance
OptionName = F_PIXEL_VIBRANCE
DefaultValue = true

[OptionRangeFloat]
GUIName = Vibrance
OptionName = A_VIBRANCE
MinValue = -0.50
MaxValue = 1.00
StepAmount = 0.01
DefaultValue = 0.38
DependentOption = F_PIXEL_VIBRANCE

[OptionRangeFloat]
GUIName = RedVibrance
OptionName = B_R_VIBRANCE
MinValue = -1.00
MaxValue = 4.00
StepAmount = 0.01
DefaultValue = 4.0
DependentOption = F_PIXEL_VIBRANCE

[OptionRangeFloat]
GUIName = GreenVibrance
OptionName = C_G_VIBRANCE
MinValue = -1.00
MaxValue = 4.00
StepAmount = 0.01
DefaultValue = 0.39
DependentOption = F_PIXEL_VIBRANCE

[OptionRangeFloat]
GUIName = BlueVibrance
OptionName = D_B_VIBRANCE
MinValue = -1.00
MaxValue = 4.00
StepAmount = 0.01
DefaultValue = 2.87
DependentOption = F_PIXEL_VIBRANCE

[OptionBool]
GUIName = FXAA
OptionName = J_FXAA_PASS
DefaultValue = True

[OptionRangeFloat]
GUIName = SubpixelMax
OptionName = A_FXAA_SUBPIX_MAX
MinValue = 0.00
MaxValue = 1.00
StepAmount = 0.01
DefaultValue = 0.68
DependentOption = J_FXAA_PASS

[OptionRangeFloat]
GUIName = EdgeThreshold
OptionName = B_FXAA_EDGE_THRESHOLD
MinValue = 0.010
MaxValue = 0.500
StepAmount = 0.001
DefaultValue = 0.18
DependentOption = J_FXAA_PASS

[OptionRangeInteger]
GUIName = ShowEdgeDetection
OptionName = C_FXAA_SHOW_EDGES
MinValue = 0
MaxValue = 1
StepAmount = 1
DefaultValue = 0
DependentOption = J_FXAA_PASS

[Pass]
EntryPoint = PS_DOF_MatsoDOF1
DependantOption = B_MATSODOF
Input0=ColorBuffer
Input0Filter=Linear
Input0Mode=Clamp
Input1=DepthBuffer
Input1Filter=Nearest
Input1Mode=Clamp
[Pass]
EntryPoint = PS_DOF_MatsoDOF2
DependantOption = B_MATSODOF
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
Input1=DepthBuffer
Input1Filter=Nearest
Input1Mode=Clamp
[Pass]
EntryPoint = PS_DOF_MatsoDOF3
DependantOption = B_MATSODOF
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
Input1=DepthBuffer
Input1Filter=Nearest
Input1Mode=Clamp
[Pass]
EntryPoint = PS_DOF_MatsoDOF4
DependantOption = B_MATSODOF
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
Input1=DepthBuffer
Input1Filter=Nearest
Input1Mode=Clamp
[Pass]
EntryPoint = PostApplyCC
DependantOption = D_COLOR_CORRECTION
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = A_ReduceSize
DependantOption = C_BLOOM
OutputScale = 0.5
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = BloomH
DependantOption = C_BLOOM
OutputScale = 0.25
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = BloomV
DependantOption = C_BLOOM
OutputScale = 0.25
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = BloomH
DependantOption = C_BLOOM
OutputScale = 0.125
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = BloomV
DependantOption = C_BLOOM
OutputScale = 0.125
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = BloomScatering
DependantOption = C_BLOOM
OutputScale = 0.03125
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = BloomMerger
DependantOption = C_BLOOM
Input0=Pass4
Input0Filter=Linear
Input0Mode=Clamp
Input1=Pass9
Input1Filter=Nearest
Input1Mode=Clamp
Input2=Pass10
Input2Filter=Linear
Input2Mode=Clamp
Input3=DepthBuffer
Input3Filter=Nearest
Input3Mode=Clamp
[Pass]
EntryPoint = SceneTonemapping
DependantOption = C_TONEMAP_PASS
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = PaintShader
DependantOption = E_PAINT_SHADING
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = PaintShader2
DependantOption = E_PAINT_SHADING_2
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = PixelVibrance
DependantOption = F_PIXEL_VIBRANCE
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = Fxaa
DependantOption = J_FXAA_PASS
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[/configuration]
*/

#define PI 		3.1415972
#define PIOVER180 	0.017453292


float4 GetMatsoDOFCA(float2 tex, float CoC)
{
	float3 ChromaPOW = float3(1, 1, 1) * GetOption(fMatsoDOFChromaPow) * CoC;
	float3 chroma = pow(float3(0.5, 1.0, 1.5), ChromaPOW) * 0.5;
	tex = (2.0 * tex - 1.0);
	float2 tr = (tex * chroma.r) + 0.5;
	float2 tg = (tex * chroma.g) + 0.5;
	float2 tb = (tex * chroma.b) + 0.5;

	float3 color = float3(SamplePrevLocation(tr).r, SamplePrevLocation(tg).g, SamplePrevLocation(tb).b) * (1.0 - CoC);

	return float4(color, 1.0);
}

float4 GetMatsoDOFBlur(int axis, float2 coord)
{
	float4 tcol = SamplePrevLocation(coord);
	float scenedepth = SampleDepth();
	float scenefocus = 0.0f;
	if (OptionEnabled(DOF_A_FOCUSPOINT_RANGE))
	{
		scenefocus = GetOption(DOF_B_FOCUSPOINT).x;
	}
	else
	{
		scenefocus = SampleDepthLocation(GetOption(DOF_B_FOCUSPOINT));
	}
	float depthdiff = abs(scenedepth - scenefocus);
	if (OptionEnabled(DOF_A_FOCUSPOINT_RANGE))
	{
		if (abs(depthdiff) < GetOption(DOF_B_FOCUSPOINT).y)
		{
			depthdiff *= 0.5f * depthdiff * depthdiff;
		}
	}
	depthdiff = (scenedepth < scenefocus) ? pow(depthdiff, GetOption(DOF_NEARBLURCURVE))*(1.0f + pow(abs(0.5f - coord.x)*depthdiff + 0.1f, 2.0)*GetOption(DOF_VIGNETTE)) : depthdiff;
	depthdiff = (scenedepth > scenefocus) ? pow(depthdiff, GetOption(DOF_FARBLURCURVE)) : depthdiff;

	float2 discRadius = depthdiff * GetOption(DOF_BLURRADIUS)*GetInvResolution()*0.5 / float(iMatsoDOFBokehQuality);

	int passnumber = 1;

	float sf = 0;

	float2 tdirs[4] = { float2(-0.306, 0.739), float2(0.306, 0.739), float2(-0.739, 0.306), float2(-0.739, -0.306) };
	float wValue = (1.0 + pow(length(tcol.rgb) + 0.1, GetOption(fMatsoDOFBokehCurve))) * (1.0 - GetOption(fMatsoDOFBokehLight));	// special recipe from papa Matso ;)

	for (int i = -iMatsoDOFBokehQuality; i < iMatsoDOFBokehQuality; i++)
	{
		float2 taxis = tdirs[axis];

		taxis.x = cos(GetOption(fMatsoDOFBokehAngle)*PIOVER180)*taxis.x - sin(GetOption(fMatsoDOFBokehAngle)*PIOVER180)*taxis.y;
		taxis.y = sin(GetOption(fMatsoDOFBokehAngle)*PIOVER180)*taxis.x + cos(GetOption(fMatsoDOFBokehAngle)*PIOVER180)*taxis.y;

		float2 tdir = float(i) * taxis * discRadius;
		float2 tcoord = coord.xy + tdir.xy;
		float4 ct;
		if (OptionEnabled(bMatsoDOFChromaEnable))
		{
			ct = GetMatsoDOFCA(tcoord.xy, discRadius.x);
		}
		else
		{
			ct = SamplePrevLocation(tcoord.xy);
		}
		float w = 0.0;
		float b = dot(ct.rgb, float3(0.333, 0.333, 0.333)) + length(ct.rgb) + 0.1;
		w = pow(b, GetOption(fMatsoDOFBokehCurve)) + abs(float(i));
		tcol += ct * w;
		wValue += w;
	}

	tcol /= wValue;

	return float4(tcol.xyz, 1.0);
}

void PS_DOF_MatsoDOF1()
{
	SetOutput(GetMatsoDOFBlur(2, GetCoordinates()));
}

void PS_DOF_MatsoDOF2()
{
	SetOutput(GetMatsoDOFBlur(3, GetCoordinates()));
}

void PS_DOF_MatsoDOF3()
{
	SetOutput(GetMatsoDOFBlur(0, GetCoordinates()));
}

void PS_DOF_MatsoDOF4()
{
	SetOutput(GetMatsoDOFBlur(1, GetCoordinates()));
}


/*------------------------------------------------------------------------------
[GLOBALS|FUNCTIONS]
------------------------------------------------------------------------------*/

#define Epsilon (1e-10)
#define lumCoeff float3(0.299f, 0.587f, 0.114f)

//Average relative luminance
float AvgLuminance(float3 color)
{
	return sqrt(dot(color * color, lumCoeff));
}

//Conversion matrices
float3 RGBtoXYZ(float3 rgb)
{
	const float3x3 m = float3x3(
		0.4124564, 0.3575761, 0.1804375,
		0.2126729, 0.7151522, 0.0721750,
		0.0193339, 0.1191920, 0.9503041);

	return mul(rgb, m);
}

float3 XYZtoRGB(float3 xyz)
{
	const float3x3 m = float3x3(
		3.2404542, -1.5371385, -0.4985314,
		-0.9692660, 1.8760108, 0.0415560,
		0.0556434, -0.2040259, 1.0572252);

	return mul(xyz, m);
}

float3 RGBtoYUV(float3 RGB)
{
	const float3x3 m = float3x3(
		0.2126, 0.7152, 0.0722,
		-0.09991, -0.33609, 0.436,
		0.615, -0.55861, -0.05639);

	return mul(RGB, m);
}

float3 YUVtoRGB(float3 YUV)
{
	const float3x3 m = float3x3(
		1.000, 0.000, 1.28033,
		1.000, -0.21482, -0.38059,
		1.000, 2.12798, 0.000);

	return mul(YUV, m);
}

//Converting XYZ to Yxy
float3 XYZtoYxy(float3 xyz)
{
	float3 Yxy;
	float w = 1.0 / (xyz.r + xyz.g + xyz.b);

	Yxy.r = xyz.g;
	Yxy.g = xyz.r * w;
	Yxy.b = xyz.g * w;

	return Yxy;
}

//Converting Yxy to XYZ
float3 YxytoXYZ(float3 Yxy)
{
	float3 xyz;
	float w = 1.0 / Yxy.b;
	xyz.g = Yxy.r;
	xyz.r = Yxy.r * Yxy.g * w;
	xyz.b = Yxy.r * (1.0 - Yxy.g - Yxy.b) * w;
	return xyz;
}

/*------------------------------------------------------------------------------
[SCENE TONE MAPPING CODE SECTION]
------------------------------------------------------------------------------*/

float3 EncodeGamma(float3 color, float gamma)
{
	color = saturate(color);
	color.r = (color.r <= 0.0404482362771082) ?
		color.r / 12.92 : pow((color.r + 0.055) / 1.055, gamma);
	color.g = (color.g <= 0.0404482362771082) ?
		color.g / 12.92 : pow((color.g + 0.055) / 1.055, gamma);
	color.b = (color.b <= 0.0404482362771082) ?
		color.b / 12.92 : pow((color.b + 0.055) / 1.055, gamma);

	return color;
}

float3 FilmicCurve(float3 color)
{
	float3 T = color;
	float tnamn = GetOption(B_TONE_AMOUNT);

	float A = 0.100;
	float B = 0.300;
	float C = 0.100;
	float D = tnamn;
	float E = 0.020;
	float F = 0.300;
	float W = 1.012;

	T.r = ((T.r*(A*T.r + C*B) + D*E) / (T.r*(A*T.r + B) + D*F)) - E / F;
	T.g = ((T.g*(A*T.g + C*B) + D*E) / (T.g*(A*T.g + B) + D*F)) - E / F;
	T.b = ((T.b*(A*T.b + C*B) + D*E) / (T.b*(A*T.b + B) + D*F)) - E / F;

	float denom = ((W*(A*W + C*B) + D*E) / (W*(A*W + B) + D*F)) - E / F;
	float3 white = float3(denom, denom, denom);

	T = T / white;
	color = saturate(T);

	return color;
}

float3 FilmicTonemap(float3 color)
{
	float3 tone = color;

	float3 black = float3(0.0, 0.0, 0.0);
	tone = max(black, tone);

	tone.r = (tone.r * (6.2 * tone.r + 0.5)) / (tone.r * (6.2 * tone.r + 1.66) + 0.066);
	tone.g = (tone.g * (6.2 * tone.g + 0.5)) / (tone.g * (6.2 * tone.g + 1.66) + 0.066);
	tone.b = (tone.b * (6.2 * tone.b + 0.5)) / (tone.b * (6.2 * tone.b + 1.66) + 0.066);

	const float gamma = 2.42;
	tone = EncodeGamma(tone, gamma);

	color = lerp(color, tone, GetOption(B_TONE_FAMOUNT));

	return color;
}

float4 TonemapPass(float4 color)
{
	float luminanceAverage = GetOption(E_LUMINANCE);
	float bmax = max(color.r, max(color.g, color.b));

	float blevel = pow(saturate(bmax), GetOption(C_BLACK_LEVELS));
	color.rgb = color.rgb * blevel;

	if (GetOption(A_TONEMAP_FILM) == 1) { color.rgb = FilmicTonemap(color.rgb); }
	if (GetOption(A_TONEMAP_TYPE) == 1) { color.rgb = FilmicCurve(color.rgb); }

	float3 XYZ = RGBtoXYZ(color.rgb);

	// XYZ -> Yxy conversion
	float3 Yxy = XYZtoYxy(XYZ);

	// (Wt) Tone mapped scaling of the initial wp before input modifiers
	float Wt = saturate(Yxy.r / AvgLuminance(XYZ));

	if (GetOption(A_TONEMAP_TYPE) == 2) { Yxy.r = FilmicCurve(Yxy).r; }

	// (Lp) Map average luminance to the middlegrey zone by scaling pixel luminance
	float Lp = Yxy.r * GetOption(D_EXPOSURE) / (luminanceAverage + Epsilon);

	// (Wp) White point calculated, based on the toned white, and input modifier
	float Wp = abs(Wt) * GetOption(F_WHITEPOINT);

	// (Ld) Scale all luminance within a displayable range of 0 to 1
	Yxy.r = (Lp * (1.0 + Lp / (Wp * Wp))) / (1.0 + Lp);

	// Yxy -> XYZ conversion
	XYZ = YxytoXYZ(Yxy);

	color.rgb = XYZtoRGB(XYZ);
	color.a = AvgLuminance(color.rgb);

	return color;
}


/*------------------------------------------------------------------------------
[PIXEL VIBRANCE CODE SECTION]
------------------------------------------------------------------------------*/

float4 VibrancePass(float4 color)
{
	float vib = GetOption(A_VIBRANCE);
	float luma = AvgLuminance(color.rgb);

	float colorMax = max(color.r, max(color.g, color.b));
	float colorMin = min(color.r, min(color.g, color.b));

	float colorSaturation = colorMax - colorMin;
	float3 colorCoeff = float3(GetOption(B_R_VIBRANCE)*
		vib, GetOption(C_G_VIBRANCE) * vib, GetOption(D_B_VIBRANCE) * vib);

	color.rgb = lerp(float3(luma, luma, luma), color.rgb, (1.0 + (colorCoeff * (1.0 - (sign(colorCoeff) * colorSaturation)))));
	color.a = AvgLuminance(color.rgb);

	return saturate(color); //Debug: return colorSaturation.xxxx;
}

/*------------------------------------------------------------------------------
[FXAA CODE SECTION]
------------------------------------------------------------------------------*/

#define FXAA_QUALITY__PS 9
#define FXAA_QUALITY__P0 1.0
#define FXAA_QUALITY__P1 1.5
#define FXAA_QUALITY__P2 2.0
#define FXAA_QUALITY__P3 2.0
#define FXAA_QUALITY__P4 2.0
#define FXAA_QUALITY__P5 2.0
#define FXAA_QUALITY__P6 2.0
#define FXAA_QUALITY__P7 4.0
#define FXAA_QUALITY__P8 8.0

float FxaaLuma(float4 rgba) { return rgba.y; }

float4 FxaaPixelShader(float4 rgbyM, float2 RcpFrame, float Subpix, float EdgeThreshold, float EdgeThresholdMin)
{
	float2 posM = GetCoordinates();
	float lumaM = FxaaLuma(rgbyM);
	float lumaS = FxaaLuma(SampleOffset(int2(0, 1)));
	float lumaE = FxaaLuma(SampleOffset(int2(1, 0)));
	float lumaN = FxaaLuma(SampleOffset(int2(0, -1)));
	float lumaW = FxaaLuma(SampleOffset(int2(-1, 0)));

	float maxSM = max(lumaS, lumaM);
	float minSM = min(lumaS, lumaM);
	float maxESM = max(lumaE, maxSM);
	float minESM = min(lumaE, minSM);
	float maxWN = max(lumaN, lumaW);
	float minWN = min(lumaN, lumaW);
	float rangeMax = max(maxWN, maxESM);
	float rangeMin = min(minWN, minESM);
	float rangeMaxScaled = rangeMax * EdgeThreshold;
	float range = rangeMax - rangeMin;
	float rangeMaxClamped = max(EdgeThresholdMin, rangeMaxScaled);
	bool earlyExit = range < rangeMaxClamped;

	if (earlyExit)
		return rgbyM;

	float lumaNW = FxaaLuma(SampleOffset(int2(-1, -1)));
	float lumaSE = FxaaLuma(SampleOffset(int2(1, 1)));
	float lumaNE = FxaaLuma(SampleOffset(int2(1, -1)));
	float lumaSW = FxaaLuma(SampleOffset(int2(-1, 1)));

	float lumaNS = lumaN + lumaS;
	float lumaWE = lumaW + lumaE;
	float subpixRcpRange = 1.0 / range;
	float subpixNSWE = lumaNS + lumaWE;
	float edgeHorz1 = (-2.0 * lumaM) + lumaNS;
	float edgeVert1 = (-2.0 * lumaM) + lumaWE;

	float lumaNESE = lumaNE + lumaSE;
	float lumaNWNE = lumaNW + lumaNE;
	float edgeHorz2 = (-2.0 * lumaE) + lumaNESE;
	float edgeVert2 = (-2.0 * lumaN) + lumaNWNE;

	float lumaNWSW = lumaNW + lumaSW;
	float lumaSWSE = lumaSW + lumaSE;
	float edgeHorz4 = (abs(edgeHorz1) * 2.0) + abs(edgeHorz2);
	float edgeVert4 = (abs(edgeVert1) * 2.0) + abs(edgeVert2);
	float edgeHorz3 = (-2.0 * lumaW) + lumaNWSW;
	float edgeVert3 = (-2.0 * lumaS) + lumaSWSE;
	float edgeHorz = abs(edgeHorz3) + edgeHorz4;
	float edgeVert = abs(edgeVert3) + edgeVert4;

	float subpixNWSWNESE = lumaNWSW + lumaNESE;
	float lengthSign = RcpFrame.x;
	bool horzSpan = edgeHorz >= edgeVert;
	float subpixA = subpixNSWE * 2.0 + subpixNWSWNESE;

	if (!horzSpan) lumaN = lumaW;
	if (!horzSpan) lumaS = lumaE;
	if (horzSpan) lengthSign = RcpFrame.y;
	float subpixB = (subpixA * (1.0 / 12.0)) - lumaM;

	float gradientN = lumaN - lumaM;
	float gradientS = lumaS - lumaM;
	float lumaNN = lumaN + lumaM;
	float lumaSS = lumaS + lumaM;
	bool pairN = abs(gradientN) >= abs(gradientS);
	float gradient = max(abs(gradientN), abs(gradientS));
	if (pairN) lengthSign = -lengthSign;
	float subpixC = saturate(abs(subpixB) * subpixRcpRange);

	float2 posB;
	posB.x = posM.x;
	posB.y = posM.y;
	float2 offNP;
	offNP.x = (!horzSpan) ? 0.0 : RcpFrame.x;
	offNP.y = (horzSpan) ? 0.0 : RcpFrame.y;
	if (!horzSpan) posB.x += lengthSign * 0.5;
	if (horzSpan) posB.y += lengthSign * 0.5;

	float2 posN;
	posN.x = posB.x - offNP.x * FXAA_QUALITY__P0;
	posN.y = posB.y - offNP.y * FXAA_QUALITY__P0;
	float2 posP;
	posP.x = posB.x + offNP.x * FXAA_QUALITY__P0;
	posP.y = posB.y + offNP.y * FXAA_QUALITY__P0;
	float subpixD = ((-2.0)*subpixC) + 3.0;
	float lumaEndN = FxaaLuma(SampleLocation(posN));
	float subpixE = subpixC * subpixC;
	float lumaEndP = FxaaLuma(SampleLocation(posP));

	if (!pairN) lumaNN = lumaSS;
	float gradientScaled = gradient * 1.0 / 4.0;
	float lumaMM = lumaM - lumaNN * 0.5;
	float subpixF = subpixD * subpixE;
	bool lumaMLTZero = lumaMM < 0.0;

	lumaEndN -= lumaNN * 0.5;
	lumaEndP -= lumaNN * 0.5;
	bool doneN = abs(lumaEndN) >= gradientScaled;
	bool doneP = abs(lumaEndP) >= gradientScaled;
	if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P1;
	if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P1;
	bool doneNP = (!doneN) || (!doneP);
	if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P1;
	if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P1;

	if (doneNP) {
		if (!doneN) lumaEndN = FxaaLuma(SampleLocation(posN.xy));
		if (!doneP) lumaEndP = FxaaLuma(SampleLocation(posP.xy));
		if (!doneN) lumaEndN = lumaEndN - lumaNN * 0.5;
		if (!doneP) lumaEndP = lumaEndP - lumaNN * 0.5;
		doneN = abs(lumaEndN) >= gradientScaled;
		doneP = abs(lumaEndP) >= gradientScaled;
		if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P2;
		if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P2;
		doneNP = (!doneN) || (!doneP);
		if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P2;
		if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P2;

		if (doneNP) {
			if (!doneN) lumaEndN = FxaaLuma(SampleLocation(posN.xy));
			if (!doneP) lumaEndP = FxaaLuma(SampleLocation(posP.xy));
			if (!doneN) lumaEndN = lumaEndN - lumaNN * 0.5;
			if (!doneP) lumaEndP = lumaEndP - lumaNN * 0.5;
			doneN = abs(lumaEndN) >= gradientScaled;
			doneP = abs(lumaEndP) >= gradientScaled;
			if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P3;
			if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P3;
			doneNP = (!doneN) || (!doneP);
			if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P3;
			if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P3;

			if (doneNP) {
				if (!doneN) lumaEndN = FxaaLuma(SampleLocation(posN.xy));
				if (!doneP) lumaEndP = FxaaLuma(SampleLocation(posP.xy));
				if (!doneN) lumaEndN = lumaEndN - lumaNN * 0.5;
				if (!doneP) lumaEndP = lumaEndP - lumaNN * 0.5;
				doneN = abs(lumaEndN) >= gradientScaled;
				doneP = abs(lumaEndP) >= gradientScaled;
				if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P4;
				if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P4;
				doneNP = (!doneN) || (!doneP);
				if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P4;
				if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P4;

				if (doneNP) {
					if (!doneN) lumaEndN = FxaaLuma(SampleLocation(posN.xy));
					if (!doneP) lumaEndP = FxaaLuma(SampleLocation(posP.xy));
					if (!doneN) lumaEndN = lumaEndN - lumaNN * 0.5;
					if (!doneP) lumaEndP = lumaEndP - lumaNN * 0.5;
					doneN = abs(lumaEndN) >= gradientScaled;
					doneP = abs(lumaEndP) >= gradientScaled;
					if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P5;
					if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P5;
					doneNP = (!doneN) || (!doneP);
					if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P5;
					if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P5;

					if (doneNP) {
						if (!doneN) lumaEndN = FxaaLuma(SampleLocation(posN.xy));
						if (!doneP) lumaEndP = FxaaLuma(SampleLocation(posP.xy));
						if (!doneN) lumaEndN = lumaEndN - lumaNN * 0.5;
						if (!doneP) lumaEndP = lumaEndP - lumaNN * 0.5;
						doneN = abs(lumaEndN) >= gradientScaled;
						doneP = abs(lumaEndP) >= gradientScaled;
						if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P6;
						if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P6;
						doneNP = (!doneN) || (!doneP);
						if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P6;
						if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P6;

						if (doneNP) {
							if (!doneN) lumaEndN = FxaaLuma(SampleLocation(posN.xy));
							if (!doneP) lumaEndP = FxaaLuma(SampleLocation(posP.xy));
							if (!doneN) lumaEndN = lumaEndN - lumaNN * 0.5;
							if (!doneP) lumaEndP = lumaEndP - lumaNN * 0.5;
							doneN = abs(lumaEndN) >= gradientScaled;
							doneP = abs(lumaEndP) >= gradientScaled;
							if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P7;
							if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P7;
							doneNP = (!doneN) || (!doneP);
							if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P7;
							if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P7;

							if (doneNP) {
								if (!doneN) lumaEndN = FxaaLuma(SampleLocation(posN.xy));
								if (!doneP) lumaEndP = FxaaLuma(SampleLocation(posP.xy));
								if (!doneN) lumaEndN = lumaEndN - lumaNN * 0.5;
								if (!doneP) lumaEndP = lumaEndP - lumaNN * 0.5;
								doneN = abs(lumaEndN) >= gradientScaled;
								doneP = abs(lumaEndP) >= gradientScaled;
								if (!doneN) posN.x -= offNP.x * FXAA_QUALITY__P8;
								if (!doneN) posN.y -= offNP.y * FXAA_QUALITY__P8;
								doneNP = (!doneN) || (!doneP);
								if (!doneP) posP.x += offNP.x * FXAA_QUALITY__P8;
								if (!doneP) posP.y += offNP.y * FXAA_QUALITY__P8;
							}
						}
					}
				}
			}
		}
	}

	float dstN = posM.x - posN.x;
	float dstP = posP.x - posM.x;
	if (!horzSpan) dstN = posM.y - posN.y;
	if (!horzSpan) dstP = posP.y - posM.y;

	bool goodSpanN = (lumaEndN < 0.0) != lumaMLTZero;
	float spanLength = (dstP + dstN);
	bool goodSpanP = (lumaEndP < 0.0) != lumaMLTZero;
	float spanLengthRcp = 1.0 / spanLength;

	bool directionN = dstN < dstP;
	float dst = min(dstN, dstP);
	bool goodSpan = directionN ? goodSpanN : goodSpanP;
	float subpixG = subpixF * subpixF;
	float pixelOffset = (dst * (-spanLengthRcp)) + 0.5;
	float subpixH = subpixG * Subpix;

	float pixelOffsetGood = goodSpan ? pixelOffset : 0.0;
	float pixelOffsetSubpix = max(pixelOffsetGood, subpixH);
	if (!horzSpan) posM.x += pixelOffsetSubpix * lengthSign;
	if (horzSpan) posM.y += pixelOffsetSubpix * lengthSign;

	if (GetOption(C_FXAA_SHOW_EDGES) == 1)
	{
		return -rgbyM;
	}
	else
	{
		return float4(SampleLocation(posM).xyz, lumaM);
	}
}

float4 FxaaPass(float4 color)
{
	return FxaaPixelShader(color, GetInvResolution(), GetOption(A_FXAA_SUBPIX_MAX), GetOption(B_FXAA_EDGE_THRESHOLD), 0.000);
}

//------------------------------------------------------------------------------
// BLOOM
//------------------------------------------------------------------------------

float4 Gauss1dPrev(float2 location, float2 baseoffset, float resolutionmultiplier)
{
	const float offset[] = { 0, 1.4, 4459.0 / 1365.0, 539.0 / 105.0 };
	const float weight[] = { 0.20947265625, 0.30548095703125, 0.08331298828125, 0.00640869140625 };
	float4 Color = SamplePrevLocation(location) * weight[0];	
	baseoffset *= GetInvResolution() * resolutionmultiplier;
	for (int i = 1; i < 4; i++)
	{
		float4 color0 = SamplePrevLocation(location + offset[i] * baseoffset);
		float4 color1 = SamplePrevLocation(location - offset[i] * baseoffset);
		Color += (color0 + color1) * weight[i];
	}
	return Color;
}

void A_ReduceSize()
{
	float3 power = float3(1, 1, 1) * GetOption(B_BLOOMPOWER);
	SetOutput(float4(pow(SamplePrev().rgb, power), 1.0));
}

void BloomH()
{
	float2 texcoord = GetCoordinates();
	SetOutput(Gauss1dPrev(texcoord, float2(1.0, 0.0), GetOption(A_BLOOMWIDTH)));
}

void BloomV()
{
	float2 texcoord = GetCoordinates();
	SetOutput(Gauss1dPrev(texcoord, float2(0.0, 1.0), GetOption(A_BLOOMWIDTH)));
}

void BloomScatering()
{
	float2 SamplePos[20] = {
		float2(0.25, 0.125),
		float2(0.375, 0.125),
		float2(0.5, 0.125),
		float2(0.625, 0.125),
		float2(0.75, 0.125),

		float2(0.25, 0.25),
		float2(0.375, 0.25),
		float2(0.5, 0.25),
		float2(0.625, 0.25),
		float2(0.75, 0.25),

		float2(0.25, 0.75),
		float2(0.375, 0.75),
		float2(0.5, 0.75),
		float2(0.625, 0.75),
		float2(0.75, 0.75),

		float2(0.25, 0.875),
		float2(0.375, 0.875),
		float2(0.5, 0.875),
		float2(0.625, 0.875),
		float2(0.75, 0.875),
	};
	float3 lumColor = GetOption(H_SCOLOR) * 3.0;
	float samplecount = 3.0;
	for (int i = 0; i < 20; i++)
	{
		float3 color = SamplePrevLocation(SamplePos[i]).rgb;
		lumColor += color;
		samplecount += 1.0;
	}
	lumColor /= samplecount;
	float luma = dot(lumColor, lumCoeff);
	float maxval = max(max(lumColor.r, lumColor.g), lumColor.b);
	lumColor /= maxval;
	SetOutput(float4(lumColor, luma * luma));
}

void BloomMerger()
{
	float4 lumColor = SampleInputLocation(2, float2(0.5, 0.5));
	float3 blur = SampleInputBicubic(1).rgb * GetOption(C_BLOOMINTENSITY);
	blur.rgb = blur.rgb * (1.0 - saturate(lumColor.a * 2.0));
	float3 basecolor = float3(0.0,0.0,0.0);
	if (!OptionEnabled(A_BLOOMONLY))
	{
		basecolor = SampleInput(0).rgb;
	}
	
	if (OptionEnabled(D_SCATTERRING))
	{
		float depth = SampleDepth();
		float linearcomponent = (GetOption(G_SEND) - depth) / (GetOption(G_SEND) - GetOption(F_SSTART));
		depth = depth * GetOption(E_SDENSITY);
		lumColor.rgb = lerp(basecolor, lumColor.rgb, saturate(GetOption(I_SINTENSITY) + lumColor.a));
		basecolor = lerp(lumColor.rgb, basecolor, clamp(linearcomponent / exp(depth * depth), 0.0, 1.0));
	}
	SetOutput(float4(basecolor + blur.rgb, 1.0));
}


/*------------------------------------------------------------------------------
[DolphinFX GLOBALS|FUNCTIONS]
------------------------------------------------------------------------------*/

#define FIX(c) max(abs(c), 1e-5)

#define Epsilon (1e-10)
#define lumCoeffDFX float3(0.2126729, 0.7151522, 0.0721750)

float ColorLuminance(float3 color)
{
	return dot(color, lumCoeffDFX);
}

//Average relative luminance
float AvgLuminanceDFX(float3 color)
{
	return sqrt(dot(color * color, lumCoeffDFX));
}

float smootherstep(float a, float b, float x)
{
	x = saturate((x - a) / (b - a));
	return x*x*x*(x*(x * 6.0 - 15.0) + 10.0);
}

/*
float4 DebugClipping(float4 color)
{
if (color.x >= 0.99999 && color.y >= 0.99999 &&
color.z >= 0.99999) color.xyz = float3(1.0f, 0.0f, 0.0f);

if (color.x <= 0.00001 && color.y <= 0.00001 &&
color.z <= 0.00001) color.xyz = float3(0.0f, 0.0f, 1.0f);

return color;
}
*/

//Conversion matrices
float3 RGBtoXYZDFX(float3 rgb)
{
	const float3x3 m = float3x3(
		0.4124564, 0.3575761, 0.1804375,
		0.2126729, 0.7151522, 0.0721750,
		0.0193339, 0.1191920, 0.9503041);

	return mul(rgb, m);
}

float3 XYZtoRGBDFX(float3 xyz)
{
	const float3x3 m = float3x3(
		3.2404542, -1.5371385, -0.4985314,
		-0.9692660, 1.8760108, 0.0415560,
		0.0556434, -0.2040259, 1.0572252);

	return mul(xyz, m);
}

float3 RGBtoYUVDFX(float3 RGB)
{
	const float3x3 m = float3x3(
		0.2126, 0.7152, 0.0722,
		-0.09991, -0.33609, 0.436,
		0.615, -0.55861, -0.05639);

	return mul(RGB, m);
}

float3 YUVtoRGBDFX(float3 YUV)
{
	const float3x3 m = float3x3(
		1.000, 0.000, 1.28033,
		1.000, -0.21482, -0.38059,
		1.000, 2.12798, 0.000);

	return mul(YUV, m);
}

//Converting XYZ to Yxy
float3 XYZtoYxyDFX(float3 xyz)
{
	float3 Yxy;
	float w = 1.0 / (xyz.r + xyz.g + xyz.b);

	Yxy.r = xyz.g;
	Yxy.g = xyz.r * w;
	Yxy.b = xyz.g * w;

	return Yxy;
}

//Converting Yxy to XYZ
float3 YxytoXYZDFX(float3 Yxy)
{
	float3 xyz;
	float w = 1.0 / Yxy.b;
	xyz.g = Yxy.r;
	xyz.r = Yxy.r * Yxy.g * w;
	xyz.b = Yxy.r * (1.0 - Yxy.g - Yxy.b) * w;
	return xyz;
}

float MidLuminance(float3 color)
{
	return sqrt(
		(color.x * color.x * 0.3333) +
		(color.y * color.y * 0.3333) +
		(color.z * color.z * 0.3333));
}


/*------------------------------------------------------------------------------
[COLOUR CORRECTION CODE SECTION]
------------------------------------------------------------------------------*/

// Converting pure hue to RGB
float3 HUEtoRGB(float H)
{
	float R = abs(H * 6.0 - 3.0) - 1.0;
	float G = 2.0 - abs(H * 6.0 - 2.0);
	float B = 2.0 - abs(H * 6.0 - 4.0);

	return saturate(float3(R, G, B));
}

// Converting RGB to hue/chroma/value
float3 RGBtoHCV(float3 RGB)
{
	float4 BG = float4(RGB.bg, -1.0, 2.0 / 3.0);
	float4 GB = float4(RGB.gb, 0.0, -1.0 / 3.0);

	float4 P = (RGB.g < RGB.b) ? BG : GB;

	float4 XY = float4(P.xyw, RGB.r);
	float4 YZ = float4(RGB.r, P.yzx);

	float4 Q = (RGB.r < P.x) ? XY : YZ;

	float C = Q.x - min(Q.w, Q.y);
	float H = abs((Q.w - Q.y) / (6.0 * C + Epsilon) + Q.z);

	return float3(H, C, Q.x);
}

// Converting RGB to HSV
float3 RGBtoHSV(float3 RGB)
{
	float3 HCV = RGBtoHCV(RGB);
	float S = HCV.y / (HCV.z + Epsilon);

	return float3(HCV.x, S, HCV.z);
}

// Converting HSV to RGB
float3 HSVtoRGB(float3 HSV)
{
	float3 RGB = HUEtoRGB(HSV.x);
	return ((RGB - 1.0) * HSV.y + 1.0) * HSV.z;
}

// Pre correction color mask
float3 PreCorrection(float3 color)
{
	float3 RGB = color;

	RGB.r = 2.0 / 3.0 * (1.0 - (RGB.r * RGB.r));
	RGB.g = 2.0 / 3.0 * (1.0 - (RGB.g * RGB.g));
	RGB.b = 2.0 / 3.0 * (1.0 - (RGB.b * RGB.b));

	RGB.r = saturate(color.r + (GetOption(B_RED_CORRECTION) / 200.0) * RGB.r);
	RGB.g = saturate(color.g + (GetOption(C_GREEN_CORRECTION) / 200.0) * RGB.g);
	RGB.b = saturate(color.b + (GetOption(D_BLUE_CORRECTION) / 200.0) * RGB.b);

	color = saturate(RGB);

	return color;
}

float3 ColorCorrection(float3 color)
{
	float X = 1.0 / (1.0 + exp(GetOption(B_RED_CORRECTION) / 2.0));
	float Y = 1.0 / (1.0 + exp(GetOption(C_GREEN_CORRECTION) / 2.0));
	float Z = 1.0 / (1.0 + exp(GetOption(D_BLUE_CORRECTION) / 2.0));

	color.r = (1.0 / (1.0 + exp(-GetOption(B_RED_CORRECTION) * (color.r - 0.5))) - X) / (1.0 - 2.0 * X);
	color.g = (1.0 / (1.0 + exp(-GetOption(C_GREEN_CORRECTION) * (color.g - 0.5))) - Y) / (1.0 - 2.0 * Y);
	color.b = (1.0 / (1.0 + exp(-GetOption(D_BLUE_CORRECTION) * (color.b - 0.5))) - Z) / (1.0 - 2.0 * Z);

	return saturate(color);
}

float4 CorrectionPass(float4 color)
{
	float3 colorspace = PreCorrection(color.rgb);

	if (GetOption(A_PALETTE) == 0) {
		colorspace = ColorCorrection(colorspace);
	}

	else if (GetOption(A_PALETTE) == 1) {
		float3 XYZ = RGBtoXYZDFX(colorspace);
		float3 Yxy = XYZtoYxyDFX(XYZ);

		Yxy = ColorCorrection(Yxy);
		XYZ = YxytoXYZDFX(Yxy);
		colorspace = XYZtoRGBDFX(XYZ);
	}

	else if (GetOption(A_PALETTE) == 2) {
		float3 XYZ = RGBtoXYZDFX(colorspace);
		float3 Yxy = XYZtoYxyDFX(XYZ);

		XYZ = YxytoXYZDFX(Yxy);
		XYZ = ColorCorrection(XYZ);
		colorspace = XYZtoRGBDFX(XYZ);
	}

	else if (GetOption(A_PALETTE) == 3) {
		float3 hsv = RGBtoHSV(colorspace);
		hsv = ColorCorrection(hsv);
		colorspace = HSVtoRGB(hsv);
	}

	else if (GetOption(A_PALETTE) == 4) {
		float3 yuv = RGBtoYUVDFX(colorspace);
		yuv = ColorCorrection(yuv);
		colorspace = YUVtoRGBDFX(yuv);
	}

	color.rgb = lerp(color.rgb, colorspace, GetOption(E_CORRECT_STR));
	color.a = AvgLuminanceDFX(color.rgb);

	return color;
}

/*------------------------------------------------------------------------------
[PAINT SHADING CODE SECTION]
------------------------------------------------------------------------------*/

float3 PaintShading(float3 color, float2 texcoord)
{
	float2 pixelSize = GetInvResolution();
	if (GetOption(PaintMethod) == 1)
	{
		float2	tex;
		int	k, j, lum, cmax = 0;

		float	C0 = 0, C1 = 0, C2 = 0, C3 = 0, C4 = 0, C5 = 0, C6 = 0, C7 = 0, C8 = 0, C9 = 0;
		float3	A = float3(0.0,0.0,0.0), B = float3(0.0, 0.0, 0.0), C = float3(0.0, 0.0, 0.0), D = float3(0.0, 0.0, 0.0), E = float3(0.0, 0.0, 0.0), F = float3(0.0, 0.0, 0.0), G = float3(0.0, 0.0, 0.0), H = float3(0.0, 0.0, 0.0), I = float3(0.0, 0.0, 0.0), J = float3(0.0, 0.0, 0.0), shade = float3(0.0, 0.0, 0.0);

		for (k = int(-PaintRadius); k < (int(PaintRadius) + 1); k++) {
			for (j = int(-PaintRadius); j < (int(PaintRadius) + 1); j++) {

				tex.x = texcoord.x + pixelSize.x * k;
				tex.y = texcoord.y + pixelSize.y * j;

				shade = SampleLocation(tex).xyz;

				lum = int(AvgLuminanceDFX(shade) * 9.0);

				C0 = (lum == 0) ? C0 + 1 : C0;
				C1 = (lum == 1) ? C1 + 1 : C1;
				C2 = (lum == 2) ? C2 + 1 : C2;
				C3 = (lum == 3) ? C3 + 1 : C3;
				C4 = (lum == 4) ? C4 + 1 : C4;
				C5 = (lum == 5) ? C5 + 1 : C5;
				C6 = (lum == 6) ? C6 + 1 : C6;
				C7 = (lum == 7) ? C7 + 1 : C7;
				C8 = (lum == 8) ? C8 + 1 : C8;
				C9 = (lum == 9) ? C9 + 1 : C9;

				A = (lum == 0) ? A + shade : A;
				B = (lum == 1) ? B + shade : B;
				C = (lum == 2) ? C + shade : C;
				D = (lum == 3) ? D + shade : D;
				E = (lum == 4) ? E + shade : E;
				F = (lum == 5) ? F + shade : F;
				G = (lum == 6) ? G + shade : G;
				H = (lum == 7) ? H + shade : H;
				I = (lum == 8) ? I + shade : I;
				J = (lum == 9) ? J + shade : J;
			}
		}

		if (C0 > cmax) { cmax = int(C0); color.xyz = A / cmax; }
		if (C1 > cmax) { cmax = int(C1); color.xyz = B / cmax; }
		if (C2 > cmax) { cmax = int(C2); color.xyz = C / cmax; }
		if (C3 > cmax) { cmax = int(C3); color.xyz = D / cmax; }
		if (C4 > cmax) { cmax = int(C4); color.xyz = E / cmax; }
		if (C5 > cmax) { cmax = int(C5); color.xyz = F / cmax; }
		if (C6 > cmax) { cmax = int(C6); color.xyz = G / cmax; }
		if (C7 > cmax) { cmax = int(C7); color.xyz = H / cmax; }
		if (C8 > cmax) { cmax = int(C8); color.xyz = I / cmax; }
		if (C9 > cmax) { cmax = int(C9); color.xyz = J / cmax; }
	}
	else
	{
		int j, i;
		float2 screenSize = GetResolution();
		float3 m0, m1, m2, m3, k0, k1, k2, k3, shade;
		float n = float((PaintRadius + 1.0) * (PaintRadius + 1.0));

		for (j = int(-PaintRadius); j <= 0; ++j) {
			for (i = int(-PaintRadius); i <= 0; ++i) {

				shade = SampleLocation(texcoord + float2(i, j) / screenSize).rgb;
				m0 += shade; k0 += shade * shade;
			}
		}

		for (j = int(-PaintRadius); j <= 0; ++j) {
			for (i = 0; i <= int(PaintRadius); ++i) {
				shade = SampleLocation(texcoord + float2(i, j) / screenSize).rgb;
				m1 += shade; k1 += shade * shade;
			}
		}

		for (j = 0; j <= int(PaintRadius); ++j) {
			for (i = 0; i <= int(PaintRadius); ++i) {
				shade = SampleLocation(texcoord + float2(i, j) / screenSize).rgb;
				m2 += shade; k2 += shade * shade;
			}
		}

		float min_sigma2 = 1e+2;
		m0 /= n; k0 = abs(k0 / n - m0 * m0);

		float sigma2 = k0.r + k0.g + k0.b;
		if (sigma2 < min_sigma2) {
			min_sigma2 = sigma2; color = m0;
		}

		m1 /= n; k1 = abs(k1 / n - m1 * m1);
		sigma2 = k1.r + k1.g + k1.b;

		if (sigma2 < min_sigma2) {
			min_sigma2 = sigma2;
			color = m1;
		}

		m2 /= n; k2 = abs(k2 / n - m2 * m2);
		sigma2 = k2.r + k2.g + k2.b;

		if (sigma2 < min_sigma2) {
			min_sigma2 = sigma2;
			color = m2;
		}
	}
	return color;
}

float4 PaintPass(float4 color, float2 texcoord)
{
	float3 paint = PaintShading(color.rgb, texcoord);
	color.rgb = lerp(color.rgb, paint, GetOption(PaintStrength));
	color.a = AvgLuminanceDFX(color.rgb);

	return color;
}

/*------------------------------------------------------------------------------
[PAINT SHADING 2 CODE SECTION]
------------------------------------------------------------------------------*/

float3 PaintShading2(float3 color, float2 texcoord)
{
	float2 pixelSize = GetInvResolution();
	if (GetOption(PaintMethod2) == 1)
	{
		float2	tex;
		int	k, j, lum, cmax = 0;

		float	C0 = 0, C1 = 0, C2 = 0, C3 = 0, C4 = 0, C5 = 0, C6 = 0, C7 = 0, C8 = 0, C9 = 0;
		float3	A = float3(0.0,0.0,0.0), B = float3(0.0, 0.0, 0.0), C = float3(0.0, 0.0, 0.0), D = float3(0.0, 0.0, 0.0), E = float3(0.0, 0.0, 0.0), F = float3(0.0, 0.0, 0.0), G = float3(0.0, 0.0, 0.0), H = float3(0.0, 0.0, 0.0), I = float3(0.0, 0.0, 0.0), J = float3(0.0, 0.0, 0.0), shade = float3(0.0, 0.0, 0.0);

		for (k = int(-PaintRadius2); k < (int(PaintRadius2) + 1); k++) {
			for (j = int(-PaintRadius2); j < (int(PaintRadius2) + 1); j++) {

				tex.x = texcoord.x + pixelSize.x * k;
				tex.y = texcoord.y + pixelSize.y * j;

				shade = SampleLocation(tex).xyz;

				lum = int(AvgLuminanceDFX(shade) * 9.0);

				C0 = (lum == 0) ? C0 + 1 : C0;
				C1 = (lum == 1) ? C1 + 1 : C1;
				C2 = (lum == 2) ? C2 + 1 : C2;
				C3 = (lum == 3) ? C3 + 1 : C3;
				C4 = (lum == 4) ? C4 + 1 : C4;
				C5 = (lum == 5) ? C5 + 1 : C5;
				C6 = (lum == 6) ? C6 + 1 : C6;
				C7 = (lum == 7) ? C7 + 1 : C7;
				C8 = (lum == 8) ? C8 + 1 : C8;
				C9 = (lum == 9) ? C9 + 1 : C9;

				A = (lum == 0) ? A + shade : A;
				B = (lum == 1) ? B + shade : B;
				C = (lum == 2) ? C + shade : C;
				D = (lum == 3) ? D + shade : D;
				E = (lum == 4) ? E + shade : E;
				F = (lum == 5) ? F + shade : F;
				G = (lum == 6) ? G + shade : G;
				H = (lum == 7) ? H + shade : H;
				I = (lum == 8) ? I + shade : I;
				J = (lum == 9) ? J + shade : J;
			}
		}

		if (C0 > cmax) { cmax = int(C0); color.xyz = A / cmax; }
		if (C1 > cmax) { cmax = int(C1); color.xyz = B / cmax; }
		if (C2 > cmax) { cmax = int(C2); color.xyz = C / cmax; }
		if (C3 > cmax) { cmax = int(C3); color.xyz = D / cmax; }
		if (C4 > cmax) { cmax = int(C4); color.xyz = E / cmax; }
		if (C5 > cmax) { cmax = int(C5); color.xyz = F / cmax; }
		if (C6 > cmax) { cmax = int(C6); color.xyz = G / cmax; }
		if (C7 > cmax) { cmax = int(C7); color.xyz = H / cmax; }
		if (C8 > cmax) { cmax = int(C8); color.xyz = I / cmax; }
		if (C9 > cmax) { cmax = int(C9); color.xyz = J / cmax; }
	}
	else
	{
		int j, i;
		float2 screenSize = GetResolution();
		float3 m0, m1, m2, m3, k0, k1, k2, k3, shade;
		float n = float((PaintRadius2 + 1.0) * (PaintRadius2 + 1.0));

		for (j = int(-PaintRadius2); j <= 0; ++j) {
			for (i = int(-PaintRadius2); i <= 0; ++i) {

				shade = SampleLocation(texcoord + float2(i, j) / screenSize).rgb;
				m0 += shade; k0 += shade * shade;
			}
		}

		for (j = int(-PaintRadius2); j <= 0; ++j) {
			for (i = 0; i <= int(PaintRadius2); ++i) {
				shade = SampleLocation(texcoord + float2(i, j) / screenSize).rgb;
				m1 += shade; k1 += shade * shade;
			}
		}

		for (j = 0; j <= int(PaintRadius2); ++j) {
			for (i = 0; i <= int(PaintRadius2); ++i) {
				shade = SampleLocation(texcoord + float2(i, j) / screenSize).rgb;
				m2 += shade; k2 += shade * shade;
			}
		}

		float min_sigma2 = 1e+2;
		m0 /= n; k0 = abs(k0 / n - m0 * m0);

		float sigma2 = k0.r + k0.g + k0.b;
		if (sigma2 < min_sigma2) {
			min_sigma2 = sigma2; color = m0;
		}

		m1 /= n; k1 = abs(k1 / n - m1 * m1);
		sigma2 = k1.r + k1.g + k1.b;

		if (sigma2 < min_sigma2) {
			min_sigma2 = sigma2;
			color = m1;
		}

		m2 /= n; k2 = abs(k2 / n - m2 * m2);
		sigma2 = k2.r + k2.g + k2.b;

		if (sigma2 < min_sigma2) {
			min_sigma2 = sigma2;
			color = m2;
		}
	}
	return color;
}

float4 PaintPass2(float4 color, float2 texcoord)
{
	float3 paint = PaintShading2(color.rgb, texcoord);
	color.rgb = lerp(color.rgb, paint, GetOption(PaintStrength2));
	color.a = AvgLuminanceDFX(color.rgb);

	return color;
}

/*------------------------------------------------------------------------------
[DOLPHINFX INTEGRATION FUNCTIONS]
------------------------------------------------------------------------------*/

void PostApplyCC()
{
	SetOutput(CorrectionPass(float4(SamplePrev())));
}

void SceneTonemapping()
{
	SetOutput(TonemapPass(float4(SamplePrev())));
}

void PaintShader()
{
	SetOutput(PaintPass(float4(SamplePrev()), GetCoordinates()));
}

void PaintShader2()
{
	SetOutput(PaintPass2(float4(SamplePrev()), GetCoordinates()));
}

void PixelVibrance()
{
	SetOutput(VibrancePass(float4(SamplePrev())));
}

void Fxaa()
{
	SetOutput(FxaaPass(float4(SamplePrev())));
}