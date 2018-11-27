/*===============================================================================*\
|########################        [Ishiiruka FX 0.8]        ######################||
|| Credist to:                                                                   ||
|| Asmodean (DolphinFX)                                                          ||
|| Matso (MATSODOF)                                                              ||
|| Gilcher Pascal aka Marty McFly (MATSODOF original port to MCFX)               ||
|| Daniel R�kos (Efficient Gaussian blur with linear sampling)                   ||
|| mudlord (FXAA)                                                                ||
|| Cozimo (Default shader settings for Twilight Princess Mod presets)            ||
|| Drakonas (Integrated DolphinFX effects for Twilight Princess Mod Team)        ||
|| Drakonas (With Tino's help, also integrated ReShade's Gaussian AnamFlare)     ||
|############################        By Tino          ############################|
\*===============================================================================*/
/*
[configuration]
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
DefaultValue = 5.7
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Bloom Intensity
OptionName = C_BLOOMINTENSITY
MinValue = 0.1
MaxValue = 1.0
StepAmount = 0.01
DefaultValue = 0.35
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
DefaultValue = 0.14
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Start
OptionName = F_SSTART
MinValue = 0.0
MaxValue = 0.5
StepAmount = 0.01
DefaultValue = 0.16
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
DefaultValue = 0.82, 0.81, 0.91
DependentOption = C_BLOOM

[OptionRangeFloat]
GUIName = Scattering Intensity
OptionName = I_SINTENSITY
MinValue = 0.0
MaxValue = 1.0
StepAmount = 0.01
DefaultValue = 0.33
DependentOption = C_BLOOM

[OptionBool]
GUIName = Gaussian AnamFlare
OptionName = C_GAUSSIAN_ANAMFLARE
DefaultValue = True
GUIDescription = Enable to apply a horizontal light beam to bright pixels.

[OptionRangeFloat]
GUIName = Threshhold
OptionName = A_ANAMFLARE_THRESHOLD
MinValue = 0.10
MaxValue = 1.00
StepAmount = 0.01
DefaultValue = 0.90
DependentOption = C_GAUSSIAN_ANAMFLARE
GUIDescription = Every pixel brighter than this value gets a flare.

[OptionRangeFloat]
GUIName = Wideness
OptionName = A_ANAMFLARE_WIDENESS
MinValue = 0.5
MaxValue = 2.0
StepAmount = 0.1
DefaultValue = 1.0
DependentOption = C_GAUSSIAN_ANAMFLARE
GUIDescription = Horizontal wideness of flare. Don't set too high, otherwise the single samples are visible

[OptionRangeFloat]
GUIName = Amount
OptionName = B_ANAMFLARE_AMOUNT
MinValue = 1.0
MaxValue = 20.0
StepAmount = 0.5
DefaultValue = 2.5
DependentOption = C_GAUSSIAN_ANAMFLARE
GUIDescription = Intensity of anamorphic flare.

[OptionRangeFloat]
GUIName = Curve
OptionName = B_ANAMFLARE_CURVE
MinValue = 1.0
MaxValue = 2.0
StepAmount = 0.1
DefaultValue = 1.0
DependentOption = C_GAUSSIAN_ANAMFLARE
GUIDescription = Intensity curve of flare with distance from source

[OptionRangeFloat]
GUIName = Colors
OptionName = C_ANAMFLARE_COLOR
MinValue = 0.0, 0.0, 0.0
MaxValue = 1.0, 1.0, 1.0
StepAmount = 0.01, 0.01, 0.01
DefaultValue = 0.010, 0.160, 0.580
DependentOption = C_GAUSSIAN_ANAMFLARE
GUIDescription = R, G and B components of anamorphic flare. Flare is always same color.

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
DefaultValue = 0.5
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = FilmStrength
OptionName = B_TONE_FAMOUNT
MinValue = 0.00
MaxValue = 1.00
StepAmount = 0.01
DefaultValue = 0.19
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = BlackLevels
OptionName = C_BLACK_LEVELS
MinValue = 0.00
MaxValue = 1.00
StepAmount = 0.01
DefaultValue = 0.06
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = Exposure
OptionName = D_EXPOSURE
MinValue = 0.50
MaxValue = 1.50
StepAmount = 0.01
DefaultValue = 1.02
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = Luminance
OptionName = E_LUMINANCE
MinValue = 0.50
MaxValue = 1.50
StepAmount = 0.01
DefaultValue = 1.00
DependentOption = C_TONEMAP_PASS

[OptionRangeFloat]
GUIName = WhitePoint
OptionName = F_WHITEPOINT
MinValue = 0.50
MaxValue = 1.50
StepAmount = 0.01
DefaultValue = 0.95
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
DefaultValue = 2.1
DependentOption = D_COLOR_CORRECTION

[OptionRangeFloat]
GUIName = Channels G|X|Y|S|U
OptionName = C_GREEN_CORRECTION
MinValue = 0.10
MaxValue = 8.00
StepAmount = 0.10
DefaultValue = 1.4
DependentOption = D_COLOR_CORRECTION

[OptionRangeFloat]
GUIName = Channels B|Y|Z|V|V
OptionName = D_BLUE_CORRECTION
MinValue = 0.10
MaxValue = 8.00
StepAmount = 0.10
DefaultValue = 0.1
DependentOption = D_COLOR_CORRECTION

[OptionRangeFloat]
GUIName = PaletteStrength
OptionName = E_CORRECT_STR
MinValue = 0.00
MaxValue = 2.00
StepAmount = 0.01
DefaultValue = 0.2
DependentOption = D_COLOR_CORRECTION

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
DefaultValue = 0.14
DependentOption = F_PIXEL_VIBRANCE

[OptionRangeFloat]
GUIName = RedVibrance
OptionName = B_R_VIBRANCE
MinValue = -1.00
MaxValue = 4.00
StepAmount = 0.01
DefaultValue = 1.00
DependentOption = F_PIXEL_VIBRANCE

[OptionRangeFloat]
GUIName = GreenVibrance
OptionName = C_G_VIBRANCE
MinValue = -1.00
MaxValue = 4.00
StepAmount = 0.01
DefaultValue = 1.00
DependentOption = F_PIXEL_VIBRANCE

[OptionRangeFloat]
GUIName = BlueVibrance
OptionName = D_B_VIBRANCE
MinValue = -1.00
MaxValue = 4.00
StepAmount = 0.01
DefaultValue = 1.00
DependentOption = F_PIXEL_VIBRANCE

[OptionBool]
GUIName = Vignette
OptionName = H_VIGNETTE_PASS
DefaultValue = True

[OptionRangeFloat]
GUIName = VignetteRatio
OptionName = A_VIG_RATIO
MinValue = 1.00
MaxValue = 2.00
StepAmount = 0.01
DefaultValue = 1.49
DependentOption = H_VIGNETTE_PASS

[OptionRangeFloat]
GUIName = VignetteRadius
OptionName = B_VIG_RADIUS
MinValue = 0.50
MaxValue = 2.00
StepAmount = 0.01
DefaultValue = 1.07
DependentOption = H_VIGNETTE_PASS

[OptionRangeFloat]
GUIName = VignetteAmount
OptionName = C_VIG_AMOUNT
MinValue = 0.10
MaxValue = 2.00
StepAmount = 0.01
DefaultValue = 0.1
DependentOption = H_VIGNETTE_PASS

[OptionRangeInteger]
GUIName = VignetteSlope
OptionName = D_VIG_SLOPE
MinValue = 2
MaxValue = 16
StepAmount = 2
DefaultValue = 14
DependentOption = H_VIGNETTE_PASS

[Pass]
EntryPoint = PreApplyCC
DependantOption = D_COLOR_CORRECTION
Input0=ColorBuffer
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
EntryPoint = BloomVwithAnamFlare
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
Input0=Pass0
Input0Filter=Linear
Input0Mode=Clamp
Input1=Pass5
Input1Filter=Nearest
Input1Mode=Clamp
Input2=Pass6
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
EntryPoint = PixelVibrance
DependantOption = F_PIXEL_VIBRANCE
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[Pass]
EntryPoint = Vignette
DependantOption = H_VIGNETTE_PASS
Input0=PreviousPass
Input0Filter=Linear
Input0Mode=Clamp
[/configuration]
*/

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
    float4 rawcolor = SamplePrev();
    float4 reducendcolor = float4(pow(rawcolor.rgb, power), 1.0);
    if(OptionEnabled(C_GAUSSIAN_ANAMFLARE))
    {
        reducendcolor.w = max(0,dot(rawcolor.xyz,float3(0.333,0.333,0.333))-GetOption(A_ANAMFLARE_THRESHOLD));
    }
    SetOutput(reducendcolor);
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
	float4 color;
	color.rgb = basecolor + blur.rgb;

	if (OptionEnabled(C_GAUSSIAN_ANAMFLARE))
	{
		float3 anamflare = SampleInputBicubic(1).w*2*GetOption(C_ANAMFLARE_COLOR);
		color.xyz += pow(anamflare.xyz,1/GetOption(B_ANAMFLARE_CURVE));
	}
	SetOutput(color);
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
[VIGNETTE CODE SECTION]
------------------------------------------------------------------------------*/

float4 VignettePass(float4 color)
{
	const float2 VignetteCenter = float2(0.500, 0.500);
	float2 tc = GetCoordinates() - VignetteCenter;

	// hardcoded pre ratio calculations, for uniform output on arbitrary resolutions.
	tc *= float2((2560.0 / 1440.0), GetOption(A_VIG_RATIO));
	tc /= GetOption(B_VIG_RADIUS);

	float v = dot(tc, tc);

	color.rgb *= (1.0 + pow(v, GetOption(D_VIG_SLOPE) * 0.25) * -GetOption(C_VIG_AMOUNT));

	return color;
}

/*------------------------------------------------------------------------------
[DOLPHINFX INTEGRATION FUNCTIONS]
------------------------------------------------------------------------------*/

float4 Gauss1dPrev2(float2 coord, float2 axis, float mult)
{
    float4 sum = 0;
    float weight[11] = { 0.082607, 0.080977, 0.076276, 0.069041, 0.060049, 0.050187, 0.040306, 0.031105, 0.023066, 0.016436, 0.011254 };
    axis *= GetInvResolution() * mult;
    for (int i = -10; i < 11; i++)
    {
        float currweight = weight[abs(i)];
        sum += SamplePrevLocation(coord.xy + axis.xy * (float)i) * currweight;
    }
    return sum;
}

void BloomVwithAnamFlare()
{
	float2 texcoord = GetCoordinates();
	float4 bloom = Gauss1dPrev(texcoord, float2(0.0, 1.0), GetOption(A_BLOOMWIDTH));
	if(OptionEnabled(C_GAUSSIAN_ANAMFLARE))
	{
		bloom.a = Gauss1dPrev2(texcoord, float2(1.0, 0.0), 2*GetOption(A_ANAMFLARE_WIDENESS)).a*2.5;
		bloom.a *= GetOption(B_ANAMFLARE_AMOUNT);
	}
	SetOutput(bloom);
}

void PreApplyCC()
{
	SetOutput(CorrectionPass(float4(SamplePrev())));
}

void SceneTonemapping()
{
	SetOutput(TonemapPass(float4(SamplePrev())));
}

void PixelVibrance()
{
	SetOutput(VibrancePass(float4(SamplePrev())));
}

void Vignette()
{
	SetOutput(VignettePass(float4(SamplePrev())));
}
