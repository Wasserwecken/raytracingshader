#version 330

vec2 vec2Zero = vec2(0.0);
vec2 vec2One = vec2(1.0);
vec3 vec3Zero = vec3(0.0);
vec3 vec3One = vec3(1.0);
vec3 vec3Left = vec3(-1.0, 0.0, 0.0);
vec3 vec3Right = vec3(1.0, 0.0, 0.0);
vec3 vec3Up = vec3(0.0, 1.0, 0.0);
vec3 vec3Down = vec3(0.0, -1.0, 0.0);
vec3 vec3Forward = vec3(0.0, 0.0, 1.0);
vec3 vec3Backward = vec3(0.0, 0.0, -1.0);

float PI = 3.14159265359;
float PI2 = 6.28318530717;
float RADPERDEG = 0.017453292519;
float DEGPERRAD = 57.29577951308;
vec3 BACKDROPCOLOR = vec3(0.0941, 0.0941, 0.0941);
vec3 INVALIDCOLOR = vec3(1.0, 0.0, 0.8314);

int MAXSTEPS = 500;
float MAXDISTANCE = 100.0;
float STEPTOLERANCE = 0.001;
float NORMALSAMPLEDISTANCE = 0.001;

vec3 COLORIRON = vec3(0.560, 0.570, 0.580);
vec3 COLORSILVER = vec3(0.972, 0.960, 0.915);
vec3 COLORALUMINIUM = vec3(0.913, 0.921, 0.925);
vec3 COLORGOLD = vec3(1.000, 0.766, 0.336);
vec3 COLORCOPPER = vec3(0.955, 0.637, 0.538);
vec3 COLORCHROM = vec3(0.550, 0.556, 0.554);








//-------------------------------------------
// Structs
//-------------------------------------------
struct Ray
{
	vec3 origin;
	vec3 direction;
};

struct DFInfo
{
    float hitDistance;
    vec3 origin;
    float materialId;
    vec2 uv;
};

struct Trace
{
    DFInfo distanceField;
    vec3 point;
    int steps;
};

struct SurfaceHit
{
    float hitDistance;
    vec3 point;
    vec3 normal;
    int steps;
    vec3 origin;
    float materialId;
    vec2 basicUV;
};

struct Material
{
    vec3 albedo;
    float metallic;
    float roughness;
    float occlusion;
    float height;
};

struct SurfaceColor
{
    vec3 color;
};

struct Light
{
    vec3 color;
    float intensity;
};

struct PointLight
{
    Light light;
    vec3 position;
    float falloff;
    float shadowStrength;
    float shadowSmoothness;
};









//-------------------------------------------
// Random
//-------------------------------------------
float RandomValue(vec2 seed)
{
    return fract(sin(dot(seed, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 RandomVector(vec3 seed) {
    return vec3(
            RandomValue(seed.xz) * 2.0 - 1.0,
            RandomValue(seed.yx) * 2.0 - 1.0,
            RandomValue(seed.zy) * 2.0 - 1.0
    );
}

vec3 RandomSpherePoint(vec3 seed) {
    return normalize(RandomVector(seed));
}

vec3 RandomHemispherePoint(vec3 normal, vec3 seed) {
    vec3 value = RandomSpherePoint(seed);
    value *= sign(dot(value, normal));
    return value;
}

vec3 RandomColor(vec2 seed)
{
    return normalize(vec3(
            RandomValue(seed * 3.0),
            RandomValue(seed * 5.0),
            RandomValue(seed * 7.0)
    ));
}

// Value noise by Inigo Quilez - iq/2013
// https://www.shadertoy.com/view/lsf3WH
float Noise(vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);
    vec2 u = f*f*(3.0-2.0*f);
    return mix( mix( RandomValue( i + vec2(0.0,0.0) ),
                     RandomValue( i + vec2(1.0,0.0) ), u.x),
                mix( RandomValue( i + vec2(0.0,1.0) ),
                     RandomValue( i + vec2(1.0,1.0) ), u.x), u.y);
}








//-------------------------------------------
// VALUE TOOLS
//-------------------------------------------
float ValueLinearStep(float edgeCenter, float edgeWidth, float value)
{
    value -= edgeCenter - (edgeWidth * 0.5);
    value /= max(0.0001, edgeWidth);
    value = clamp(value, 0.0, 1.0);

    return value;
}

vec2 ValueLinearStep(vec2 edgeCenter, vec2 edgeWidth, vec2 value)
{
    value -= edgeCenter - (edgeWidth * 0.5);
    value /= max(vec2(0.0001), edgeWidth);
    value = clamp(value, 0.0, 1.0);

    return value;
}

float EasePowerOut(float power, float value) {
    value = pow(--value, power);

    float corrector = (floor(mod(power + 0.1, 2.0)) * 2.0) - 1.0;
    value = (value + corrector) * corrector;

    return value;
}








//-------------------------------------------
// UV TOOLS
//-------------------------------------------
vec2 UVRotate(vec2 uv, vec2 origin, float angle)
{	
    angle *= RADPERDEG;

    float sinus = sin(angle);
	float cosinus = cos(angle);
	mat2 matrix = mat2(cosinus, -sinus, sinus, cosinus);

    uv -= origin;
    uv = matrix * uv;
    uv += origin;

	return uv;
}

vec2 UVTiling(vec2 uv, vec2 tiles, float offsetX, out vec2 tileIndex)
{
    uv = uv * tiles;
    uv.x += floor(uv).y * offsetX;

    tileIndex = floor(uv);
    return fract(uv);
}

vec2 UVWarp(vec2 uv, float noise, float strength)
{
    vec2 rotationOrigin = uv - vec2(1.0);
    return UVRotate(uv, rotationOrigin, 360.0 * noise * strength);
}

vec2 UVSplitSpaceRadial(vec2 uv, float splits, float originDistance, out float spaceIndex)
{
    float rotationStep = 1.0 / splits;
    float rotationOffset = rotationStep * 0.5;
    float angle = 1.0 - (fract(((atan(uv.y, uv.x) + PI) / PI2) - rotationOffset));

    spaceIndex = floor(angle * splits); 
    float rotation = spaceIndex * rotationStep;
    
    uv = UVRotate(uv, vec2(0.0), -rotation * PI2 * DEGPERRAD);
    uv.x += originDistance;

    return uv;
}









//-------------------------------------------
// SHAPE TOOLS
//-------------------------------------------
float ShapeRectangle(vec2 uv, vec2 size, vec2 border)
{
    size *= 0.5;
    vec2 isRectangle = vec2(1.0);
    isRectangle *= ValueLinearStep(vec2(0.5) - size, vec2(border), uv);
    isRectangle *= 1.0 - ValueLinearStep(vec2(0.5) + size, vec2(border), uv);

    return min(isRectangle.x,  isRectangle.y);
}

float ShapeCircle(vec2 uv, float radius, float blur)
{
    float isCircle = ValueLinearStep(radius, blur, length(uv));
    return 1.0 - isCircle;
}










//-------------------------------------------
// MATERIALS
//-------------------------------------------
Material Tiles(vec2 uv)
{
    vec2 tileIndex;
    vec2 tiledUV = UVTiling(uv, vec2(10.0), 0.0, tileIndex);

    float randomRotation = RandomValue(tileIndex) * 2.0 - 1.0;
    randomRotation *= 0.5;
    tiledUV = UVRotate(tiledUV, vec2(0.5), randomRotation);

    float tiles = ShapeRectangle(tiledUV, vec2(0.9), vec2(0.06));
    float tileMask = tiles;

    vec2 crackUV1 = UVRotate(uv, vec2(0.5), RandomValue(tileIndex * 3.0) * 180.0);
    vec2 crackUV2 = UVRotate(uv, vec2(0.5), RandomValue(tileIndex * 5.0) * 180.0);
    vec2 crackUV3 = UVRotate(uv, vec2(0.5), RandomValue(tileIndex * 7.0) * 180.0);
    float cracks1 = Noise((crackUV1 + tileIndex) * 25.0);
    float cracks2 = Noise(vec2(10.0, 1.0) * (crackUV2 - tileIndex));
    float cracks3 = Noise(vec2(10.0, 1.0) * (crackUV3 + tileIndex));
    float cracks = min(cracks2, cracks3) * cracks1 * cracks1;
    cracks = step(0.2, cracks) * step(cracks, 0.25);

    float probability = 0.8 + ((tileIndex.y + 10.0) * 0.01);
    float crackredTiles = step(probability, RandomValue(tileIndex));
    tiles -= crackredTiles * cracks * tiles;
    float crackedTileMask = tiles;

    float roughness = RandomValue(tileIndex) * 0.01 + 0.8;
    roughness = mix(roughness, 0.5, step(tileIndex.y, -14.0));
    roughness = mix(roughness, 0.5, step(-3.0, tileIndex.y) * step(tileIndex.y, -2.0));
    roughness *= step(0.95, crackedTileMask);
    roughness = 1.0 - roughness;

    vec3 tileAlbedo = vec3(0.9);
    tileAlbedo = mix(tileAlbedo, vec3(0.1), step(tileIndex.y, -14.0));
    tileAlbedo = mix(tileAlbedo, vec3(0.0, 0.0, 0.5), step(-3.0, tileIndex.y) * step(tileIndex.y, -2.0));

    vec3 gapAlbedo = mix(vec3(0.0), vec3(0.2), Noise(uv * 25.0));
    vec3 albedo = mix(gapAlbedo, tileAlbedo, crackedTileMask);

    albedo *= mix(vec3(1.0, 1.0, 0.97), vec3(1.0), RandomValue(tileIndex));
    albedo *= mix(vec3(0.9), vec3(1.0), RandomValue(++tileIndex));

    float height = crackedTileMask * 0.2;
    float tilt = (UVRotate(tiledUV, vec2(0.5), RandomValue(tileIndex) * 360.0).x - 0.5) * 0.2;
    float wobble = Noise((uv + tileIndex) * 120.0) * 0.01;
    height = clamp(height + tilt - wobble, 0.0, 1.0);


    return Material(albedo, 0.0, roughness, 1.0, height);
}

Material Floor(vec2 uv)
{
    uv *= 0.5;

    vec2 tileIndex;
    vec2 tiledUV = UVTiling(uv, vec2(2.25), 0.0, tileIndex);
    float tileArea = step(-2.5, tileIndex.x) * step(tileIndex.x, 1.5);

    float tile = ShapeRectangle(tiledUV, vec2(0.99), vec2(0.02)) * tileArea;
    float wobble = Noise(tiledUV * 2.0 + tileIndex.yx) * 0.05;
    float tilt = (UVRotate(tiledUV, vec2(0.5), RandomValue(tileIndex) * 360.0).x - 0.5);

    vec2 gumIndex;
    vec2 gumUV = UVTiling(tiledUV, vec2(6.0), 0.0, gumIndex);
    gumUV = UVWarp(gumUV, Noise(uv * 100.0), 0.02);
    vec2 gumPosition = -vec2(0.5) + (RandomValue(gumIndex + tileIndex * 10.0) * 0.5) - 0.25;
    float hasGum = step(0.98, RandomValue(gumIndex + tileIndex * 10.0));
    float gum = ShapeCircle(gumUV + gumPosition, 0.2, 0.05) * hasGum;

    float height = (tile * 0.9) - wobble - (tilt * 0.2) + (gum * 0.1);
    height += 1.0 - tileArea;
    height = clamp(height, 0.0, 1.0);
    float waterMask = step(height, 0.82);
    float heightWW = clamp(height, 0.82, 1.0) * 0.1;

    vec2 roughnessUV = UVWarp(uv, Noise(uv * 3.0), 0.01);
    float roughness = Noise(roughnessUV.yx * 15.0) * 0.8;
    roughness += Noise(roughnessUV.xy * 40.0) * 0.2;
    roughness *= roughness;
    //roughness = roughness * 0.5 + 0.2;
    roughness *= (1.0 - waterMask * 0.9);
    roughness *= (1.0 - gum * 0.5);

    float albedoMix = Noise(tiledUV * 25.0 + tileIndex.yx) * 0.8
            + Noise(tiledUV.yx * 50.0 + tileIndex) * 0.2;
    albedoMix = (albedoMix * albedoMix);
    vec3 albedo = mix(vec3(0.8392, 0.8157, 0.7098), vec3(0.2549, 0.2549, 0.2549), albedoMix);
    vec3 gumColor = mix(vec3(0.0), vec3(0.4196, 0.4078, 0.302), RandomValue(gumIndex + tileIndex));
    albedo = mix(albedo, gumColor, gum);
    float dirtMix = Noise(UVWarp(uv * 7.0, Noise(uv * 1.0), 0.3));
    albedo = mix(albedo, vec3(0.1216, 0.1098, 0.1098), dirtMix);
    albedo *= 1.0 - (waterMask * 0.5);
    albedo *= tile;

    return Material(albedo, 0.0, roughness, 1.0, heightWW);
}

Material EvaluateMaterial(float id, vec2 uv)
{
    switch(int(id))
    {
        case 0:
            return Tiles(uv);
        case 1:
            return Material(vec3(0.0471, 0.2314, 0.0471), 0.0, 0.8, 1.0, 0.0);
        case 2:
            return Floor(uv);
        case 3:
            return Material(COLORALUMINIUM, 1.0, 0.1, 1.0, 0.0);
        case 4:
            return Material(vec3(0.5), 0.0, 1.0, 1.0, 0.0);
    }
}










//-------------------------------------------
// DF Tools
//-------------------------------------------
float DFSphere(vec3 point, vec3 origin, float radius)
{
    return length(point - origin) - radius;
}

float DFBox(vec3 point, vec3 origin, vec3 size, float roundness)
{
    vec3 d = abs(point - origin) - size;
    return length(max(d, 0.0)) - roundness
        + min(max(d.x,max(d.y,d.z)),0.0);
}

DFInfo DFMatSphere(vec3 point, vec3 origin, float radius, vec3 materials)
{
    vec3 diff = (point - origin);
    vec3 absDiff = abs(diff);
    vec3 dominantAxis = step(absDiff.yxx, absDiff.xyz) * step(absDiff.zzy, absDiff.xyz);

    float dist = length(diff) - radius
    ;

    materials *= dominantAxis;
    float materialId = materials.x + materials.y + materials.z
    ;

    vec2 uvX = (vec2(atan(diff.z, diff.x), atan(diff.x, diff.y)) + PI) / PI2;
    vec2 uvY = (vec2(atan(diff.y, diff.x), atan(diff.z, diff.y)) + PI) / PI2;
    vec2 uvZ = (vec2(atan(diff.z, diff.x), atan(diff.z, diff.y)) + PI) / PI2;
    vec2 uv = dominantAxis.x * uvX
            + dominantAxis.y * uvY
            + dominantAxis.z * uvZ
    ;

    return DFInfo(dist, origin, materialId, uv);
}

DFInfo DFMatBox(vec3 point, vec3 origin, vec3 size, float roundness, vec3 materials)
{
    vec3 diff = abs(point - origin) - size;
    vec3 dominantAxis = step(diff.yxx, diff.xyz) * step(diff.zzy, diff.xyz);

    float dist = length(max(diff, 0.0)) - roundness
            + min(max(diff.x, max(diff.y, diff.z)), 0.0)
    ;

    materials *= dominantAxis;
    float materialId = materials.x + materials.y + materials.z
    ;

    vec3 ref = point - origin;
    vec2 uv = dominantAxis.x * ref.zy
            + dominantAxis.y * ref.xz
            + dominantAxis.z * ref.xy
    ;

    return DFInfo(dist, origin, materialId, uv);
}

DFInfo MinDFInfo(DFInfo a, DFInfo b)
{
    float isA = step(a.hitDistance, b.hitDistance);
    a.hitDistance = isA * a.hitDistance + (1.0 - isA) * b.hitDistance;
    a.materialId = isA * a.materialId + (1.0 - isA) * b.materialId;
    a.uv = isA * a.uv + (1.0 - isA) * b.uv;

    return a;
}











//-------------------------------------------
// Distancefield tracing
//-------------------------------------------
DFInfo ComplexDF(vec3 point)
{
    vec3 roomSize = vec3(2.0, 1.5, 10.0);
    vec3 roomOrigin = vec3(0.0, 1.5, 0.0);
    vec3 roomMaterials = vec3(0.0, 2.0, 0.0);

    //defining stair
    float stepCount = 7.0;
    float stepDepth = 0.30;
    float stepHeight = 0.18;
    float stairStart = 0.0 - (stepCount * stepDepth * 0.5);

    // splitting space for repeat (global space)
    float spaceIndex;
    point.xz = UVSplitSpaceRadial(point.xz, 8.0, 15.0, spaceIndex);
    float spaceModifier = spaceIndex * stepHeight * stepCount;
    point.y -= spaceModifier;

    // manipulating local space (room space)
    float stairIndex = (1.0 / stepDepth) * (point.z - stairStart);
    float roomModifier = clamp(stairIndex * stepHeight, 0.0, stepCount * stepHeight);
    vec3 roomPoint = vec3(point.x, point.y - roomModifier, point.z);


    //// calc Room -> origin / materialId / distance
    DFInfo info = DFMatBox(roomPoint, roomOrigin, roomSize, 0.0, roomMaterials);
    float stairSegment = step(stairStart, point.z)
            * step(point.z, stairStart + stepDepth * stepCount);
    info.materialId = (1.0 - stairSegment) * info.materialId + stairSegment * 1.0;
    info.origin = roomOrigin + vec3(0.0, roomModifier + spaceModifier, 0.0);
    info.hitDistance *= -1.0;


    //// calc Stairs -> origin / materialId / distance
    vec3 stepPoint = point;
    vec3 next = vec3(0.0, stepHeight, stepDepth);
    vec3 stairPosition = vec3(0.0, 0.0, stairStart) + next * 0.5;
    for(float index = 0.0; index < stepCount; index++)
    {
        vec3 stepPosition = stairPosition + next * index;
        DFInfo stairStepInfo = DFMatBox(
            stepPoint,
            stepPosition,
            vec3(4.0, stepHeight, stepDepth) * 0.5,
            0.0,
            vec3(4.0)
        );
        info = MinDFInfo(info, stairStepInfo);
    }


    //// calc Railing -> origin / materialId / distance
    vec3 railPoint = vec3(abs(roomPoint.x), roomPoint.y, 0.0);
    DFInfo railingInfo = DFMatSphere(railPoint, vec3(1.85, 1.0, 0.0), 0.035, vec3(3.0));
    info = MinDFInfo(info, railingInfo);


    //// calc pipes -> origin / materialId / distance
    vec3 pipePoint = vec3(mod(clamp(abs(roomPoint.x - 1.0), 0.0, 0.6), 0.3), roomPoint.y, 0.0);
    DFInfo pipeInfo = DFMatSphere(pipePoint, vec3(0.15, 2.9, 0.0), 0.1, vec3(1.0));
    info = MinDFInfo(info, pipeInfo);



    //// Detailed surface height
    if (info.hitDistance < 0.01)
        info.hitDistance -= EvaluateMaterial(info.materialId, info.uv).height * 0.01;

    return info;
}

Trace TraceComplexDF(Ray ray)
{    
    vec3 point = ray.origin;
    float traceDistance = 0.0;
    float nextDistance = STEPTOLERANCE + 1.0;
    int steps = 0;
    DFInfo dfInfo;

    while(nextDistance > STEPTOLERANCE && steps < MAXSTEPS && traceDistance < MAXDISTANCE)
    {
        dfInfo = ComplexDF(point);
        nextDistance = dfInfo.hitDistance;
        point += ray.direction * nextDistance;
        traceDistance += nextDistance;
        steps += 1;
    }

    dfInfo.hitDistance = traceDistance;
    return Trace(dfInfo, point, steps);
}

vec3 AproxNormalComplexDF(vec3 point)
{
    vec2 dir = vec2(1,-1);
    return normalize(
            dir.xyy * ComplexDF(point + dir.xyy * NORMALSAMPLEDISTANCE).hitDistance + 
            dir.yyx * ComplexDF(point + dir.yyx * NORMALSAMPLEDISTANCE).hitDistance + 
            dir.yxy * ComplexDF(point + dir.yxy * NORMALSAMPLEDISTANCE).hitDistance + 
            dir.xxx * ComplexDF(point + dir.xxx * NORMALSAMPLEDISTANCE).hitDistance
    );
}

SurfaceHit EvaluateComplexSDF(Ray ray)
{
    Trace trace = TraceComplexDF(ray);
    vec3 normal = AproxNormalComplexDF(trace.point);

    return SurfaceHit(trace.distanceField.hitDistance,
            trace.point,
            normal,
            trace.steps,
            trace.distanceField.origin,
            trace.distanceField.materialId,
            trace.distanceField.uv
    );
}

float ShadowComplexDF(Ray ray, float maxDistance, float smoothness)
{
    vec3 point = ray.origin;
    float stepDistance = STEPTOLERANCE + 1.0;
    float traceDistance = 0.0;
    float shadow = 1.0;
    
    while(traceDistance < maxDistance)
    {
        if (stepDistance < STEPTOLERANCE)
            return 1.0;

        stepDistance = ComplexDF(point).hitDistance;
        point += ray.direction * stepDistance;
        traceDistance += stepDistance;
        shadow = min(shadow, smoothness * stepDistance / traceDistance);
    }

    return 1.0 - shadow;
}












//-------------------------------------------
// Direct Illumination
//-------------------------------------------
float Flicker(vec2 seed, float strength, float intervall)
{
    float noise = Noise(seed * intervall) * 2.0 - 0.5;
    noise *= 5.0;
    noise = clamp(noise, 0.0, 0.9);
    noise += 0.01;
    noise *= strength;

    return noise;
}

void ProvideLights(out Light ambient, out PointLight point[2])
{
    float time = iTime * .1;
    float lightStrength = Flicker(vec2(time, 0.0), 20.0, 40.0);
    
    point[0] = PointLight(
            Light(vec3(1.0, 1.0, 1.0), Flicker(vec2(time, 0.0), 20.0, 40.0)),
            vec3(-15.0 + 5.0, 4.5, 11.0),
            1.0,
            1.0,
            32.0
    );
    point[1] = PointLight(
            Light(vec3(1.0, 1.0, 1.0), 20.0),
            vec3(-15.0, 3.7, 2.0),
            1.0,
            1.0,
            32.0
    );

    ambient = Light(vec3One, 0.02);
}

vec3 FresnelSchlick(float nDot, vec3 reflectionColor, float roughness)
{
    return reflectionColor
            + (max(vec3(1.0 - roughness), reflectionColor) - reflectionColor)
            * pow(1.0 - nDot, 5.0);
}

float DistributionTrowbridgeReitzGGX(vec3 normal, vec3 halfway, float roughness)
{
    roughness = pow(roughness, 4.0);
    
    float nDotH = max(0.0, dot(normal, halfway));    
    float denom = (nDotH * nDotH * (roughness - 1.0) + 1.0);
	
    return roughness / (PI * denom * denom);
}

float GeometrySchlickGGX(float nDotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;
    float denom = nDotV * (1.0 - k) + k;
	
    return nDotV / denom;
}

float GeometrySmith(vec3 normal, vec3 viewDirection, vec3 lightDirection, float roughness)
{
    float nDotV = max(dot(normal, viewDirection), 0.0);
    float nDotL = max(dot(normal, lightDirection), 0.0);
    float ggxView  = GeometrySchlickGGX(nDotV, roughness);
    float ggxLight  = GeometrySchlickGGX(nDotL, roughness);
	
    return ggxView * ggxLight;
}

SurfaceColor DirectPBRBRDF(Ray ray, SurfaceHit hit, Material material)
{
    Light ambientLight;
    PointLight[2] PointLights;
    ProvideLights(ambientLight, PointLights);


    // Ambient light
    vec3 evaluatedColor = ambientLight.color
            * ambientLight.intensity
            * material.albedo
            * material.occlusion;
            
    // Constants per pointlight
    vec3 viewDirection = -ray.direction;
    vec3 fresnelReflectionColor = mix(vec3(0.04), material.albedo, material.metallic);
    float nDotView = max(0.0, dot(hit.normal, -ray.direction));

    // Direct lighting per light
    for(int index = 0; index < 2; index++)
    {
        // PREPERATIONS
        PointLight pointLight = PointLights[index];
        vec3 lightVector = pointLight.position - hit.point;
        vec3 lightDirection = normalize(lightVector);
        vec3 halfway = normalize(viewDirection + lightDirection);

        // SPECULAR COLOR
        // Fresnel
        float hDotView = max(0.0, dot(hit.normal, viewDirection));
        vec3 fresnelColor = FresnelSchlick(
                hDotView,
                fresnelReflectionColor,
                material.roughness
        );

        // Microsurface normal distribution
        float normalDistribution = DistributionTrowbridgeReitzGGX(
                hit.normal,
                halfway,
                material.roughness
        );

        // Microsurface selfshadowing
        float selfShadowing = GeometrySmith(
                hit.normal,
                viewDirection,
                lightDirection,
                material.roughness
        );

        // Specular aggregation
        float nDotLight = max(0.0, dot(hit.normal, lightDirection));
        vec3 specularColor = fresnelColor * normalDistribution * selfShadowing;
        specularColor /= max(0.0001, 4.0 * nDotView * nDotLight);

        // DIFFUSE COLOR
        vec3 diffuseColor = material.albedo
                * (vec3(1.0) - fresnelColor)
                * (1.0 - material.metallic);

        // RADIANCE COLOR
        float lightDistance = length(lightVector);
        float distanceAttenuation = mix(1.0, 1.0 / pow(1.0 + lightDistance, 2.0), pointLight.falloff);
        vec3 radiance = pointLight.light.color
                * pointLight.light.intensity
                * distanceAttenuation;

        // SHADOW
        vec3 safetyOffset = hit.normal * 0.001;
        Ray shadowRay = Ray(hit.point + safetyOffset, lightDirection);
        float hasShadow = ShadowComplexDF(shadowRay, length(lightVector), pointLight.shadowSmoothness);
        hasShadow *= pointLight.shadowStrength;

        // FINAL MIX
        evaluatedColor += ((diffuseColor / PI) + specularColor)
                * radiance
                * nDotLight
                * (1.0 - hasShadow)
                * material.occlusion
        ;
    }

    return SurfaceColor(evaluatedColor);
}










//-------------------------------------------
// General programm / coordination
//-------------------------------------------

Ray CameraRay(vec2 pixel, vec3 position, vec3 direction, float fov)
{
    // ray spread
    float fovRad = fov * RADPERDEG;
    float stretch = tan(fovRad / 2.0) / iResolution.y;
    vec2 camPlane = (pixel * 2.0 - iResolution.xy) * stretch;
    vec3 rawDirection = normalize(vec3(camPlane.x, camPlane.y, 1.0));

    // rotate towards view direction
    vec3 f = -direction;
    vec3 s = normalize(cross(f, vec3Up));
    vec3 u = cross(s, f);
    mat3 rotationMatrix = mat3(s, u, -f);

    return Ray(position, rotationMatrix * rawDirection);
}

vec3 GammaCorrection(vec3 linearColor)
{
    linearColor = linearColor / (linearColor + vec3(1.0));
    return pow(linearColor, vec3(0.45454545)); // (1.0 / 2.2) --> 0.45454545...
}

void main() {
    vec3 pixelColor = BACKDROPCOLOR;

    //// Camera setup and ray creation
    vec3 cameraPosition = vec3(-16.5, 1.7, -5.0 + iTime * 0.1);
    vec3 viewTarget = cameraPosition + vec3(0.3, 0.03, 1.0);
    vec3 viewDirection = normalize(viewTarget - cameraPosition);
    Ray ray = CameraRay(gl_FragCoord.xy, cameraPosition, viewDirection, 60.0);


    //// Check if anything got hit
    SurfaceHit hit = EvaluateComplexSDF(ray);
    if (hit.steps < MAXSTEPS || hit.hitDistance < MAXDISTANCE)
    {
        //// Get surface composition
        Material material = EvaluateMaterial(hit.materialId, hit.basicUV);

        //// Calculate direct lighting
        vec3 directColor = DirectPBRBRDF(ray, hit, material).color;

        // gamma corection
        pixelColor = GammaCorrection(directColor);
    }

    gl_FragColor = vec4(pixelColor, 1.0);
}