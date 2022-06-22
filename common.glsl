/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}


vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    //Calculate eye_offset and ray direction
    vec3 eye_offset = cam.eye + cam.u * ls.x + cam.v * ls.y;

    float x = cam.focusDist * cam.width * (pixel_sample.x  / iResolution.x - 0.5f);
    float y = cam.focusDist * cam.height * (pixel_sample.y  / iResolution.y - 0.5f);
    float z = -cam.focusDist * cam.planeDist;

    vec3 p = vec3(x, y, z) - vec3(ls, 0.0f);

    vec3 ray_dir = p.x * cam.u + p.y * cam.v + p.z * cam.n;
    ray_dir = normalize(ray_dir);
    
    return createRay(eye_offset, normalize(ray_dir), time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for dialectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDialectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIALECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};


float schlick(float cosine, float refIdx)
{
    float R0 = pow((refIdx - 1.0f) / (refIdx + 1.0f), 2.0f);
    float kr = R0 + (1.0f - R0) * pow(1.0f - cosine, 5.0f);

    return kr;
}

vec3 metalSchlick(float cosine, vec3 F0) {
    return F0 + (1.0f - F0) * pow(clamp(1.0f - cosine, .0f, 1.0f), 5.0f); // PROF CODE
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if(rec.material.type == MT_DIFFUSE)
    {
        vec3 s = rec.pos + rec.normal + normalize(randomInUnitSphere(gSeed));
        rScattered = createRay(rec.pos + epsilon * rec.normal, normalize(s - rec.pos), rIn.t);
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
        vec3 refletedRayDirection = reflect(rIn.d, rec.normal);
        vec3 fuzzyDirection = rec.pos + refletedRayDirection + rec.material.roughness * randomInUnitSphere(gSeed);
        rScattered = createRay(rec.pos + epsilon * rec.normal, normalize(fuzzyDirection - rec.pos), rIn.t);
        atten = metalSchlick(-dot(rIn.d, rec.normal), rec.material.specColor);
        return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {

        atten = vec3(1.0f);
        vec3 outwardNormal;
        float niOverNt, cosine;

        if (dot(rIn.d, rec.normal) > .0f)  // From inside
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            atten = exp(-rec.material.refractColor * rec.t);  // Beer's Law
        }
        else  // From outside
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0f / rec.material.refIdx;
        }

        float reflectProb, cosT;
        vec3 v = rIn.d * -1.0;
        vec3 vt = (outwardNormal * dot(v, outwardNormal)) - v;
        float sinI = length(vt);
        float cosI = sqrt(1.0f - pow(sinI, 2.0f));
        float sinT = niOverNt * sinI;

        if (sinT < 1.0f)  // non-total reflection
        {
            cosT = sqrt(1.0f - pow(sinT, 2.0f));
            cosine = (dot(rIn.d, outwardNormal) > 0.0f) ? cosT : cosI;
            reflectProb = schlick(cosine, rec.material.refIdx);
        }
        else  // total reflection
        {
            reflectProb = 1.0f;
        }

        if(hash1(gSeed) < reflectProb)  //Reflection
        {
            vec3 refletedRayDirection = reflect(rIn.d, rec.normal);
            vec3 fuzzyDirection = rec.pos + refletedRayDirection + rec.material.roughness * randomInUnitSphere(gSeed);
            rScattered = createRay(rec.pos + rec.normal * epsilon, normalize(fuzzyDirection - rec.pos), rIn.t);
            // atten *= vec3(reflectProb); not necessary since we are only scattering reflectProb rays and not all reflected rays
        }
        else  //Refraction
        {
            vec3 refractedRayDirection = normalize(vt) * sinT - outwardNormal * cosT;
            vec3 fuzzyDirection = rec.pos + refractedRayDirection + rec.material.roughness * randomInUnitSphere(gSeed);
            rScattered = createRay(rec.pos - outwardNormal * epsilon, normalize(fuzzyDirection - rec.pos), rIn.t);
            // atten *= vec3(1.0 - reflectProb); not necessary since we are only scattering 1-reflectProb rays and not all refracted rays
        }

        return true;
    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle tri, Ray r, float tmin, float tmax, out HitRecord rec)
{
    // Tomas Möller, Ben Trumbore “Fast, Minimum Storage Ray/Triangle Intersection”, 1997
    vec3 v1 = tri.b - tri.a;
	vec3 v2 = tri.c - tri.a;
	vec3 v3 = r.d * -1.0;
	vec3 vs = r.o - tri.a;

	float det = dot(v1, cross(v2, v3));
	if (det == 0.0f) return false;

	float beta = dot(vs, cross(v2, v3)) / det;
	float gamma = dot(v1, cross(vs, v3)) / det;
	float t = dot(v1, cross(v2, vs)) / det;

	if (beta < 0.0f || beta > 1.0f || gamma < 0.0f || gamma > 1.0f) return false;
	if ((beta + gamma) > 1.0f) return false;

    vec3 normal = cross((tri.b - tri.a), (tri.c - tri.a));
	normal = normalize(normal);

    //calculate a valid t and normal
    if(t < tmax && t > tmin)
    {
        rec.t = t;
        rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}

struct Quad {vec3 a; vec3 b; vec3 c; vec3 d;};

Quad createQuad(vec3 v0, vec3 v1, vec3 v2, vec3 v3)
{
    Quad q;
    q.a = v0; q.b = v1; q.c = v2; q.d = v3;
    return q;
}


bool hit_quad(Quad quad, Ray r, float tmin, float tmax, out HitRecord rec)
{
    // Retrieved from https://www.shadertoy.com/view/tsBBWW
    // Calculates normal and flip vertices order, if needed
    vec3 normal = normalize(cross(quad.b - quad.a, quad.c - quad.a));
    if (dot(normal, r.d) > 0.0f)
    {
        normal *= -1.0f;
		vec3 temp = quad.d;
        quad.d = quad.a;
        quad.a = temp;
        temp = quad.b;
        quad.b = quad.c;
        quad.c = temp;
    }
    
    // -----------------------------------------------------------------

    // Adapted for Quads from Tomas Möller, Ben Trumbore
    // “Fast, Minimum Storage Ray/Triangle Intersection”, 1997
    vec3 v1 = quad.b - quad.a;
	vec3 v2 = quad.d - quad.a;
	vec3 v3 = r.d * -1.0;
	vec3 vs = r.o - quad.a;

	float det = dot(v1, cross(v2, v3));
	if (det == 0.0f) return false;

	float beta = dot(vs, cross(v2, v3)) / det;
	float gamma = dot(v1, cross(vs, v3)) / det;
	float t = dot(v1, cross(v2, vs)) / det;

	if (beta < 0.0f || beta > 1.0f || gamma < 0.0f || gamma > 1.0f) return false;

    //calculate a valid t and normal
    if(t < tmax && t > tmin)
    {
        rec.t = t;
        rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}

struct Box {
    vec3 minPoint, maxPoint;
};

Box createBox(vec3 minPoint, vec3 maxPoint){
    Box b;
    b.minPoint = minPoint;
    b.maxPoint = maxPoint;
    return b;
}

// Kay & Kajiya algorithm
bool hit_box(Box box, Ray r, float tmin, float tmax, out HitRecord rec) {

	float ox = r.o.x;
	float oy = r.o.y;
	float oz = r.o.z;

	float dx = r.d.x;
	float dy = r.d.y;
	float dz = r.d.z;

	float tx_min, ty_min, tz_min, tx_max, ty_max, tz_max;

	// CALCULATE INTERSECTIONS
	float a = 1.0f / dx;
	if (a >= 0.0f) {
		tx_min = (box.minPoint.x - ox) * a;
		tx_max = (box.maxPoint.x - ox) * a;
	}
	else { 
		tx_min = (box.maxPoint.x - ox) * a;
		tx_max = (box.minPoint.x - ox) * a;
	}

	float b = 1.0f / dy;
	if (b >= 0.0f) {
		ty_min = (box.minPoint.y - oy) * b;
		ty_max = (box.maxPoint.y - oy) * b;
	}
	else {
		ty_min = (box.maxPoint.y - oy) * b;
		ty_max = (box.minPoint.y - oy) * b;
	}

	float c = 1.0f / dz;
	if (c >= 0.0f) {
		tz_min = (box.minPoint.z - oz) * c;
		tz_max = (box.maxPoint.z - oz) * c;
	}
	else {
		tz_min = (box.maxPoint.z - oz) * c;
		tz_max = (box.minPoint.z - oz) * c;
	}

	float tE, tL; // Entering and leaving t values
	vec3 face_in, face_out; // Normals

	// FIND LARGEST ENTERING T VALUE
	if (tx_min > ty_min) {
		tE = tx_min;
		face_in = (a >= 0.0f) ? vec3(-1.0f, 0.0f, 0.0f) : vec3(1.0f, 0.0f, 0.0f);
	}
	else {
		tE = ty_min;
		face_in = (b >= 0.0f) ? vec3(0.0f, -1.0f, 0.0f) : vec3(0.0f, 1.0f, 0.0f);
	}
	if (tz_min > tE) {
		tE = tz_min;
		face_in = (c >= 0.0f) ? vec3(0.0f, 0.0f, -1.0f) : vec3(0.0f, 0.0f, 1.0f);
	}

	// FIND SMALLEST LEAVING T VALUE
	if (tx_max < ty_max) {
		tL = tx_max;
		face_out = (a >= 0.0f) ? vec3(1.0f, 0.0f, 0.0f) : vec3(-1.0f, 0.0f, 0.0f);
	}
	else {
		tL = ty_max;
		face_out = (b >= 0.0f) ? vec3(0.0f, 1.0f, 0.0f) : vec3(0.0f, -1.0f, 0.0f);
	}
	if (tz_max < tL) {
		tL = tz_max;
		face_out = (c >= 0.0f) ? vec3(0.0f, 0.0f, 1.0f) : vec3(0.0f, 0.0f, -1.0f);
	}


    float t;
    vec3 normal;

	// Condition to hit
	if (tE < tL && tL > 0.0f) {
		if (tE > 0.0f) {
		// ray hits outside surface
			t = tE;
			normal = face_in;
		}
		else {
		// ray hits inside surface
			t = tL;
			normal = face_out;
		}

        //calculate a valid t and normal
        if(t < tmax && t > tmin)
        {
            rec.t = t;
            rec.normal = normal;
            rec.pos = pointOnRay(r, rec.t);
            return true;
        }
	}
	return false;
}

struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    vec3 moving_center = mvsphere.center0 + (mvsphere.center1 - mvsphere.center0) * time;
    return moving_center;
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{

    //INSERT YOUR CODE HERE
    vec3 rayDirection = normalize(r.d);
    vec3 oc = s.center - r.o;

    float b = dot(oc, rayDirection);
    float c = dot(oc, oc) - pow(s.radius, 2.0f);
    float discriminant = pow(b, 2.0f) - c;
    float t;

	if (c > .0f) {
		if (b <= .0f) return false;
		if (discriminant <= .0f) return false;
		t = b - sqrt(discriminant);
	}
	else {
		t = b + sqrt(discriminant);
	}
    
    //calculate a valid t and normal
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normalize(rec.pos - s.center);
        return true;
    }
    return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    vec3 movingCenter = center(s, r.t);
    vec3 rayDirection = normalize(r.d);
    vec3 oc = movingCenter - r.o;

    float b = dot(oc, rayDirection);
    float c = dot(oc, oc) - pow(s.radius, 2.0f);
    float discriminant = pow(b, 2.0f) - c;
    float t;

	if (c > .0f) {
		if (b <= .0f) return false;
		if (discriminant <= .0f) return false;
		t = b - sqrt(discriminant);
	}
	else {
		t = b + sqrt(discriminant);
	}
    
    //calculate a valid t and normal
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normalize(rec.pos - movingCenter);
        return true;
    }
    return false;
}

struct AAC {
    vec3 base, apex, normal;
    float baseR, apexR;
};

AAC createAAC(vec3 base, vec3 apex, float baseR, float apexR) {
    AAC c;
    c.base = base;
    c.apex = apex;
    c.baseR = baseR;
    c.apexR = apexR;

    return c;
}

bool hit_AAC(AAC c, Ray r, float tmin, float tmax, out HitRecord rec) {
    vec3 base = c.base;
    vec3 apex = c.apex;
    float baseR = c.baseR;
    float apexR = c.apexR;

    float t;

    bool x = base.x != apex.x;
    bool y = base.y != apex.y;
    bool z = base.z != apex.z;

    if((x || y) && (x || z) && (y || z)) return false;

    vec3 d = normalize(r.d);
    vec3 o = r.o;

    if (c.baseR == c.apexR) {
        float a = x ? pow(d.z, 2.0f) + pow(d.y, 2.0f) : y ? pow(d.x, 2.0f) + pow(d.z, 2.0f) : pow(d.x, 2.0f) + pow(d.y, 2.0f);
		float b = x ? d.z * (o.z - base.z) + d.y * (o.y - base.y) : y ? d.x * (o.x - base.x) + d.z * (o.z - base.z) : d.x * (o.x - base.x) + d.y * (o.y - base.y);
		float c = x ? pow(o.z - base.z, 2.0f) + pow(o.y - base.y, 2.0f) - pow(baseR, 2.0f) : y ? pow(o.x - base.x, 2.0f) + pow(o.z - base.z, 2.0f) - pow(baseR, 2.0f) : pow(o.x - base.x, 2.0f) + pow(o.y - base.y, 2.0f) - pow(baseR, 2.0f);

		float discriminant = pow(b, 2.0f) - a * c;

        if (discriminant <= .0f) return false;  // No intersection

		float sol1 = (-1.0f * b - sqrt(discriminant)) / a;
		float sol2 = (-1.0f * b + sqrt(discriminant)) / a;

		if (sol2 < .0f) return false;  // Cylinder behind ray origin

		// ----------------------------------------------------------------------------------------------------------------------------

		// Minimum and maximum coordenates along axis (coordenates of base/apex planes)
		float coordMin = x ? min(base.x, apex.x) : y ? min(base.y, apex.y) : min(base.z, apex.z);
		float coordMax = x ? max(base.x, apex.x) : y ? max(base.y, apex.y) : max(base.z, apex.z);

		// Coordenates of both intersection points along axis
		float coord1 = x ? (o + d * sol1).x : y ? (o + d * sol1).y : (o + d * sol1).z;
		float coord2 = x ? (o + d * sol2).x : y ? (o + d * sol2).y : (o + d * sol2).z;

		// Swap base and apex in order to ensure that base.x < apex.x, base.y < apex.y or base.z < apex.z
		if (coordMin != (x ? base.x : y ? base.y : base.z)) {
			vec3 temp = base;
			base = apex;
			apex = temp;
		}

		bool baseIntersection = false, apexIntersection = false, sideIntersection = false;

		if (sol1 < .0f)  // Ray origin inside infinite cylinder
		{
			float originCoord = x ? o.x : y ? o.y : o.z;
			bool isInside = (originCoord > coordMin) && (originCoord < coordMax);

			if (isInside && (coord2 < coordMin)) baseIntersection = true;
			if (!isInside && (coord2 > coordMin) && (length(o - base) < length(o - apex))) baseIntersection = true;

			if (isInside && (coord2 > coordMax)) apexIntersection = true;
			if (!isInside && (coord2 < coordMax) && (length(o - base) > length(o - apex))) apexIntersection = true;

			if (isInside && (coord2 > coordMin) && (coord2 < coordMax)) {
				sideIntersection = true;
				t = sol2;
			}
		}
		else if (sol1 > .0f)  // Cylinder in front of ray origin
		{
			// Intersection above and below base/apex planes
			if (coord1 < coordMin && coord2 < coordMin) return false;
			if (coord1 > coordMax && coord2 > coordMax) return false;

			t = sol1;

			if (coord1 > coordMin && coord1 < coordMax) sideIntersection = true;
			if (coord1 < coordMin && coord2 > coordMin) baseIntersection = true;
			if (coord1 > coordMax && coord2 < coordMax) apexIntersection = true;
		}

		// Intersection with base plane
		if (baseIntersection && t < tmax && t > tmin) {
			rec.normal = normalize(base - apex);
			rec.t = dot((base - o), rec.normal) / dot(rec.normal, d);
            rec.pos = pointOnRay(r, rec.t);
			return true;
		}

		// Intersection with apex plane
		if (apexIntersection && t < tmax && t > tmin) {
			rec.normal = normalize(apex - base);
			rec.t = dot((apex - o), rec.normal) / dot(rec.normal, d);
            rec.pos = pointOnRay(r, rec.t);
			return true;
		}

		// Intersection with side
		if (sideIntersection && t < tmax && t > tmin) {
			vec3 axis = apex - base;
			vec3 closestOnAxis = base + axis * ((((o + d * t) * axis) - dot(base, axis)) / dot(axis, axis));
			rec.normal = normalize((o + d * t) - closestOnAxis);
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);
			return true;
		}
    }
    else if (apexR == .0f)  // Cone
	{
		// Minimum and maximum coordenates along axis (coordenates of base/apex planes)
		float coordMin = x ? min(base.x, apex.x) : y ? min(base.y, apex.y) : min(base.z, apex.z);
		float coordMax = x ? max(base.x, apex.x) : y ? max(base.y, apex.y) : max(base.z, apex.z);

		float slope = baseR / (coordMax - coordMin);

		float a = x ? pow(d.y, 2.0f) + pow(d.z, 2.0f) - slope * pow(d.x, 2.0f) : y ? pow(d.x, 2.0f) + pow(d.z, 2.0f) - slope * pow(d.y, 2.0f) : pow(d.x, 2.0f) + pow(d.y, 2.0f) - slope * pow(d.z, 2.0f);
		float b = x ? d.y * (o.y - apex.y) + d.z * (o.z - apex.z) - (slope * d.x * (o.x - apex.x)) : y ? d.x * (o.x - apex.x) + d.z * (o.z - apex.z) - (slope * d.y * (o.y - apex.y)) : d.x * (o.x - apex.x) + d.y * (o.y - apex.y) - (slope * d.z * (o.z - apex.z));
		float c = x ? pow((o.y - apex.y), 2.0f) + pow((o.z - apex.z), 2.0f) - slope * pow((o.x - apex.x), 2.0f) : y ? pow((o.x - apex.x), 2.0f) + pow((o.z - apex.z), 2.0f) - slope * pow((o.y - apex.y), 2.0f) : pow((o.x - apex.x), 2.0f) + pow((o.y - apex.y), 2.0f) - slope * pow((o.z - apex.z), 2.0f);

		float discriminant = pow(b, 2.0f) - a * c;

		if (discriminant <= .0f) return false;  // No intersection

		float sol1 = (-1.0f * b - sqrt(discriminant)) / a;
		float sol2 = (-1.0f * b + sqrt(discriminant)) / a;

		// Coordenates of both intersection points along axis
		float coord1 = x ? (o + d * sol1).x : y ? (o + d * sol1).y : (o + d * sol1).z;
		float coord2 = x ? (o + d * sol2).x : y ? (o + d * sol2).y : (o + d * sol2).z;

		vec3 axis = (apex - base);
		float height = length(axis);
		axis = normalize(axis);

		// distance between base and projection of intersections in axis
		float c1 = dot(((o + d * sol1) - base), axis);
		float c2 = dot(((o + d * sol2) - base), axis);

		float originCoord = x ? o.x : y ? o.y : o.z;
		bool isInside = (originCoord > coordMin) && (originCoord < coordMax);

		if (sol2 < .0f && (c1 < height || c2 > height || c2 < .0f)) return false;

		bool baseIntersection = false, sideIntersection = false;

		if (sol2 < .0f) {
			if (isInside && c1 > height && c2 > .0f && c2 < height) baseIntersection = true;  // DENTRO
		}
		if (sol2 > .0f && sol1 < .0f) {
			if (c1 > height && c2 > .0f && c2 < height) sideIntersection = true;  // t = sol2;  FORA
			if (isInside && c2 > .0f && c2 < height) sideIntersection = true;  // t = sol2;  DENTRO
			if (isInside && c1 > .0f && c1 < height && c2 < .0f) baseIntersection = true;  // DENTRO
			if (!isInside && c2 > .0f && c2 < height && c1 < .0f) baseIntersection = true;  // FORA
		}
		if (sol1 > .0f) {

			if (isInside && c2 > .0f && c2 < height) sideIntersection = true;  // t = sol2;  DENTRO
			if (c1 > height && c2 > .0f && c2 < height) sideIntersection = true; // t = sol2;
			if (c1 > .0f && c1 < height) sideIntersection = true; // t = sol1;

			if (c1 < .0f && c2 > .0f && c2 < height) baseIntersection = true;  // FORA
			if (!isInside && c1 > .0f && c1 < height && c2 > height) baseIntersection = true;  // FORA (DENTRO DO CIL INF)
		}

		if (sideIntersection) t = (sol1 > .0f && c1 > .0f && c1 < height) ? sol1 : sol2;

		if (baseIntersection && t < tmax && t > tmin) {
			rec.normal = normalize(base - apex);
			rec.t = dot((base - o), rec.normal) / dot(rec.normal, d);
            rec.pos = pointOnRay(r, rec.t);

			return true;
		}

		// Intersection with side surface
		if (sideIntersection && t < tmax && t > tmin) {
			axis = axis * height;
			vec3 closestOnAxis = base + axis * ((((o + d * t) * axis) - dot(base, axis)) / dot(axis, axis));
			vec3 cylinderNormal = normalize((o + d * t) - closestOnAxis);
			vec3 toApex = normalize(apex - (o + d * t));
			rec.normal = normalize(cross(cross(toApex, cylinderNormal), toApex));
            rec.t = t;
            rec.pos = pointOnRay(r, rec.t);

			return true;
		}
	}
	return false;

}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}