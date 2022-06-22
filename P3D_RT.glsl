/**
* ver hash functions em
* https://www.shadertoy.com/view/XlGcRh hash functions GPU
* http://www.jcgt.org/published/0009/03/02/
 */

 #include "./common.glsl"
 #iChannel0 "self"
 #iKeyboard


// ------------------------------------------------

// Scenes (only one should be selected)
bool SHIRLEY = false;
bool CORNELL_BOX_BOXES = false;
bool CORNELL_BOX_SPHERES = true;
bool CAUSTICS = false;

// Extra: Soft Shadows (only one should be selected)
bool SOFT_SHADOWS = true;
bool SOFT_SHADOWS_GRID  = false;

// ------------------------------------------------

bool hit_world(Ray r, float tmin, float tmax, out HitRecord rec)
{
    bool hit = false;
    rec.t = tmax;

    if (SHIRLEY) {

        if(hit_triangle(createTriangle(vec3(-10.0, -0.05, 10.0), vec3(10.0, -0.05, 10.0), vec3(-10.0, -0.05, -10.0)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.9));
        }

        if(hit_triangle(createTriangle(vec3(-10.0, -0.05, -10.0), vec3(10.0, -0.05, 10), vec3(10.0, -0.05, -10.0)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.9));
        }

        if(hit_sphere(
            createSphere(vec3(-4.0, 1.0, 0.0), 1.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2, 0.95, 0.1));
            //rec.material = createDiffuseMaterial(vec3(0.4, 0.2, 0.1));
        }

        if(hit_sphere(
            createSphere(vec3(4.0, 1.0, 0.0), 1.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createMetalMaterial(vec3(0.7, 0.6, 0.5), 0.0);
        }

        if(hit_sphere(
            createSphere(vec3(-1.5, 1.0, 0.0), 1.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDialectricMaterial(vec3(0.0), 1.3, 0.0);
        }

        if(hit_sphere(
            createSphere(vec3(-1.5, 1.0, 0.0), -0.55),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDialectricMaterial(vec3(0.0), 1.3, 0.0);
        }

        if(hit_sphere(
            createSphere(vec3(1.5, 1.0, 0.0), 1.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDialectricMaterial(vec3(0.0, 0.7, 0.9), 1.05, 0.05);
        }
    

        int numxy = 5;
        
        for(int x = -numxy; x < numxy; ++x)
        {
            for(int y = -numxy; y < numxy; ++y)
            {
                float fx = float(x);
                float fy = float(y);
                float seed = fx + fy / 1000.0;
                vec3 rand1 = hash3(seed);
                vec3 center = vec3(fx + 0.9 * rand1.x, 0.2, fy + 0.9 * rand1.y);
                float chooseMaterial = rand1.z;
                if(distance(center, vec3(4.0, 0.2, 0.0)) > 0.9)
                {
                    if(chooseMaterial < 0.3)
                    {
                        vec3 center1 = center + vec3(0.0, hash1(gSeed) * 0.5, 0.0);
                        // diffuse
                        if(hit_movingSphere(
                            createMovingSphere(center, center1, 0.2, 0.0, 1.0),
                            r,
                            tmin,
                            rec.t,
                            rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.5)
                    {
                        // diffuse
                        if(hit_sphere(
                            createSphere(center, 0.2),
                            r,
                            tmin,
                            rec.t,
                            rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.7)
                    {
                        // metal
                        if(hit_sphere(
                            createSphere(center, 0.2),
                            r,
                            tmin,
                            rec.t,
                            rec))
                        {
                            hit = true;
                        // rec.material.type = MT_METAL;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, 0.0);
                        }
                    }
                    else if(chooseMaterial < 0.9)
                    {
                        // metal
                        if(hit_sphere(
                            createSphere(center, 0.2),
                            r,
                            tmin,
                            rec.t,
                            rec))
                        {
                            hit = true;
                        // rec.material.type = MT_METAL;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, hash1(seed));
                        }
                    }
                    else
                    {
                        // glass (dialectric)
                        if(hit_sphere(
                            createSphere(center, 0.2),
                            r,
                            tmin,
                            rec.t,
                            rec))
                        {
                            hit = true;
                            rec.material.type = MT_DIALECTRIC;
                            rec.material = createDialectricMaterial(hash3(seed), 1.5, 0.0);
                        }
                    }
                }
            }
        }
    }

    float yShift = -3.5f;
    float zShift = -8.0f;
    float margin = 0.05f;
    float boxSize = 6.0f;
    float boxHeight = 8.0f;

    if (CORNELL_BOX_BOXES || CORNELL_BOX_SPHERES) {

        yShift = -3.5f;
        zShift = -8.0f;
        margin = 0.05f;
        boxSize = 6.0f;
        boxHeight = 8.0f;

        // Floor
        if(hit_quad(createQuad(vec3((boxSize+margin), -margin+yShift, (boxSize+margin)+zShift), vec3((boxSize+margin), -margin+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), -margin+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), -margin+yShift, (boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(1.0));
        }

        // Ceiling
        if(hit_quad(createQuad(vec3((boxSize+margin), boxHeight+yShift, (boxSize+margin)+zShift), vec3((boxSize+margin), boxHeight+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), boxHeight+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), boxHeight+yShift, (boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(1.0));
        }
        
        // Left Wall (Red)
        if(hit_quad(createQuad(vec3(-boxSize, -margin+yShift, (boxSize+margin)+zShift), vec3(-boxSize, -margin+yShift, -(boxSize+margin)+zShift), vec3(-boxSize, (boxHeight+margin)+yShift, -(boxSize+margin)+zShift), vec3(-boxSize, (boxHeight+margin)+yShift, (boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.95, margin, margin));
        }
        
        // Right Wall (Green)
        if(hit_quad(createQuad(vec3(boxSize, -margin+yShift, -(boxSize+margin)+zShift), vec3(boxSize, -margin+yShift, (boxSize+margin)+zShift), vec3(boxSize, boxHeight+yShift, (boxSize+margin)+zShift), vec3(boxSize, boxHeight+yShift, -(boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(margin, 0.95, margin));
        }

        // Front Wall
        if(hit_quad(createQuad(vec3(-(boxSize+margin), -margin+yShift, -(boxSize+margin)+zShift), vec3((boxSize+margin), -margin+yShift, -boxSize+zShift), vec3((boxSize+margin), boxHeight+yShift, -boxSize+zShift), vec3(-(boxSize+margin), boxHeight+yShift, -boxSize+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(1.0));
        }

        // Sphere
        if(hit_sphere(
            createSphere(vec3(0.0, 1.0, -3.0), 1.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDialectricMaterial(vec3(0.0, 1.0, 0.0), 1.10, 0.0);
        }
    }

    if (CAUSTICS) {

        yShift = -1.5f;
        zShift = -3.5f;
        margin = 0.05f;
        boxSize = 6.0f;

        // Floor
        if(hit_quad(createQuad(vec3((boxSize+margin), -margin+yShift, (boxSize+margin)+zShift), vec3((boxSize+margin), -margin+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), -margin+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), -margin+yShift, (boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.95));
        }

        // Ceiling
        /*
        if(hit_quad(createQuad(vec3((boxSize+margin), boxSize+yShift, (boxSize+margin)+zShift), vec3((boxSize+margin), boxSize+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), boxSize+yShift, -(boxSize+margin)+zShift), vec3(-(boxSize+margin), boxSize+yShift, (boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.95));
        }*/

        // Piece of ceiling containing the light
        float s = 5.0f; // Diminuicao em relação ao plano de cima 
        float z = -2.5f; // Aumento no z
        float y = 0.6f; // Aumento no y
        if(hit_quad(createQuad(vec3((boxSize+margin)/s, boxSize+yShift+y, ((boxSize+margin)+zShift)/s+z), vec3((boxSize+margin)/s, boxSize+yShift+y, (-(boxSize+margin)+zShift)/s+z), vec3(-(boxSize+margin)/s, boxSize+yShift+y, (-(boxSize+margin)+zShift)/s+z), vec3(-(boxSize+margin)/s, boxSize+yShift+y, ((boxSize+margin)+zShift)/s+z)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.95));
            rec.material.emissive = vec3(1.0) * 10.0f;
        }
        
        // Left Wall (Red)
        
        if(hit_quad(createQuad(vec3(-boxSize, -margin+yShift, (boxSize+margin)+zShift), vec3(-boxSize, -margin+yShift, -(boxSize+margin)+zShift), vec3(-boxSize, (boxSize+margin)+yShift, -(boxSize+margin)+zShift), vec3(-boxSize, (boxSize+margin)+yShift, (boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.95, margin, margin));
        }
        
        // Right Wall (Green)
        if(hit_quad(createQuad(vec3(boxSize, -margin+yShift, -(boxSize+margin)+zShift), vec3(boxSize, -margin+yShift, (boxSize+margin)+zShift), vec3(boxSize, boxSize+yShift, (boxSize+margin)+zShift), vec3(boxSize, boxSize+yShift, -(boxSize+margin)+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(margin, 0.95, margin));
        }

        // Front Wall
        if(hit_quad(createQuad(vec3(-(boxSize+margin), -margin+yShift, -(boxSize+margin)+zShift), vec3((boxSize+margin), -margin+yShift, -boxSize+zShift), vec3((boxSize+margin), boxSize+yShift, -boxSize+zShift), vec3(-(boxSize+margin), boxSize+yShift, -boxSize+zShift)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.95));
        }

        // Sphere
        if(hit_sphere(
            createSphere(vec3(-1.0, 0.0, -3.0), 1.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDialectricMaterial(vec3(0.0, 1.0, 0.0), 1.1, 0.05);
        }

        // Sphere 2
        if(hit_sphere(
            createSphere(vec3(1.0, 0.0, -3.0), 1.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDialectricMaterial(vec3(1.0, 1.0, 0.0), 1.15, 0.05);
        }
    }

    if (CORNELL_BOX_BOXES) {

        if (hit_box(createBox(vec3(-2.5, 0.0+yShift, 0.0+zShift), vec3(0.5f, 5.0+yShift , 3.0+zShift)), r, tmin, rec.t, rec)) {
            hit = true;
            rec.material = createDialectricMaterial(vec3(0.0, 0.2, 0.8), 1.5, 0.05);
        }

        if (hit_box(createBox(vec3(-5.0f, 0.0+yShift, -5.0+zShift), vec3(5.5f, 2.0+yShift ,-3.0+zShift)), r, tmin, rec.t, rec)) {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.3, 0.1, 0.95));
        }

        if (hit_box(createBox(vec3(5.0f, 1.0+yShift, -1.0+zShift), vec3(5.5f, 7.0+yShift , 5.0+zShift)), r, tmin, rec.t, rec)) {
            hit = true;
            rec.material = createMetalMaterial(vec3(0.7, 0.6, 0.5), 0.0);
        }

    }

    if (CORNELL_BOX_SPHERES) {

        if(hit_sphere(
            createSphere(vec3(4.0f, 2.0+yShift, -2.0+zShift), 2.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createMetalMaterial(vec3(0.8, 0.1, 0.8), 0.05);
        }

        if(hit_sphere(
            createSphere(vec3(-.5f, 2.0+yShift, 2.0+zShift), 2.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDialectricMaterial(vec3(0.0, 0.2, 0.9), 1.5, 0.0);
        }

        if(hit_sphere(
            createSphere(vec3(-3.0f, 5.0+yShift, -1.0+zShift), 2.0),
            r,
            tmin,
            rec.t,
            rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.5, 0.6, 0.7));
        }

    }

    return hit;
}

vec3 directlighting(pointLight pl, Ray r, HitRecord rec){
    vec3 diffuseColor, specularColor;
    vec3 colorOut = vec3(0.0, 0.0, 0.0);
    float shininess;
    HitRecord dummy;
    vec3 emissive;

   vec3 lightPos;

    if (SOFT_SHADOWS)  // Soft shadows: random light source point from area light
    {
        vec3 a,b;
        if (CAUSTICS){
            a = vec3(0.7f, 0.0f, 0.0f);
            b = vec3(0.0f, 0.0f, 0.7f);
        }
        else{
            a = vec3(5.0f, 0.0f, 0.0f);
            b = vec3(0.0f, 0.0f, 5.0f);
        }

        lightPos = pl.pos + (hash1(gSeed) * a) + (hash1(gSeed) * b);
    }
    else if (SOFT_SHADOWS_GRID)  // Soft shadows: random light source point from N by N grid
    {
        float N = 5.0f;
        float q = floor(hash1(gSeed) * N);
        float p = floor(hash1(gSeed) * N);

        lightPos = pl.pos + vec3(1.0f, 0.0f, 0.0f) * p + vec3(0.0f, 0.0f, 1.0f) * q;
    }
    else {
        lightPos = pl.pos;
    }

    vec3 l = normalize(lightPos - rec.pos);
    float intensity = max(dot(rec.normal, l), .0f);
    Ray shadowRay = createRay(rec.pos + epsilon * rec.normal, l);
    float dist = length(pl.pos - rec.pos);
    
    if(!hit_world(shadowRay, 0.0f, dist, dummy) && intensity > .0f) {

        if(rec.material.type == MT_DIFFUSE) {
            diffuseColor = rec.material.albedo / pi * intensity;
            specularColor = rec.material.specColor;
            shininess = 4.0f / (pow(rec.material.roughness, 4.0f) + epsilon) - 2.0f;
            emissive = rec.material.emissive;
        }
        else if(rec.material.type == MT_METAL) {
            diffuseColor = rec.material.albedo;
            specularColor = metalSchlick(intensity, rec.material.specColor);
            shininess = 8.0f / (pow(rec.material.roughness, 4.0f) + epsilon) - 2.0f;
            emissive = rec.material.emissive;
        }
        else if(rec.material.type == MT_DIALECTRIC) {
            diffuseColor = rec.material.albedo;
            specularColor = rec.material.specColor;
            shininess = 500.0f;  
            emissive = rec.material.emissive;
        }

        vec3 h = normalize(l - r.d);
        colorOut = pl.color * diffuseColor + pl.color * specularColor * pow(max(dot(h, rec.normal), .0f), shininess) + emissive;
        
    }

	return colorOut;
}

#define MAX_BOUNCES 10

vec3 rayColor(Ray r)
{
    HitRecord rec;
    vec3 col = vec3(0.0);
    vec3 throughput = vec3(1.0f, 1.0f, 1.0f);

    for(int i = 0; i < MAX_BOUNCES; ++i)
    {
        if(hit_world(r, 0.001, 10000.0, rec))
        {
            //calculate direct lighting with 3 white point lights:
            if (dot(r.d, rec.normal) < 0.0f)
            {
                if (SHIRLEY) {
                    col += directlighting(createPointLight(vec3(-10.0, 15.0, 0.0), vec3(1.0, 1.0, 1.0)), r, rec) * throughput;
                    col += directlighting(createPointLight(vec3(8.0, 15.0, 3.0), vec3(1.0, 1.0, 1.0)), r, rec) * throughput;
                    col += directlighting(createPointLight(vec3(1.0, 15.0, -9.0), vec3(1.0, 1.0, 1.0)), r, rec) * throughput;
                }
                if (CORNELL_BOX_BOXES || CORNELL_BOX_SPHERES) {
                    float yShift = -3.5f;
                    float zShift = -8.0f;
                    float boxHeight = 8.0f;
                    col += directlighting(createPointLight(vec3(-0.1, boxHeight+yShift-.5, 0.1+zShift), vec3(1.0)), r, rec) * throughput;
                    col += directlighting(createPointLight(vec3(0.1, boxHeight+yShift-.5, -0.1+zShift), vec3(1.0)), r, rec) * throughput;
                    col += directlighting(createPointLight(vec3(0.1, boxHeight+yShift-.5, 0.1+zShift), vec3(1.0)), r, rec) * throughput;
}
                if (CAUSTICS){
                    float boxSize = 6.0f;
                    float yShift = -1.5f;
                    float zShift = -3.5f;
                    col += directlighting(createPointLight(vec3(0.0, boxSize+yShift, 0.0+zShift), vec3(1.0)), r, rec) * throughput;
                }

            }
           
            //calculate secondary ray and update throughput
            Ray scatterRay;
            vec3 atten;
            if(scatter(r, rec, atten, scatterRay))
            {  
                r = scatterRay;
                throughput *= atten;

            }
            else {
                return vec3(.0f);
            }
        
        }
        else  //background
        {
            float t = 0.8 * (r.d.y + 1.0);
            col += throughput * mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
            break;
        }
    }
    return col;
}

#define MAX_SAMPLES 10000.0

mat3 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

float isToggled(float keyCode) {
    return texelFetch(iKeyboard, ivec2(keyCode, 2), 0).x;
}

const float KEY_LEFT  = 37.0f;
const float KEY_UP    = 38.0f;
const float KEY_RIGHT = 39.0f;
const float KEY_DOWN  = 40.0f;

const ivec2 camTransforn = ivec2(0, 0);


void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

    vec2 mouse = iMouse.xy / iResolution.xy;
    mouse.x = mouse.x * 2.0 - 1.0;
    vec3 camPos = vec3(mouse.x * 10.0, mouse.y * 5.0, 8.0); 
    vec3 camTarget = vec3(0.0, 0.0, -1.0);
    float fovy = 60.0;
    float aperture = 0.0;
    float distToFocus = 1.0;
    float time0 = 0.0;
    float time1 = 1.0;
    
    vec3 pos2Target = camTarget - camPos;
    float zoom = texelFetch(iChannel0, camTransforn, 0).z; 
    float theta = texelFetch(iChannel0, camTransforn, 0).w;

    if(isToggled(KEY_UP) != .0f || isToggled(KEY_DOWN) != .0f) {
        
        zoom += iTimeDelta * ((isToggled(KEY_UP) != .0f)? 1.0f : -1.0f); // make zoom dependent on fps and not absolute time
        
        if (int(gl_FragCoord.x) == camTransforn.x && int(gl_FragCoord.y) == camTransforn.y) {
            gl_FragColor = vec4(texelFetch(iChannel0, camTransforn, 0).xy, zoom, texelFetch(iChannel0, camTransforn, 0).w);
        }
    }

    if(isToggled(KEY_LEFT) != .0f || isToggled(KEY_RIGHT) != .0f) {
        theta += sin(iTimeDelta) / 10.f * pi * ((isToggled(KEY_RIGHT) != .0f)? 1.0f : -1.0f); // make theta dependent on fps and not absolute time

        if (int(gl_FragCoord.x) == camTransforn.x && int(gl_FragCoord.y) == camTransforn.y) {
            gl_FragColor = vec4(texelFetch(iChannel0, camTransforn, 0).xyz, theta);
        }
    }
    
    camPos += zoom * normalize(pos2Target);

    pos2Target = camTarget - camPos; // changed due to zoom 

    vec3 rotated = pos2Target * rotateY(theta); // rotate pos2Target so cam orbits around the target and not the scene origin

    camPos = camTarget - rotated;

    Camera cam = createCamera(
        camPos,
        camTarget,
        vec3(0.0, 1.0, 0.0),    // world up vector
        fovy,
        iResolution.x / iResolution.y,
        aperture,
        distToFocus,
        time0,
        time1);

    //usa-se o 4 canal de cor para guardar o numero de samples e não o iFrame pois quando se mexe o rato faz-se reset

    vec4 prev = texture(iChannel0, gl_FragCoord.xy / iResolution.xy);
    vec3 prevLinear = toLinear(prev.xyz);  

    vec2 ps = gl_FragCoord.xy + hash2(gSeed);
    //vec2 ps = gl_FragCoord.xy;
    vec3 color = rayColor(getRay(cam, ps));


    if (int(gl_FragCoord.x) == camTransforn.x && int(gl_FragCoord.y) == camTransforn.y) {
        gl_FragColor = vec4(texelFetch(iChannel0, camTransforn, 0).xy, zoom, theta);
        return;
    }

    if(isToggled(KEY_LEFT) != .0f || isToggled(KEY_UP) != .0f || isToggled(KEY_RIGHT) != .0f || isToggled(KEY_DOWN) != .0f) {
        // clear frames so there's no distortion/blurriness from zoom and orbital rotation
        gl_FragColor = vec4(toGamma(color), 1.0);
        return;
    }    
    if(iMouseButton.x != 0.0 || iMouseButton.y != 0.0)
    {
        gl_FragColor = vec4(toGamma(color), 1.0);  //samples number reset = 1
        return;
    }
    if(prev.w > MAX_SAMPLES)   
    {
        gl_FragColor = prev;
        return;
    }

    float w = prev.w + 1.0;
    color = mix(prevLinear, color, 1.0/w);
    gl_FragColor = vec4(toGamma(color), w);
}
