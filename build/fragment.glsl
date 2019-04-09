#version 330 core

#define float2 vec2
#define float3 vec3
#define float4 vec4
#define float4x4 mat4
#define float3x3 mat3

#define MAX_SCENE_DISTANCE 1e2
#define MAX_MARCHES_NUM 100
#define LIGHTS_NUM 2
#define GROUND_LEVEL -1.01f
#define MAX_REFLECTIONS 1

in float2 fragmentTexCoord;

layout(location = 0) out vec4 fragColor;

uniform int g_screenWidth;
uniform int g_screenHeight;

uniform float4x4 g_rayMatrix;

uniform int g_softShadow;

vec3 lights[2];

struct Hit {
    float dist;
    vec3 amb, diff, spec;
    float spec_exp;
};

float sdSphere(vec3 p, float s) {
    return length(p)-s;
}

float sdBox(vec3 p, vec3 b) {
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float maxXY (vec2 p) {
    return max(p.x,p.y);
}

float sdCross (vec3 p) {
    float da = maxXY(abs(p.xy));
    float db = maxXY(abs(p.yz));
    float dc = maxXY(abs(p.zx));
    return min(da,min(db,dc))-1.0;
}

float Menger (vec3 p) {
    float d = sdBox(p,vec3(1.0));

    float m = 1.0;
    for (int it = 0; it < 4; it++) {
        vec3 a = mod(p * m, 2.0) -1.0;
        m *= 3.0;
        vec3 r = 1.0 - 3.0 * abs(a);

        float c = sdCross(r) / m;
        d = max(d,c);
    }

    return d;
}

float Mandelbulb (vec3 p) {
    vec3 cur = p;
    float rad1 = 0.0;
    float hit = 0.0;
    float c = 16.0;
    float d = 1.0;
    for (int i = 0; i < 3; i++) {
        rad1 = length(cur);

        if (rad1 > 3.0) {   
            hit = 0.25 * log(rad1) * rad1 / d;
        }
        else {
            float th = atan(length(cur.xy), cur.z);
            float phi = atan(cur.y, cur.x);     
            float rad2 = pow(rad1,8.0);
            d = pow(rad1, 7.0) * 7.0 * d + 1.0;
            


            float sint = sin(th * c);
            cur.x = rad2 * sint * cos(phi * c);
            cur.y = rad2 * sint * sin(phi * c);
            cur.z = rad2 * cos(th * c) ;
            cur += p;
        }
    }
    return hit;
}

vec4 squareQuaternion (vec4 a) {
    return vec4( a.x*a.x - a.y*a.y - a.z*a.z - a.w*a.w,
                 2.0*a.x*a.y,
                 2.0*a.x*a.z,
                 2.0*a.x*a.w );
}

float JuliaQuaternion (vec3 p, vec4 c) {
    vec4 z = vec4(p,0.0);
    float md2 = 1.0;
    float mz2 = dot(z,z);

    for( int i=0; i<5; i++ )
    {
        md2 *= 4.0*mz2;
        z = squareQuaternion(z) + c;

        mz2 = dot(z,z);
        if(mz2>4.0) break;
    }

    return 0.25*sqrt(mz2/md2)*log(mz2);
}

Hit leastHit (Hit a, Hit b) {
    return a.dist < b.dist ? a : b;
}

Hit leastDistance (vec3 pos) {
    Hit res = Hit (MAX_SCENE_DISTANCE, vec3(0), vec3(0), vec3(0), 0.0);

    res = leastHit(res, Hit (sdSphere(pos-vec3(1.25,-0.5,-1.5), 0.5),
                             vec3(0.24725,0.1995,0.0745), vec3(0.75164,0.60648,0.22648), vec3(0.628281,0.555802,0.366065), 0.4)); //gold sphere
    res = leastHit(res, Hit (sdBox(pos-vec3(-1.25,-0.5,-1.5), vec3(0.5,0.5,0.2)),
                             vec3(0.0215,0.1745,0.0215), vec3(0.07568,0.61424,0.07568), vec3(0.633,0.727811,0.633), 0.6)); //emerald box
    res = leastHit(res, Hit (Menger(pos - vec3(3.0,0.0,2.0)),
                             vec3(0.25), vec3(0.4), vec3(0.774597), 0.6));
    res = leastHit(res, Hit (Mandelbulb(pos - vec3(-2.5,0.1,2.0)),
                             vec3(0.05375,0.05,0.06625), vec3(0.18275,0.17,0.22525), vec3(0.332741,0.328634,0.346435), 0.6));
    res = leastHit(res, Hit (JuliaQuaternion(pos - vec3(0.0,0.15,5.5), vec4(0.2,0.3,0.4,0.5)),
                             vec3(0.1745,0.01175,0.01175), vec3(0.61424,0.04136,0.04136), vec3(0.727811,0.626959,0.626959), 0.6));

    return res;
}

Hit castRay (vec3 ray_pos, vec3 ray_dir) {
    Hit res = Hit (-1, vec3(0.5, 0.0, 0.5), vec3(0), vec3(0), 0.0);
    
    float tmin = 1e-2;
    float tmax = MAX_SCENE_DISTANCE;

    float gl = (GROUND_LEVEL - ray_pos.y) / ray_dir.y;
    if (gl > 0.0) {
        tmax = min (tmax, gl); // Эта строчка должна уменьшать область поиска пересечения, но она портит тени
        res = Hit (-1, vec3(0.3,0.4,0.5), vec3(0.3), vec3(0.0), 0.0);
    }

    float t = max(leastDistance(ray_pos).dist, tmin);

    for (int i = 0; i < MAX_MARCHES_NUM && t < tmax; i++) {
        Hit cur = leastDistance(ray_pos + t*ray_dir);

        if  (abs(cur.dist) < 1e-4*t) {
            res = cur;
            res.dist += t;
            break;
        }

        t += cur.dist;
    }

    return res;
}

vec3 normalToPoint(vec3 pos) {
    vec2 my_eps = vec2(1.0,-1.0)*0.5683*0.0005;
    return normalize( my_eps.xyy*leastDistance( pos + my_eps.xyy ).dist + 
                      my_eps.yyx*leastDistance( pos + my_eps.yyx ).dist + 
                      my_eps.yxy*leastDistance( pos + my_eps.yxy ).dist + 
                      my_eps.xxx*leastDistance( pos + my_eps.xxx ).dist );
}

float mySmoothSpec (float a, float spec_exp) {
    return (a > 0) ? pow (a, 128 * spec_exp) : 0.0;
}

float mySmoothDiff (float a) {
    return (a > 0) ? pow (a, 2) : 0.0;
}

/*bool isFurtherThenLight (float hit_dist, float light_dist) {
    return hit_dist < 0 || hit_dist > light_dist;
}*/

vec3 diffAndSpecLightFrom (vec3 hit_pos, vec3 light_pos, vec3 ray_dir, vec3 n, vec3 diff, vec3 spec, float spec_exp) {
    vec3 dir_to_light = normalize(light_pos - hit_pos);
    vec3 r = 2 * dot(dir_to_light, n) * n - dir_to_light;
    return diff * mySmoothDiff(dot(dir_to_light, n)) + spec * mySmoothSpec(dot(r, (-1) * ray_dir), spec_exp);

    return vec3(0.0);
}

float castShadow (vec3 ray_pos, vec3 ray_dir, float tmin, float tmax) {
    for (float t = tmin; t < tmax; ) {
        float d = leastDistance(ray_pos + ray_dir * t).dist;
        if(d < 0.001)
            return 0.0;

        t += d;
    }
    return 1.0;
}

float castSoftShadow (vec3 ray_pos, in vec3 ray_dir, float tmin, float tmax, float p) {
    float res = 1.0;
    for (float t = tmin; t < tmax; ) {
        float d = leastDistance(ray_pos + ray_dir * t).dist;
        if(d < 0.001)
            return 0.0;
        
        res = min(res, p * d / t);
        t += d;
    }
    return clamp (res, 0.0, 1.0);
}

vec3 constrainToLegitColor (vec3 v) {
    return vec3(clamp(v.x, 0.0, 1.0), clamp(v.y, 0.0, 1.0), clamp(v.z, 0.0, 1.0));
}

vec3 render (vec3 ray_pos, vec3 ray_dir) {
    
    Hit first_hit = castRay(ray_pos, ray_dir);
    vec3 light = first_hit.amb;

    vec3 first_hit_point = ray_pos + ray_dir * first_hit.dist;

    float gl = (GROUND_LEVEL - ray_pos.y) / ray_dir.y;
    if (gl > 0.0 && first_hit.dist == -1) {
        first_hit.dist = gl;
        first_hit_point = ray_pos + ray_dir * gl;
        for (int i = 0; i < LIGHTS_NUM; i++) {
            if (g_softShadow == 0) 
                light += first_hit.diff * castShadow(first_hit_point, normalize(lights[i] - first_hit_point), 1e-2, 2.5);
            else 
                light += first_hit.diff * castSoftShadow(first_hit_point, normalize(lights[i] - first_hit_point), 1e-2, 2.5, 8.0);
        }
        return vec3(constrainToLegitColor(light));
    }
    else if (first_hit.dist == -1) {
        return first_hit.amb;
    }

    for (int i = 0; i < LIGHTS_NUM; i++) {
        light += diffAndSpecLightFrom(first_hit_point, lights[i], ray_dir, normalToPoint(first_hit_point), first_hit.diff, first_hit.spec, first_hit.spec_exp) 
                 * castSoftShadow(first_hit_point, normalize(lights[i] - first_hit_point), 1e-2, 2.5, 8.0);
    }

    vec3 reflect_pos = first_hit_point, reflect_dir = reflect(ray_dir, normalToPoint(reflect_pos)); 
    for (int i = 0; i < MAX_REFLECTIONS; i++) {
        Hit reflect_hit = castRay(reflect_pos, reflect_dir);
        vec3 reflect_hit_point = reflect_pos + reflect_dir * reflect_hit.dist;
        gl = (GROUND_LEVEL - reflect_pos.y) / reflect_dir.y;

        for (int j = 0; j < LIGHTS_NUM; j++) {
            if (gl > 0.0 && reflect_hit.dist == -1)
                light += diffAndSpecLightFrom(reflect_hit_point, lights[j], reflect_dir, vec3(0.0,1.0,0.0), reflect_hit.diff, reflect_hit.spec, reflect_hit.spec_exp);
            else
                light += diffAndSpecLightFrom(reflect_hit_point, lights[j], reflect_dir, normalToPoint(reflect_hit_point), reflect_hit.diff, reflect_hit.spec, reflect_hit.spec_exp);
        }
    }

    return constrainToLegitColor(light);
}

float3 EyeRayDir(float x, float y, float w, float h) {
    float fov = 3.141592654f/(2.0f); 
    float3 ray_dir;
  
    ray_dir.x = x+0.5f - (w/2.0f);
    ray_dir.y = y+0.5f - (h/2.0f);
    ray_dir.z = -(w)/tan(fov/2.0f);
    
    return normalize(ray_dir);
}

void main(void) {

    float w = float(g_screenWidth);
    float h = float(g_screenHeight);
  
    // get curr pixelcoordinates
    //
    float x = fragmentTexCoord.x*w; 
    float y = fragmentTexCoord.y*h;
  
    // generate initial ray
    //
    float3 ray_pos = float3(0,0,0); 
    float3 ray_dir = EyeRayDir(x,y,w,h);
 
    // transorm ray with matrix
    //
    ray_pos = (g_rayMatrix*float4(ray_pos,1)).xyz;
    ray_dir = float3x3(g_rayMatrix)*ray_dir;

    lights[0] = vec3(10,12,0);
    lights[1] = vec3(0,14,8);
  
    //float tmin = 0, tmax = 0;
    
    // //Начало самодеятельности

    fragColor = vec4(render (ray_pos, ray_dir), 1.0);
    //fragColor = vec4(abs(ray_dir.x),abs(ray_dir.y),abs(ray_dir.x),1) ;
}