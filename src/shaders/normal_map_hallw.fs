#version 330 core
out vec4 FragColor;

in VS_OUT {
    vec3 FragPos;
    vec2 TexCoords;
    vec3 TangentLightPos;
    vec3 TangentViewPos;
    vec3 TangentFragPos;
} fs_in;

uniform sampler2D diffuseMap;
uniform float shininess;

//uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec2 scrRes;
uniform vec2 texRes;

#define NR_POINT_LIGHTS_MAX 70
#define OFFSET_X 1
#define OFFSET_Y 1
#define DEPTH	 15.5

uniform int nr_point_lights;
uniform int utime;

struct Light {
    vec3 position;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    float constant;
    float linear;
    float quadratic;
};

struct PointLight {
    vec3 position;

    float constant;
    float linear;
    float quadratic;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

struct SpotLight {
    vec3 position;
    vec3 direction;
    float cutOff;
    float outerCutOff;

    float constant;
    float linear;
    float quadratic;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

struct DirLight {
    vec3 direction;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

uniform PointLight pointLights[NR_POINT_LIGHTS_MAX];
uniform SpotLight spotLight;
uniform DirLight dirLight;

// function prototypes
vec3 CalcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir);
vec3 CalcSpotLight(SpotLight light, vec3 normal, vec3 fragPos, vec3 viewDir);
vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir);

// The dX AND DY WHEN DIFFERENTIATING HEIGHT VALUES.
#define DIFF 0.001
// how FAST THE LIGHT MOVES.
#define TIMESCALE 1.5
// the RADIUS OF THE LIGHT'S PATH.
#define LIGHTPATHRADIUS 0.35
// The center of the light's path.
#define LIGHTPATHCENTER vec3(0.5, 0.5, 0.125)
// Color of the light.
#define LIGHTCOLOR vec4(1.0, 1.0, 1.0, 1.0)
// Light strength multiplier.
#define LIGHTSTRENGTH 2.0
// The implied Z position of the lit surface.
#define SURFACEZDEPTH 1.0
// The light's ambient term.
#define AMBIENTCOLOR vec4(0.0, 0.0, 0.0, 1.0)
// The falloff factor of the specular lighting term.
#define SPECULARFACTOR 8.0
// The coefficient to the specular factor when negotiating
// the brightness of distant, but satisfactorily oriented bumps.
#define SPECULARRATIO 5.0
// The coefficient of the values given by the specular map.
#define SPECULARMAPRATIO 8.0
// Whether or not to use a texture as the base value.

float getHeightValue(sampler2D tex, vec2 coord)
{
	return texture(tex, coord).r;
}


vec2 getLocalDiff(sampler2D tex, vec2 coord)
{
	// Get the local difference of height along the X axis.
	float diffX = getHeightValue(tex, vec2(coord.x+DIFF, coord.y) )
		-getHeightValue(tex, vec2(coord.x-DIFF, coord.y) );

	// Do the same along the Y axis.
	float diffY = getHeightValue(tex, vec2(coord.x, coord.y+DIFF) )
		-getHeightValue(tex, vec2(coord.x, coord.y-DIFF) );

	// Return the two values as a 2D vector.
	return vec2(diffX, diffY);
}


vec3 getSurfaceNormal(sampler2D tex, vec2 coord)
{
	// Get the local difference in height about the coordinate given.
	vec2 localDiff = getLocalDiff(tex, coord);

	// Remember that the surface normal is a negative reciprocal of
	// the surface tangent (which is what the local difference really is).
	// This step does half that job, negating the local difference.
	localDiff *= -1.0;

	// Remember that this is to be stored in a pixel, so we have to
	// fit it to the range [0..1].
	localDiff = (localDiff/2.0)+.5;
//	localDiff = normalize(localDiff * 2.0 - 1.0);

	float localDiffMag = length(localDiff);
	float z = sqrt(1.0-pow(localDiffMag, 2.0));

	return vec3(localDiff, z);
}

float getSpecularity(sampler2D tex, vec2 coords)
{
	return texture(tex, coords).b*.5;
}

//vec3 genLightCoords()
//{
//	// Use simple trig to rotate the light position around a point.
//	vec3 lightCoords = vec3(LIGHTPATHCENTER.x + (sin((utime/10)*TIMESCALE)*LIGHTPATHRADIUS),
//				LIGHTPATHCENTER.y + (cos((utime/10)*TIMESCALE)*LIGHTPATHRADIUS),
//				LIGHTPATHCENTER.z);
//	return lightCoords;
//}

//the vector of incidence between a light position and a surface position.
vec3 getIncidence(vec3 lightPos, vec2 coord)
{
	// To get the incidence vector we subtract the final position from the original
	// position. This gives us a vector pointing into the surface.
	return lightPos - vec3(coord, SURFACEZDEPTH);
}

//the cosine of the angle of incidence of our light and the *flat* surface.
float getAngle(vec3 lightIncidence, vec3 normal)
{
	// We have to unpack the normal vector.
	normal.xy -= .5;
	normal.xy *= 2.0;

	// Normalize the two participating vectors so we don't get
	// strange results.
	normal = normalize(normal);
	lightIncidence = normalize(lightIncidence);

	// Return the dot product of the two, which represents the cosine of the angle
	// between the two vectors.
	return dot(lightIncidence, normal);
}

float getDist(vec3 light, vec2 coord)
{
	return distance(light, vec3(coord, 0.0));
}

vec4 getLighting(sampler2D heightTex, sampler2D specTex, vec2 coord)
{
	// Get the current light position.
	//vec3 lightPos = genLightCoords();
    vec3 lightPos = pointLights[0].position;

	// Get the vector of incidence the light has with the curren texel.
	vec3 lightIncidence = getIncidence(lightPos, coord);

	// Also get the surface normal of the current texel.
	vec3 surfaceNormal = getSurfaceNormal(heightTex, coord);

	// Determine the cosine of the angle between the incident and normal vectors.
	float cosine = getAngle(lightIncidence, surfaceNormal);

	// Also get the distance from the light to the current texel, for
	// distance falloff.
	float dist = getDist(lightPos, coord);


	vec4 ambient = AMBIENTCOLOR;

	// Create a linear-falloff diffuse light term.
	vec4 diffuse = vec4(1.0);
	diffuse *= LIGHTSTRENGTH;
	diffuse *=  (1.0-dist);
	diffuse *= cosine;
	diffuse *= LIGHTCOLOR;

	// Get the local specularity (shininess of the material.
	float spec = getSpecularity(specTex, coord);

	// Create a powered-falloff specular term.
	vec4 specular = vec4(1.0);
	specular *= LIGHTSTRENGTH;
	specular *=  pow((1.0-dist), SPECULARFACTOR);
	specular *= pow(cosine, SPECULARFACTOR*SPECULARRATIO);
	specular *= LIGHTCOLOR;
	specular *= spec*SPECULARMAPRATIO;

	return ambient+diffuse+specular;
}



void main()
{

        FragColor = vec4(1.0); // set all 4 vector values to 1.0
        // obtain normal from normal map in range [0,1]
        vec3 normal = getSurfaceNormal(diffuseMap, fs_in.TexCoords);
        // transform normal vector to range [-1,1]
        normal = normalize(normal * 2.0 - 1.0);  // this normal is in tangent space

        vec3 viewDir = normalize(fs_in.TangentViewPos - fs_in.TangentFragPos);

        vec3 result = vec3(0.0);

        result = CalcDirLight(dirLight, normal, viewDir);
        //for(int i = 0; i < nr_point_lights; i++)
        for(int i = 0; i < nr_point_lights; i++){
           result += CalcPointLight(pointLights[i], normal, fs_in.TangentFragPos, viewDir);
        }
        FragColor = vec4(result, 1.0);

}

// calculates the color when using a directional light.
vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir)
{
    vec3 lightDir = normalize(-light.direction);
    // combine results
    vec3 color = texture(diffuseMap, fs_in.TexCoords).rgb;

    vec3 spec_color = normal;
    // ambient
    vec3 ambient = light.ambient * color;

    //vec3 ambient = light.ambient * vec3(texture(light.diffuse, fs_in.TexCoords));
    //vec3 diffuse = light.diffuse * diff * vec3(texture(light.diffuse, fs_in.TexCoords));
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = light.diffuse * diff * color;

    //vec3 specular = light.specular * spec * vec3(texture(light.specular, fs_in.TexCoords));
    vec3 reflectDir = reflect(-lightDir, normal);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), shininess);
    vec3 specular = light.specular * spec_color * spec;

    return (ambient + diffuse + specular);
}
// calculates the color when using a point light.
vec3 CalcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir)
{

    // get diffuse color
    vec3 color = texture(diffuseMap, fs_in.TexCoords).rgb;

    vec3 spec_color = normal;
    // ambient
    vec3 ambient = light.ambient * color;
    // diffuse
    //fs_in.TangentLightPos = light.position;
    //vec3 lightDir = normalize(fs_in.TangentLightPos - fs_in.TangentFragPos);
    vec3 lightDir = normalize(light.position - fs_in.TangentFragPos);
    float diff = max(dot(lightDir, normal), 0.0);
    vec3 diffuse = light.diffuse * diff * color;
    // specular
    //vec3 viewDir = normalize(fs_in.TangentViewPos - fs_in.TangentFragPos);
    vec3 reflectDir = reflect(-lightDir, normal);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), shininess);
    vec3 specular = light.specular * spec_color * spec;

    float distance    = length(light.position - fs_in.TangentFragPos);
    float attenuation = 1.0 / (light.constant + light.linear * distance + light.quadratic * (distance * distance));

    ambient  *= attenuation;
    diffuse   *= attenuation;
    specular *= attenuation;

    return (ambient + diffuse + specular);
}

// calculates the color when using a spot light.
vec3 CalcSpotLight(SpotLight light, vec3 normal, vec3 fragPos, vec3 viewDir)
{
     // get diffuse color
    vec3 color = texture(diffuseMap, fs_in.TexCoords).rgb;
    vec3 spec_color = normal;
    // ambient
    vec3 ambient = light.ambient * color;
    // diffuse
    //fs_in.TangentLightPos = light.position;
    //vec3 lightDir = normalize(fs_in.TangentLightPos - fs_in.TangentFragPos);
    vec3 lightDir = normalize(light.position - fs_in.TangentFragPos);
    float diff = max(dot(lightDir, normal), 0.0);
    vec3 diffuse = light.diffuse * diff * color;
    // specular
    //vec3 viewDir = normalize(fs_in.TangentViewPos - fs_in.TangentFragPos);
    vec3 reflectDir = reflect(-lightDir, normal);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), shininess);
    vec3 specular = light.specular * spec_color * spec;

    // attenuation
    float distance = length(light.position - fragPos);
    float attenuation = 1.0 / (light.constant + light.linear * distance + light.quadratic * (distance * distance));
    // spotlight intensity
    float theta = dot(lightDir, normalize(-light.direction));
    float epsilon = light.cutOff - light.outerCutOff;
    float intensity = clamp((theta - light.outerCutOff) / epsilon, 0.0, 1.0);
    // combine results
   // vec3 ambient = light.ambient * vec3(texture(light.diffuse, fs_in.TexCoords));
   // vec3 diffuse = light.diffuse * diff * vec3(texture(light.diffuse, fs_in.TexCoords));
   // vec3 specular = light.specular * spec * vec3(texture(light.specular, fs_in.TexCoords));


    ambient *= attenuation * intensity;
    diffuse *= attenuation * intensity;
    specular *= attenuation * intensity;
    return (ambient + diffuse + specular);
}
