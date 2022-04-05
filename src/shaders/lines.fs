#version 330 core
out vec4 FragColor;

in vec3 fColor;
uniform bool trace_path;

void main()
{
    if(trace_path){
       FragColor = vec4(fColor, 1.0);
    } else {
       FragColor = vec4(vec3(0.7,0.85,1.0), 1.0);
    }
}