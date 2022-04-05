#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoords;

uniform mat4 projection;
uniform mat4 view;
uniform vec2 move_to;
uniform mat4 model;

out vec2 TexCoords;

void main()
{
    TexCoords = aTexCoords;
    gl_Position = model * view * vec4(aPos.x+move_to.x, aPos.y+move_to.y, 0.0, 1.0);
}