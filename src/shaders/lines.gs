#version 330 core
layout (points) in;
layout (triangle_strip, max_vertices = 24) out;

in VS_OUT {
    vec3 color;
} gs_in[];

out vec3 fColor;
#define NR_CELLS_MAX 500
uniform float cell_size;
uniform float line_thickness;
uniform bool N; //0
uniform bool S; //1
uniform bool E; //2
uniform bool W; //3
uniform vec2 crd;
uniform bool trace_path;
uniform int drc; //03,02,01,12,13,23
uniform vec2 scrRes;
uniform vec2 texRes;

uniform int nr_cells;

void build_grid(vec4 position)
{
        fColor = gs_in[0].color;
        float dx = crd.x;
        float dy = crd.y;
        float start_x = cell_size*dx;
        float start_y = -cell_size*dy;
        float start_y_b = -cell_size*(dy+1);
        float start_x_r = cell_size*(dx+1);
        if(!trace_path){
            if(N){
                gl_Position = position + vec4(start_x/scrRes.x, (-(line_thickness/2)+start_y)/scrRes.y, 0.0, 0.0); // 1:bottom-left
                EmitVertex();
                gl_Position = position + vec4(start_x_r, -(line_thickness/2)+start_y, 0.0, 0.0); // 2:bottom-right
                EmitVertex();
                gl_Position = position + vec4(start_x, (line_thickness/2)+start_y, 0.0, 0.0); // 3:top-left
                EmitVertex();
                gl_Position = position + vec4(start_x_r, (line_thickness/2)+start_y, 0.0, 0.0); // 4:top-right
                EmitVertex();
                EndPrimitive();
            }
            if(S){
                gl_Position = position + vec4(start_x, -(line_thickness/2)+start_y_b, 0.0, 0.0); // 1:bottom-left
                EmitVertex();
                gl_Position = position + vec4(start_x_r, -(line_thickness/2)+start_y_b, 0.0, 0.0); // 2:bottom-right
                EmitVertex();
                gl_Position = position + vec4(start_x, (line_thickness/2)+start_y_b, 0.0, 0.0); // 3:top-left
                EmitVertex();
                gl_Position = position + vec4(start_x_r, (line_thickness/2)+start_y_b, 0.0, 0.0); // 4:top-right
                EmitVertex();
                EndPrimitive();
            }
            if(E){
                gl_Position = position + vec4(-(line_thickness/2)+start_x_r, start_y_b, 0.0, 0.0); // 1:bottom-left
                EmitVertex();
                gl_Position = position + vec4( (line_thickness/2)+start_x_r, start_y_b, 0.0, 0.0); // 2:bottom-right
                EmitVertex();
                gl_Position = position + vec4(-(line_thickness/2)+start_x_r, start_y, 0.0, 0.0); // 3:top-left
                EmitVertex();
                gl_Position = position + vec4( (line_thickness/2)+start_x_r, start_y, 0.0, 0.0); // 4:top-right
                EmitVertex();
                EndPrimitive();
            }
            if(W){
                gl_Position = position + vec4(-(line_thickness/2)+start_x, start_y_b, 0.0, 0.0); // 1:bottom-left
                EmitVertex();
                gl_Position = position + vec4( (line_thickness/2)+start_x, start_y_b, 0.0, 0.0); // 2:bottom-right
                EmitVertex();
                gl_Position = position + vec4(-(line_thickness/2)+start_x, start_y, 0.0, 0.0); // 3:top-left
                EmitVertex();
                gl_Position = position + vec4( (line_thickness/2)+start_x, start_y, 0.0, 0.0); // 4:top-right
                EmitVertex();
                EndPrimitive();
            }
        } else {
            switch (drc) {
                case 3: //NW,WN
                    gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y, 0.0, 0.0); // 1:top-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y, 0.0, 0.0); // 2:top-right
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 3:bottom-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 4:bottom-right
                    EmitVertex();
                    EndPrimitive();

                    gl_Position = position + vec4(start_x, start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 1:top-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x, start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 2:bottom-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 3:top-right
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 4:bottom-right
                    EmitVertex();
                    EndPrimitive();

                    break;
                case 2: //NE,EN
                    gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y, 0.0, 0.0); // 1:top-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y, 0.0, 0.0); // 2:top-right
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 3:bottom-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 4:bottom-right
                    EmitVertex();
                    EndPrimitive();

                    gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 1:top-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 2:bottom-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x_r, start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 3:top-right
                    EmitVertex();
                    gl_Position = position + vec4(start_x_r, start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 4:bottom-right
                    EmitVertex();
                    EndPrimitive();

                    break;
                case 1: //NS,SN
                    gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y, 0.0, 0.0); // 1:top-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y, 0.0, 0.0); // 2:top-right
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y_b, 0.0, 0.0); // 3:bottom-left
                    EmitVertex();
                    gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y_b, 0.0, 0.0); // 4:bottom-right
                    EmitVertex();

                    break;
                case 12: //SE,ES
                     gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 1:top-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 2:top-right
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y_b, 0.0, 0.0); // 3:bottom-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y_b, 0.0, 0.0); // 4:bottom-right
                     EmitVertex();
                     EndPrimitive();

                     gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 1:top-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 2:bottom-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x_r, start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 3:top-right
                     EmitVertex();
                     gl_Position = position + vec4(start_x_r, start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 4:bottom-right
                     EmitVertex();
                     EndPrimitive();

                     break;
                case 13: //SW,WS
                     gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 1:top-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y-(cell_size/2), 0.0, 0.0); // 2:top-right
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2)-(line_thickness/2), start_y_b, 0.0, 0.0); // 3:bottom-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2)+(line_thickness/2), start_y_b, 0.0, 0.0); // 4:bottom-right
                     EmitVertex();
                     EndPrimitive();

                     gl_Position = position + vec4(start_x, start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 1:top-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x, start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 2:bottom-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 3:top-right
                     EmitVertex();
                     gl_Position = position + vec4(start_x+(cell_size/2), start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 4:bottom-right
                     EmitVertex();
                     EndPrimitive();

                     break;
                case 23: //EW,WE
                     gl_Position = position + vec4(start_x, start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 1:top-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x, start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 2:bottom-left
                     EmitVertex();
                     gl_Position = position + vec4(start_x_r, start_y-(cell_size/2)+(line_thickness/2), 0.0, 0.0); // 3:top-right
                     EmitVertex();
                     gl_Position = position + vec4(start_x_r, start_y-(cell_size/2)-(line_thickness/2), 0.0, 0.0); // 4:bottom-right
                     EmitVertex();
                     EndPrimitive();

                     break;
            }
        }

   fColor = vec3(1.0, 1.0, 1.0);
}

void main() {
    build_grid(gl_in[0].gl_Position);
}