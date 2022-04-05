#define GLM_ENABLE_EXPERIMENTAL
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "filesystem.h"
#include "shader.h"
#include "camera.h"
#include "model.h"

#include <iostream>
#include <fstream>
#include <numeric>
#include <chrono>
#include <set>
#include <algorithm>
#include <random>

enum EdgeType{
    leftedge = 0,
    rightedge = 1,
    topedge = 2,
    bottomedge = 3,
    corner_tl = 6,
    corner_tr = 7,
    corner_bl = 8,
    corner_br = 9,
    notedge = 10,
    unset = 11
};

enum Shape{
    Triangle = 0,
    Quad = 1
};

typedef std::set<Direction> Cell;
typedef std::vector<std::vector<Cell>> Grid;
typedef std::pair<int ,int> Coord;
typedef std::map<Coord,EdgeType> EdgeTMap;
typedef std::pair<std::vector<std::pair<Coord,Direction>>,int> PathTo;
typedef std::map<Coord, PathTo> PathMp;


void get_texture_dimensions(const std::string& file_path, unsigned int tex, std::map<unsigned int, glm::vec2>& tex_dims)
{
    unsigned char buf[8];

    std::ifstream in(FileSystem::getPath(file_path));
    in.seekg(16);
    in.read(reinterpret_cast<char*>(&buf), 8);

    unsigned int width; unsigned int height;
    width = (buf[0] << 24) + (buf[1] << 16) + (buf[2] << 8) + (buf[3] << 0);
    height = (buf[4] << 24) + (buf[5] << 16) + (buf[6] << 8) + (buf[7] << 0);

    tex_dims[tex] = glm::vec2(width, height);
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

bool mouse_disabled = false;

//void processInput(GLFWwindow *window, PathTo& path, Grid& grid);
void processInput(GLFWwindow *window);
unsigned int loadTexture(const char *path);
unsigned int renderQuad();
unsigned int loadCubemap(vector<std::string> faces);

std::map<Direction,Direction> opposite_drc{{N,S},{S,N},{E,W},{W,E}};

template<typename T>
void shuffle(std::vector<T>& v){
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(v.begin(), v.end(), g);
}

void make_path(Coord& this_crd, std::vector<std::pair<Coord, Direction>>& next,
               Grid& grid, std::set<Coord>& added, PathMp& pathmap) {
    shuffle(next);
    bool dead_end = true;
    for(auto& c : next) {
        if(!added.contains(c.first)) {
            if(pathmap[c.first].first.empty()) {
            //    pathmap[c.first].first = std::vector<Coord>(pathmap[this_crd].first);
                pathmap[c.first].first = pathmap[this_crd].first;
               // pathmap[c.first].first.push_back(c);
               // pathmap[c.first].second = pathmap[this_crd].second+1;
            } //else {
            pathmap[c.first].first.push_back(c);
            pathmap[c.first].second = pathmap[this_crd].second+1;
            //}
            added.insert(c.first);
            switch(c.second){
                case N : {
                    grid[this_crd.first][this_crd.second].erase(N);
                    grid[c.first.first][c.first.second].erase(S);
                    break;
                }
                case S : {
                    grid[this_crd.first][this_crd.second].erase(S);
                    grid[c.first.first][c.first.second].erase(N);
                    break;
                }
                case E : {
                    grid[this_crd.first][this_crd.second].erase(E);
                    grid[c.first.first][c.first.second].erase(W);
                    break;
                }
                case W : {
                    grid[this_crd.first][this_crd.second].erase(W);
                    grid[c.first.first][c.first.second].erase(E);
                    break;
                }
            }
            this_crd = c.first;
            dead_end = false;
            break;
        }
    }
    if(dead_end){
        if(pathmap[this_crd].first.size() >= 2){
           int last_idx = pathmap[this_crd].first.size() - 2;
           this_crd = pathmap[this_crd].first[last_idx].first;
        } else {
            int foo = 2;
        }
    }
}

void find_next_path(Coord& this_crd, PathMp& pathmap, Grid& grid, std::set<Coord>& added,
                    EdgeTMap& edge_map, int size){

    int col = std::get<0>(this_crd);
    int row = std::get<1>(this_crd);
    std::vector<std::pair<Coord,Direction>> next;
    EdgeType edge_type = edge_map[this_crd];
    switch(edge_type) {
        case leftedge : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row+1},S));
            next.emplace_back(std::make_pair<Coord,Direction>({col,row-1},N));
            next.emplace_back(std::make_pair<Coord,Direction>({col+1,row},E));
            make_path(this_crd, next,grid, added, pathmap);
            break;
        }
        case rightedge : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row+1},S));
            next.emplace_back(std::make_pair<Coord,Direction>({col,row-1},N));
            next.emplace_back(std::make_pair<Coord,Direction>({col-1,row},W));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        case topedge : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row+1},S));
            next.emplace_back(std::make_pair<Coord,Direction>({col-1,row},W));
            next.emplace_back(std::make_pair<Coord,Direction>({col+1,row},E));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        case bottomedge : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row-1},N));
            next.emplace_back(std::make_pair<Coord,Direction>({col-1,row},W));
            next.emplace_back(std::make_pair<Coord,Direction>({col+1,row},E));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        case corner_tl : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row+1},S));
            next.emplace_back(std::make_pair<Coord,Direction>({col+1,row},E));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        case corner_tr : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row+1},S));
            next.emplace_back(std::make_pair<Coord,Direction>({col-1,row},W));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        case corner_bl : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row-1},N));
            next.emplace_back(std::make_pair<Coord,Direction>({col+1,row},E));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        case corner_br : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row-1},N));
            next.emplace_back(std::make_pair<Coord,Direction>({col-1,row},W));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        case notedge : {
            next.emplace_back(std::make_pair<Coord,Direction>({col,row-1},N));
            next.emplace_back(std::make_pair<Coord,Direction>({col,row+1},S));
            next.emplace_back(std::make_pair<Coord,Direction>({col-1,row},W));
            next.emplace_back(std::make_pair<Coord,Direction>({col+1,row},E));
            make_path(this_crd, next, grid, added, pathmap);
            break;
        }
        default:;
    }

    if(added.size() < size){
        find_next_path(this_crd, pathmap, grid, added, edge_map, size);
    }
}

void getEdges(EdgeTMap& edge_map, int size){
    for(int r = 0; r < size; ++r){
        for(int c = 0; c < size; ++c){
            Coord this_crd{c,r};
            bool notedge_h = false;
            bool notedge_v = false;

            EdgeType edge_type = unset;
            if(c == 0) {
                if(r == 0){
                    edge_type = corner_tl;
                } else if(r == size-1) {
                    edge_type = corner_bl;
                } else {
                    edge_type = leftedge;
                }
            } else if(c == size-1) {
                if(r == 0){
                    edge_type = corner_tr;
                } else if(r == size-1) {
                    edge_type = corner_br;
                } else {
                    edge_type = rightedge;
                }
            } else {
                notedge_h = true;
            }
            if(r == 0) {
                if(edge_type == unset){
                    edge_type = topedge;
                }
            } else if(r == size-1) {
                if(edge_type == unset){
                    edge_type = bottomedge;
                }
            } else {
                notedge_v = true;
            }
            if(notedge_h && notedge_v) {
                edge_type = notedge;
            }
            edge_map[this_crd] = edge_type;
        }
    }
}

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

const int num_grid_cells = 20;
const int num_levels = 10;
const int level_height = 2;

typedef std::map<unsigned int, glm::vec2> TexDims;

struct RoomSettings{
    Shader shader;
    unsigned int fbTexture;
    unsigned int lrTexture;
    unsigned int floorTex;
    unsigned int ceilTex;
    Shader mirror_shader;
    int num_mirrors;
    unsigned int mirrorTex;
    unsigned int innerDoorTex;
    int num_lights;
    unsigned int lightCubeVAO;
    unsigned int lightTex;
    Shader lightShader;
    TexDims tex_dims;
    glm::mat4 view;
    glm::vec3 viewPos;
    glm::mat4 projection;
    float room_scale;
    float sz;
};

struct GridLevelSettings{
    Shader shader;
    unsigned int fbTexture;
    unsigned int lrTexture;
    unsigned int floorTex;
    unsigned int ceilTex;
    Shader mirror_shader;
    int num_mirrors;
    unsigned int mirrorTex;
    unsigned int innerDoorTex;
    TexDims tex_dims;
    float room_scale;
    int height;
    float sz;
};


class GridObj{
public:
    PathTo longest_path;
    Grid grid;
    int level;
    std::unique_ptr<GridLevelSettings> grid_settings;
    GridObj(int size){
        for(int c=0; c < size; ++c){
            std::vector<Cell> col;
            grid.emplace_back(col);
            for(int r=0; r < size; ++r){
                grid[grid.size()-1].push_back(Cell{N,S,E,W});
            }
        }
    }
};

std::vector<std::unique_ptr<GridObj>> all_grids;


unsigned int quadVAO = 0;
unsigned int quadVBO;

void _renderQuad(Shape shape);

// camera
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
float lastX = (float)SCR_WIDTH / 2.0;
float lastY = (float)SCR_HEIGHT / 2.0;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

unsigned int renderRoom(RoomSettings& rsts, int height, int m, int n, glm::vec3 move_to, bool& pass_door_walls, bool pass_doors);

void renderGrids(int level, int grid_map_size, glm::vec3 move_to,
                 glm::vec3 viewPos, glm::mat4 view, glm::mat4 projection,
                 bool& pass_door_walls, bool pass_doors);

unsigned int _renderRoom(Shader& shader, int height, int m, int n, unsigned int lr_dMap,
                         unsigned int fb_dMap, unsigned int c_dMap,
                         unsigned int f_dMap, glm::vec3 move_to,

                         float sz, Shader mirror_shader, int num, unsigned int textureColorbuffer,
                         unsigned int innerDoorTex,

                         int num_lights, unsigned int lightCubeVAO, unsigned int lightTex,
                         Shader& lightShader, TexDims& tex_dims,
                         glm::mat4 view, glm::mat4 projection, glm::vec3 viewPos, float scale, bool& pass_door_walls, bool pass_doors);

unsigned int renderCube(Shader& shader, int height, int m, int n, unsigned int lr_dMap, unsigned int fb_dMap,
                        unsigned int f_dMap, unsigned int c_dMap,
                        TexDims& tex_dims,
                        glm::vec3 move_to, glm::vec3 viewPos, glm::mat4 view, float scale, bool& pass_door_walls);

void renderMirrors(float sz, Shader& shader, int num, int n, int m, unsigned int textureColorbuffer,
                   glm::vec3 move_to,
                   glm::mat4 view, glm::vec3 viewPos,  glm::mat4 projection,
                   bool in_grid, Direction drc, int level);

void renderDoors(float sz, Shader& shader, Shader& normal_shader, int num_doors, int n, int m, int h,
                 unsigned int mirrorTex,
                 unsigned int inner_brickTex, TexDims& tex_dims, bool grid_doors,
                 glm::vec3 move_to, glm::mat4 view, glm::mat4 projection, glm::vec3 viewPos);

void renderLamps(Shader& shader, unsigned int lightCubeVAO, int num_lights, float light_scale, glm::vec3 move_to, float len_m, float wid_n,
                 glm::mat4 view, glm::mat4 projection);

void make_grid(Shader& gridShader, unsigned int grid_VAO, Grid& grid, PathMp& path_map){
    int grid_size = grid[0].size();
    EdgeTMap edge_map;
    getEdges(edge_map, grid_size);
    Coord start{0,0};
    std::set<Coord> added;
    find_next_path(start,path_map,grid,added,edge_map,std::pow(grid_size,2));
    int foo = 3;
    std::cout << path_map.size() << '\n';
    std::cout << foo;
}

void draw_cell(Grid& grid, Coord& crd, Shader& gridShader, unsigned int grid_VAO,
               int grid_size, float scale){
    gridShader.use();
    gridShader.setFloat("cell_size", 0.1f*scale);
    gridShader.setFloat("line_thickness", 0.005f*scale);
    gridShader.setFloat("nr_cells", std::pow(grid_size,2));
    int c = crd.first;
    int r = crd.second;
    gridShader.setVec2("crd", c, r);
    gridShader.setBool("N", grid[c][r].contains(N));
    gridShader.setBool("S", grid[c][r].contains(S));
    gridShader.setBool("E", grid[c][r].contains(E));
    gridShader.setBool("W", grid[c][r].contains(W));
    glBindVertexArray(grid_VAO);
    glDrawArrays(GL_POINTS, 0, 1);
}

//void close_short_paths(){
//
//}

bool not_in_path(Coord& crd, PathTo& path){
    bool valid = true;
    for(auto& p : path.first){
        if(p.first == crd){
            valid = false;
        }
    }
    return valid;
}

void close_n_paths(int level, PathTo& path, Grid& grid){
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<int> gen(0, grid.size()-1);
    for(int i = 0; i < level*2; ++i){
        for(;;){
            bool closed = false;
            Coord crd{gen(g),gen(g)};
            if(not_in_path(crd,path)) {
                std::set<Direction> walls = grid[crd.first][crd.second];
                if(walls.size() >= 3){ break; }
                std::set<Direction> all_walls{N,S,E,W};
                std::vector<Direction> no_walls;
                std::set_difference(all_walls.begin(),all_walls.end(),walls.begin(),walls.end(),
                                    std::back_inserter(no_walls));
                shuffle(no_walls);
                if(no_walls.size() <= 1){
                    break;
                }
                for(Direction drc : no_walls) {
                    switch(drc) {
                        case N : {
                            Coord adj{crd.first,crd.second-1};
                            if(not_in_path(adj,path)){
                                grid[crd.first][crd.second].erase(N);
                                grid[adj.first][adj.second].erase(S);
                                closed = true;
                            }
                            break;
                        }
                        case S : {
                            Coord adj{crd.first,crd.second+1};
                            if(not_in_path(adj,path)){
                                grid[crd.first][crd.second].erase(S);
                                grid[adj.first][adj.second].erase(N);
                                closed = true;
                            }
                            break;
                        }
                        case E : {
                            Coord adj{crd.first+1,crd.second};
                            if(not_in_path(adj,path)){
                                grid[crd.first][crd.second].erase(E);
                                grid[adj.first][adj.second].erase(W);
                                closed = true;
                            }
                            break;
                        }
                        case W : {
                            Coord adj{crd.first-1,crd.second};
                            if(not_in_path(adj,path)){
                                grid[crd.first][crd.second].erase(W);
                                grid[adj.first][adj.second].erase(E);
                                closed = true;
                            }
                            break;
                        }
                    }
                    if(closed)
                         break;
                }
            }
            if(closed)
                 break;
        }
       // std::cout << "c,r" << crd.first << ',' << crd.second << '\n';
    }
}

PathTo find_longest_path(PathMp& paths, Grid& grid){
    int curr_max = 0;
    auto m_iter = paths.begin();
    for(auto iter = paths.begin(); iter != paths.end(); ++iter){
        if(iter->second.second > curr_max){
            Coord end = iter->second.first[iter->second.first.size()-1].first;
            Direction end_drc = iter->second.first[iter->second.first.size()-1].second;
            Direction end_p_drc = iter->second.first[iter->second.first.size()-2].second;
            bool valid = false;
            if(end_drc == end_p_drc){
                if(end_drc == N || end_drc == S){
                    if(grid[end.first][end.second].contains(E) &&
                    grid[end.first][end.second].contains(W) &&
                    grid[end.first][end.second].contains(end_drc)){
                        valid = true;
                    }
                } else {
                    if(grid[end.first][end.second].contains(N) &&
                    grid[end.first][end.second].contains(S) &&
                    grid[end.first][end.second].contains(end_drc)){
                        valid = true;
                    }
                }
            }
            if(valid){
                curr_max = iter->second.second;
                m_iter = iter;
            }
        }
    }
    return m_iter->second;
}

PathTo find_path_to_corner(PathMp& paths){
    int curr_max = 0;
    auto m_iter = paths.begin();
    int max_dist = 0;
    Coord start = paths.begin()->first;
    for(auto iter = paths.begin(); iter != paths.end(); ++iter){
        Coord end = iter->second.first[iter->second.first.size()-1].first;
        int dist = std::pow(std::abs(end.first-start.first),2) +
                std::pow(std::abs(end.second-start.second),2);
        if(dist > max_dist){
            curr_max = iter->second.second;
            m_iter = iter;
        }
    }

    ///Coord prev = m_iter->second.first[0].first;
    ///for(auto p : m_iter->second.first){
    ///    bool lk = true;
    ///    //if((p.first.first != prev.first &&
    ///    //   p.first.first != prev.first+1 &&
    ///    //   p.first.first != prev.first-1) || (
    ///    //   p.first.second != prev.second &&
    ///    //   p.first.second != prev.second+1 &&
    ///    //   p.first.second != prev.second-1)){
    ///    //    int foo = 2;
    ///    //    std::cout << "cell not linked" << '\n';
    ///    //    std::cout << "prev " << prev.first << ',' << prev.second << '\n';
    ///    //    std::cout << "this " << p.first.first << ',' << p.first.second << '\n';
    ///    //} else {
    ///    //    prev = p.first;
    ///    //}
    ///}
    return m_iter->second;
}

void trace_path(Grid& grid, PathTo& longest_path, Shader& gridShader,
                unsigned int grid_VAO, float scale){
    int grid_size = grid[0].size();
    gridShader.use();
    gridShader.setFloat("cell_size", 0.1f*scale);
    gridShader.setFloat("line_thickness", 0.02f*scale);
    gridShader.setFloat("nr_cells", std::pow(grid_size,2));
    gridShader.setBool("trace_path", true);

    //std::vector<std::pair<Coord,Direction>> path = find_longest_path(paths).first;
    std::vector<std::pair<Coord,Direction>> path = longest_path.first;
    int num_n = 0;
    int dups = 0;
    std::set<Coord> v;
    for(auto it = path.begin()+1; it != path.end()-1; ++it){
        int drc_enter = static_cast<int>(opposite_drc[it->second]);
        int drc_exit = static_cast<int>((it+1)->second);
        int drc = (10*std::min(drc_enter,drc_exit))+std::max(drc_enter,drc_exit);
        //if(!std::set{3,2,1,12,13,23}.contains(drc)){
            //std::cout << "DRC " << drc << '\n';
        //}
        if(drc == 1){
            num_n++;
        }
        if(v.contains(it->first)){
            dups++;
        } else {
            v.insert(it->first);
        }

        gridShader.setVec2("crd", it->first.first, it->first.second);
        gridShader.setInt("drc", drc);
        gridShader.setBool("trace_path", true);
        glBindVertexArray(grid_VAO);
        glDrawArrays(GL_POINTS, 0, 1);
    }
    std::cout << "dups " << dups << '\n';
    std::cout << "numn " << num_n << '\n';
    //std::cout << "end crd " << path[path.size()-1].first.first << ',' << path[path.size()-1].first.second << '\n';
}

void draw_grid(Grid& grid, Shader& gridShader, unsigned int grid_VAO, float scale){
    int grid_size = grid[0].size();
    gridShader.use();
    gridShader.setFloat("cell_size", 0.1f*scale);
    gridShader.setFloat("line_thickness", 0.005f*scale);
    gridShader.setFloat("nr_cells", std::pow(grid_size,2));
    gridShader.setBool("trace_path", false);
    for(int c = 0; c < grid_size; ++c){
       for(int r = 0; r < grid_size; ++r){
           gridShader.setVec2("crd", c, r);
           gridShader.setBool("N", grid[c][r].contains(N));
           gridShader.setBool("S", grid[c][r].contains(S));
           gridShader.setBool("E", grid[c][r].contains(E));
           gridShader.setBool("W", grid[c][r].contains(W));
           glBindVertexArray(grid_VAO);
           glDrawArrays(GL_POINTS, 0, 1);
       }
    }
}


int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // lighting info
    // -------------
    glm::vec3 lightPos(20.2f, 0.7f, 2.9f);

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);
    glStencilFunc(GL_NOTEQUAL, 1, 0xFF);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

    glDepthFunc(GL_LEQUAL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //glBlendFunc(GL_ZERO, GL_ONE_MINUS_SRC_ALPHA);

    // build and compile shaders
    // -------------------------

    //    Shader shader("normal_map.vs", "normal_map.fs");
    Shader shader("normal_map_hallw.vs", "normal_map_hallw.fs");
    Shader screenShader("hud.vs", "hud.fs");
    Shader basicShader("framebuffers.vs", "framebuffers.fs");
    //Shader skyboxShader("skybox_hallw.vs", "skybox_hallw.fs");
    Shader lightShader("light_cube.vs", "light_cube.fs");
    Shader mirrorShader("mirror_shader.vs", "mirror_shader.fs");

    Shader gridShader("lines.vs", "lines.fs", "lines.gs");


    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    float grid_points[] = {
            -0.95f,  0.95f, 1.0f, 0.0f, 0.0f, // top-left
    };
    unsigned int grid_VBO, grid_VAO;
    glGenBuffers(1, &grid_VBO);
    glGenVertexArrays(1, &grid_VAO);
    glBindVertexArray(grid_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, grid_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(grid_points), &grid_points, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), 0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);


    for(int i = 0; i < num_levels; ++i){
        all_grids.push_back(std::make_unique<GridObj>(num_grid_cells));
        PathMp path_map;
        make_grid(gridShader, grid_VAO, all_grids[i]->grid, path_map);
        all_grids[i]->longest_path = (find_longest_path(path_map, all_grids[i]->grid));
        close_n_paths(i,all_grids[i]->longest_path,all_grids[i]->grid);
        //close_short_paths(i,all_grids[i]->longest_path,all_grids[i]->grid);
        all_grids[i]->level = i;
    }

    float lightVertices[] = {
            // positions          // normals           // texture coords
            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  0.0f,
            0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  0.0f,
            0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
            0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
            -0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  1.0f,
            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  0.0f,

            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,
            0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  0.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
            -0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  1.0f,
            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,

            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
            -0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            -0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

            0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
            0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
            0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
            0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,
            0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  1.0f,
            0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
            0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
            -0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,

            -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f,
            0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  1.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
            -0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,
            -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f
    };


    float doorVertices[] = {
            9.5f, 0.2f, 2.0f,  0.0f, 0.0f,
            9.5f, 0.2f, 2.0f,  1.0f, 0.0f,
            9.5f, 0.2f, 2.0f,  1.0f, 1.0f,
            9.5f, 0.2f, 2.0f,  1.0f, 1.0f,
            9.5f, 0.2f, 2.0f,  0.0f, 1.0f,
            9.5f, 0.2f, 2.0f,  0.0f, 0.0f
    };

    //float skyboxVertices[] = {
    //        // positions
    //        -500.0f,  500.0f, -500.0f,
    //        -500.0f, -500.0f, -500.0f,
    //         500.0f, -500.0f, -500.0f,
    //         500.0f, -500.0f, -500.0f,
    //         500.0f,  500.0f, -500.0f,
    //        -500.0f,  500.0f, -500.0f,

    //        -500.0f, -500.0f,  500.0f,
    //        -500.0f, -500.0f, -500.0f,
    //        -500.0f,  500.0f, -500.0f,
    //        -500.0f,  500.0f, -500.0f,
    //        -500.0f,  500.0f,  500.0f,
    //        -500.0f, -500.0f,  500.0f,

    //         500.0f, -500.0f, -500.0f,
    //         500.0f, -500.0f,  500.0f,
    //         500.0f,  500.0f,  500.0f,
    //         500.0f,  500.0f,  500.0f,
    //         500.0f,  500.0f, -500.0f,
    //         500.0f, -500.0f, -500.0f,

    //        -500.0f, -500.0f,  500.0f,
    //        -500.0f,  500.0f,  500.0f,
    //         500.0f,  500.0f,  500.0f,
    //         500.0f,  500.0f,  500.0f,
    //         500.0f, -500.0f,  500.0f,
    //        -500.0f, -500.0f,  500.0f,

    //        -500.0f,  500.0f, -500.0f,
    //         500.0f,  500.0f, -500.0f,
    //         500.0f,  500.0f,  500.0f,
    //         500.0f,  500.0f,  500.0f,
    //        -500.0f,  500.0f,  500.0f,
    //        -500.0f,  500.0f, -500.0f,

    //        -500.0f, -500.0f, -500.0f,
    //        -500.0f, -500.0f,  500.0f,
    //         500.0f, -500.0f, -500.0f,
    //         500.0f, -500.0f, -500.0f,
    //        -500.0f, -500.0f,  500.0f,
    //         500.0f, -500.0f,  500.0f
    //};



    // shader configuration
    // --------------------
    screenShader.use();
    screenShader.setInt("screenTexture", 0);

    // framebuffer configuration
    // -------------------------
    unsigned int framebuffer;
    glGenFramebuffers(1, &framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
    // create a color attachment texture
    unsigned int textureColorbuffer;
    glGenTextures(1, &textureColorbuffer);
    glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);
    // create a renderbuffer object for depth and stencil attachment (we won't be sampling these)
    // unsigned int rbo;
    // glGenRenderbuffers(1, &rbo);
    // glBindRenderbuffer(GL_RENDERBUFFER, rbo);
    // glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, SCR_WIDTH, SCR_HEIGHT); // use a single renderbuffer object for both a depth AND stencil buffer.
    // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo); // now actually attach it
    // now that we actually created the framebuffer and added all attachments we want to check if it is actually complete now
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << endl;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    shader.use();
    shader.setInt("brickTex", 0);

    unsigned int VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    unsigned int lightCubeVBO;
    glGenBuffers(1, &lightCubeVBO);

    glBindBuffer(GL_ARRAY_BUFFER, lightCubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(lightVertices), lightVertices, GL_STATIC_DRAW);

    // second, configure the light's VAO (VBO stays the same; the vertices are the same for the light object which is also a 3D cube)
    unsigned int lightCubeVAO;
    glGenVertexArrays(1, &lightCubeVAO);
    glBindVertexArray(lightCubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, lightCubeVBO);
    // note that we update the lamp's position attribute's stride to reflect the updated buffer data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // skybox VAO
    //unsigned int skyboxVAO, skyboxVBO;
    //glGenVertexArrays(1, &skyboxVAO);
    //glGenBuffers(1, &skyboxVBO);
    //glBindVertexArray(skyboxVAO);
    //glBindBuffer(GL_ARRAY_BUFFER, skyboxVBO);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(skyboxVertices), &skyboxVertices, GL_STATIC_DRAW);
    //glEnableVertexAttribArray(0);
    //glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

    //vector<std::string> faces
    //{
    //    FileSystem::getPath("resources/textures/skybox2/negx.jpg"),
    //    FileSystem::getPath("resources/textures/skybox2/posx.jpg"),
    //    FileSystem::getPath("resources/textures/skybox2/posy.jpg"),
    //    FileSystem::getPath("resources/textures/skybox2/negy.jpg"),
    //    FileSystem::getPath("resources/textures/skybox2/negz.jpg"),
    //    FileSystem::getPath("resources/textures/skybox2/posz.jpg")
    //};
    //unsigned int cubemapTexture = loadCubemap(faces);
    //skyboxShader.use();
    //skyboxShader.setInt("skybox", 0);

    unsigned int lastQuad = 0;
    unsigned int total_rooms = 4;

    glm::vec3 origin(0.0f,0.0f,0.0f);
    glm::vec3 move1(-20.0f,0.0f,10.0f);
    glm::vec3 move2(3.05f,0.0f,0.0f);
    glm::vec3 move3(3.0f,0.0f,-30.0f);
    glm::vec3 move4(3.05f,0.0f,-60.0f);
    glm::vec3 move5(20.0f,0.0f,-63.0f);

    camera.AddRoom(2,80,3, origin);
    camera.AddRoom(2,20,20,move1);
    camera.AddRoom(3,15,30,move2);
    camera.AddRoom(2,10,20,move3);
    camera.AddRoom(2,3,70,move4);
    camera.AddRoom(2,30,3,move5);

    camera.Position = glm::vec3(1.0f, 1.0f, 1.0f);

    int num_lights = 8;
    int num_mirrors = 20;
    float room_scale = 1.0f;

    // load textures
    // -------------
    std::map<unsigned int, glm::vec2> texture_dims; //map of texture IDs to width/height pairs
    unsigned int purpleglassTex  = loadTexture(FileSystem::getPath("resources/textures/window2_nmp.png").c_str());
    get_texture_dimensions("resources/textures/window2.png", purpleglassTex, texture_dims);
    unsigned int mirrorTex  = loadTexture(FileSystem::getPath("resources/textures/glass3.png").c_str());
    get_texture_dimensions("resources/textures/glass3.png", mirrorTex, texture_dims);
    unsigned int frameWindowTex  = loadTexture(FileSystem::getPath("resources/textures/window.png").c_str());
    get_texture_dimensions("resources/textures/window.png", frameWindowTex, texture_dims);
    unsigned int crystalTex  = loadTexture(FileSystem::getPath("resources/textures/crystal.png").c_str());
    get_texture_dimensions("resources/textures/crystal.png", crystalTex, texture_dims);

    for(int i = 0; i < num_levels; ++i){
        unsigned int wallTex = loadTexture(FileSystem::getPath("resources/textures/wall"+std::to_string(i)+".png").c_str());
        unsigned int floorTex = loadTexture(FileSystem::getPath("resources/textures/floor"+std::to_string(i)+".png").c_str());

        all_grids[i]->grid_settings = std::make_unique<GridLevelSettings>(
                GridLevelSettings{shader,wallTex,wallTex,floorTex,wallTex,
                                  basicShader, num_mirrors, textureColorbuffer, purpleglassTex,
                                  texture_dims, room_scale, level_height, 0.3});

        get_texture_dimensions("resources/textures/wall"+std::to_string(i)+".png", wallTex, texture_dims);
        get_texture_dimensions("resources/textures/wall"+std::to_string(i)+".png", floorTex, texture_dims);
    }


    auto time_start = std::chrono::system_clock::now();
    int current_level = 0;

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        //processInput(window, longest_paths[0], all_grids[0]->grid);
        processInput(window);

        // render
        // ------
        // bind to framebuffer and draw scene as we normally would to color texture
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
        //glEnable(GL_DEPTH_TEST); // enable depth testing (is disabled for rendering screen-space quad)

        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        camera.Yaw   += 180.0f; // rotate the camera's yaw 180 degrees around
        camera.ProcessMouseMovement(0, 0, false); // call this to make sure it updates its camera vectors, note that we disable pitch constrains for this specific case (otherwise we can't reverse camera's pitch values)
        glm::mat4 view = camera.GetViewMatrix();
        camera.Yaw   -= 180.0f; // reset it back to its original orientation
        camera.ProcessMouseMovement(0, 0, true);
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom-10), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 1000.0f);


        shader.use();
        shader.setMat4("projection", projection);
        shader.setMat4("view", view);
        auto time_end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = time_end - time_start;
        std::chrono::milliseconds integral_diff =
                std::chrono::duration_cast<std::chrono::milliseconds>(diff);
        shader.setInt("utime", integral_diff.count());
        // render normal-mapped quad
        glm::mat4 model = glm::mat4(1.0f);

        glm::vec3 origin(0.0f,0.0f,0.0f);

        glm::vec3 room_left(20.0f, 0.0f, -40.0f);

        glm::vec3 viewPos = camera.Position;

        bool pass_door_walls = false;



        unsigned int bTex = all_grids[0]->grid_settings->lrTexture;
        unsigned int fTex = all_grids[0]->grid_settings->floorTex;
        RoomSettings rsts{shader,bTex,bTex,fTex,bTex,
                          basicShader, num_mirrors, textureColorbuffer, purpleglassTex,
                          num_lights, lightCubeVAO, crystalTex, lightShader, texture_dims, view, viewPos, projection, room_scale, 0.3};

        //std::vector<Grid*> gsts_grids = std::vector<Grid*>{&all_grids[0]->grid};
        //GridLevelSettings gsts{shader,brickTex,brickTex,floorTex,brickTex,
        //                  basicShader, num_mirrors, textureColorbuffer, purpleglassTex,
        //                  texture_dims, view, viewPos, projection, room_scale, 0.3};

        int ww, wh;
        glfwGetWindowSize(window, &ww, &wh);
        shader.setVec2("scrRes", ww, wh);



        for(auto i = 0; i < 1; ++i){
            renderGrids(i, 300,origin,
                        camera.Position, view, projection,
                        pass_door_walls,false);
            renderGrids(i, 300,origin,
                        camera.Position, view, projection,
                        pass_door_walls,true);
        }
        /////lastQuad += renderRoom(rsts,2,80,3,origin,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts,2,20,20,move1,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 3,15,30,move2,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,10,20,move3,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,3,70,move4,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,30,3,move5,pass_door_walls, false);

        /////lastQuad += renderRoom(rsts,2,80,3,origin,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts,2,20,20,move1,pass_door_walls, true);
        ///////glClear(GL_STENCIL_BUFFER_BIT);
        /////lastQuad += renderRoom(rsts, 3,15,30,move2,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts, 2,10,20,move3,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts, 2,3,70,move4,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts, 2,30,3,move5,pass_door_walls, true);


        pass_door_walls = true;


        /////lastQuad += renderRoom(rsts,2,80,3,origin,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts,2,20,20,move1,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 3,15,30,move2,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,10,20,move3,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,3,70,move4,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,30,3,move5,pass_door_walls, false);

        // second render pass: draw as normal
        // ----------------------------------
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDepthFunc(GL_LEQUAL);
        view = camera.GetViewMatrix();

        shader.use();
        shader.setMat4("projection", projection);
        shader.setMat4("view", view);

        rsts.viewPos = camera.Position;
        rsts.view = view;
        rsts.projection = projection;

        pass_door_walls = false;

        for(auto i = 0; i < 1; ++i){
            renderGrids(i, 300, origin,
                        camera.Position, view, projection,
                        pass_door_walls,false);
            renderGrids(i, 300,origin,
                        camera.Position, view, projection,
                        pass_door_walls,true);
        }

        glm::vec3 move_to_r = origin;
       // move_to_r.z += 30.0f;
       // lastQuad += renderRoom(rsts,2,80,3,origin,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts,2,20,20,move1,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 3,15,30,move2,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,10,20,move3,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,3,70,move4,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,30,3,move5,pass_door_walls, false);

        /////lastQuad += renderRoom(rsts,2,80,3,origin,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts,2,20,20,move1,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts, 3,15,30,move2,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts, 2,10,20,move3,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts, 2,3,70,move4,pass_door_walls, true);
        /////lastQuad += renderRoom(rsts, 2,30,3,move5,pass_door_walls, true);


        pass_door_walls = true;

        /////lastQuad += renderRoom(rsts,2,80,3,origin,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts,2,20,20,move1,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 3,15,30,move2,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,10,20,move3,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,3,70,move4,pass_door_walls, false);
        /////lastQuad += renderRoom(rsts, 2,30,3,move5,pass_door_walls, false);


        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, purpleglassTex);	// use the color attachment texture as the texture of the quad plane
        glm::mat4 view_2 = lookAt(glm::vec3(0,0,1),glm::vec3(0,0,0),glm::vec3(0,1,0));

        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        glClear(GL_STENCIL_BUFFER_BIT);
        glDepthFunc(GL_ALWAYS);
        screenShader.use();
        screenShader.setMat4("view", view_2);
        float hud_alpha = 0.15f;


        screenShader.setVec2("scrRes", ww, wh);
     //   std::cout << "ww: " << ww << '\n';
      //  std::cout << "wh: " << wh << '\n';

        _renderQuad(Quad);
        screenShader.setVec2("move_to", glm::vec2(-1.5f,1.5f));
        glm::mat4 model0 = glm::mat4(1.0f);
        model0 = glm::scale(model0, glm::vec3(0.5f));
        screenShader.setMat4("model", model0);
        screenShader.setFloat("alpha", hud_alpha);
        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);


        //screenShader.setVec2("move_to", glm::vec2(1.2f,1.2f));
        //glm::mat4 model_cnr = glm::mat4(1.0f);
        //model_cnr = glm::scale(model_cnr, glm::vec3(0.3f));
        //model_cnr = glm::rotate(model_cnr, glm::radians(45.0f), glm::vec3(0.0f, 0.0f, 1.0f));
        //screenShader.setMat4("model", model_cnr);
        //glBindVertexArray(quadVAO);
        //glDrawArrays(GL_TRIANGLES, 0, 6);

        glStencilMask(0xFF);
        glStencilFunc(GL_ALWAYS, 1, 0xFF);

        _renderQuad(Triangle);
        screenShader.setVec2("move_to", glm::vec2(3.0f,0.5f));
        glm::mat4 model3 = glm::mat4(1.0f);
        model3 = glm::scale(model3, glm::vec3(0.2f));
        model3 = glm::rotate(model3, glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
        screenShader.setMat4("model", model3);
        screenShader.setFloat("alpha", 0.0f);
        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glStencilFunc(GL_NOTEQUAL, 1, 0xFF);
        glStencilMask(0x00);

        _renderQuad(Quad);
        screenShader.setVec2("move_to", glm::vec2(-0.5f,1.5f));
        glm::mat4 model1 = glm::mat4(1.0f);
        model1 = glm::scale(model1, glm::vec3(0.5f));
        screenShader.setMat4("model", model1);
        screenShader.setFloat("alpha", hud_alpha);
        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);

        glBindVertexArray(0);
        glStencilMask(0xFF);
        glStencilFunc(GL_ALWAYS, 0, 0xFF);
        glEnable(GL_DEPTH_TEST);

        //glClear(GL__BUFFER_BIT);
        glDepthFunc(GL_ALWAYS);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        float grid_scale = 0.2f;
        draw_grid(all_grids[current_level]->grid, gridShader, grid_VAO, grid_scale);
        trace_path(all_grids[current_level]->grid, all_grids[current_level]->longest_path,
                   gridShader,grid_VAO,grid_scale);
       // Coord o = Coord{0,0};
      //  std::cout << grid[0][0].contains(N) << ' ' << grid[0][0].contains(S) << ' '
      //  << grid[0][0].contains(E) << ' ' << grid[0][0].contains(W) << '\n';
      //  draw_cell(grid, o , gridShader, grid_VAO, grid[0].size());

        //// draw skybox as last
        //glDepthFunc(GL_LEQUAL);  // change depth function so depth test passes when values are equal to depth buffer's content
        //skyboxShader.use();
        //view = glm::mat4(glm::mat3(camera.GetViewMatrix())); // remove translation from the view matrix
        //skyboxShader.setMat4("view", view);
        //skyboxShader.setMat4("projection", projection);
        //// skybox cube
        //glBindVertexArray(skyboxVAO);
        //glActiveTexture(GL_TEXTURE0);
        //glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
        //glDrawArrays(GL_TRIANGLES, 0, 36);
        //glBindVertexArray(0);
        //glDepthFunc(GL_LESS); // set depth function back to default
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &lightCubeVAO);
    // glDeleteVertexArrays(1, &screenquadVAO);
    glDeleteBuffers(1, &VBO);
    //glDeleteBuffers(1, &screenquadVBO);
    glfwTerminate();
    return 0;
}

void renderLamps(Shader& shader, unsigned int lightCubeVAO, int num_lights, float light_scale, glm::vec3 move_to, float len_m, float wid_n,
                 glm::mat4 view, glm::mat4 projection){
    // draw the lamp object
    float il = len_m / 4;
    float il2 = il/2;
    shader.use();
    for(int k = 0; k < num_lights/2; ++k){
        float il_pos = (il*k)+il2;
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::scale(model, glm::vec3(light_scale)); // a smaller cube
        shader.setMat4("projection", projection);
        shader.setMat4("view", view);
        model = glm::translate(model, glm::vec3(-5*move_to.z, move_to.y, 5*move_to.x));

        model = glm::translate(model, glm::vec3( il_pos*(1/light_scale),  1.6f, 0.3f));
        shader.setMat4("model", model);
        //shader.setBool("isLightSrc", true);
        glBindVertexArray(lightCubeVAO);
        glDrawArrays(GL_TRIANGLES, 0, 36);

        glm::mat4 model_2 = glm::mat4(1.0f);
        model_2 = glm::scale(model_2, glm::vec3(light_scale)); // a smaller cube
        //lightShader.setMat4("projection", projection);
        //lightShader.setMat4("view", view);
        model_2 = glm::translate(model_2, glm::vec3(-5*move_to.z, move_to.y, 5*move_to.x));
        model_2 = glm::translate(model_2, glm::vec3( il_pos*(1/light_scale),  1.6f, (1/light_scale)*static_cast<float>(wid_n-0.01)));
        shader.setMat4("model", model_2);
        glDrawArrays(GL_TRIANGLES, 0, 36);
        //if(!add_to_dcs){
        //    glDrawArrays(GL_TRIANGLES, 0, 36);
        //} else {
        //    DrawCall dc(36, shader, lightCubeVAO, 0, 0, false, false, model_2);
        //    dcs.push_back(dc);
        //}
    }
}

void initShader(Shader& shader, int num_lights, glm::mat4 view, bool is_grid, int m, Grid& grid){
    shader.use();
    shader.setVec3("viewPos", camera.Position);
    shader.setMat4("view", view);
    unsigned int lastQuad = 0;

    shader.setVec3("dirLight.direction", -0.2f, -1.0f, -0.3f);
    shader.setVec3("dirLight.ambient", 0.05f, 0.05f, 0.05f);
    shader.setVec3("dirLight.diffuse", 0.4f, 0.4f, 0.4f);
    shader.setVec3("dirLight.specular", 0.5f, 0.5f, 0.5f);

    shader.setInt("nr_point_lights", num_lights);
    for(int i = 0; i < num_lights; ++i){
        shader.setVec3("pointLights["+std::to_string(i)+"].ambient", 0.05f, 0.05f, 0.05f);
        shader.setVec3("pointLights["+std::to_string(i)+"].diffuse", 0.8f, 0.8f, 0.8f);
        shader.setVec3("pointLights["+std::to_string(i)+"].specular", 1.0f, 1.0f, 1.0f);
        shader.setFloat("pointLights["+std::to_string(i)+"].constant", 1.0f);
        shader.setFloat("pointLights["+std::to_string(i)+"].linear", 0.09f);
        shader.setFloat("pointLights["+std::to_string(i)+"].quadratic", 0.032f);
    }
    // spotLight
    shader.setVec3("spotLight.position", camera.Position);
    shader.setVec3("spotLight.direction", camera.Front);
    shader.setVec3("spotLight.ambient", 0.3f, 0.3f, 0.3f);
    shader.setVec3("spotLight.diffuse", 1.0f, 1.0f, 1.0f);
    shader.setVec3("spotLight.specular", 1.0f, 1.0f, 1.0f);
    shader.setFloat("spotLight.constant", 1.0f);
    shader.setFloat("spotLight.linear", 0.09f);
    shader.setFloat("spotLight.quadratic", 0.032f);
    shader.setFloat("spotLight.cutOff", glm::cos(glm::radians(12.5f)));
    shader.setFloat("spotLight.outerCutOff", glm::cos(glm::radians(15.0f)));

    // material properties
    shader.setFloat("shininess", 52.0f);
    float il = m / (num_lights/2);
    float il2 = il/2;

    if(!is_grid){
        for(int i = 0; i < num_lights/2; ++i){
            float il_pos = (il*i)+il2;
            shader.setVec3("pointLights["+std::to_string(i)+"].position", il_pos,  0.3f, 0.3f);
        }
        for(int i = 0; i < num_lights/2; ++i){
            float il_pos = (il*i)+il2;
            shader.setVec3("pointLights["+std::to_string((num_lights/2)+i)+"].position", il_pos,  0.3f, 2.9f);
        }
    } else {

    }
}

void renderGrids(int level, int grid_map_size,
                 glm::vec3 move_to,
                 glm::vec3 viewPos, glm::mat4 view, glm::mat4 projection,
                 bool& pass_door_walls, bool pass_doors){

     initShader(all_grids[level]->grid_settings->shader,2,view, true, std::floor(grid_map_size/num_grid_cells), all_grids[level]->grid);

     unsigned int lastQuad = 0;
     float quad_sz = 1.5;
     glm::vec3 lightPos(10.5f, 1.0f, 0.3f);
     glm::mat4 model = glm::mat4(1.0f);
     auto gsts = all_grids[level]->grid_settings.get();
     gsts->shader.setMat4("model", model);
     gsts->shader.setMat4("view", view);
     gsts->shader.setVec3("viewPos", viewPos);
     gsts->shader.setVec3("lightPos", lightPos);
     int grid_sz = all_grids[level]->grid.size();

     float wall_scale = 0.155;
     int m = 15; int n = 1;
     int wall_height = 14;

     glm::vec3 level_move_to = move_to;
     level_move_to.y += (level*(gsts->height+0.01));
     glm::vec3 move_to_walls = level_move_to;
     move_to_walls.z += 2.7*wall_scale;
     move_to_walls.y -= 2.8*wall_scale;

     int num_paths = all_grids[level]->longest_path.first.size();
     Coord start = all_grids[level]->longest_path.first[0].first;
     Coord end = all_grids[level]->longest_path.first[num_paths-1].first;
     for(int r = 0; r < grid_sz; ++r){
         for(int c = 0; c < grid_sz; ++c){

             if(pass_doors){
                 glm::vec3 move_to_d = move_to_walls;
                 move_to_d.x += (r+0.5)*m*wall_scale;
                 move_to_d.z -= (c+0.25)*m*wall_scale;
                 if(c == start.first && r == start.second){
                     renderDoors(10, gsts->mirror_shader, gsts->shader, 9, n, m, gsts->height, gsts->mirrorTex, gsts->innerDoorTex,
                                 gsts->tex_dims, true,
                                 move_to_d, view, projection, viewPos);
                 }
                 if(c == end.first && r == end.second){
                     renderDoors(10, gsts->mirror_shader, gsts->shader, 9, n, m, gsts->height, gsts->mirrorTex, gsts->innerDoorTex,
                                 gsts->tex_dims, true,
                                 move_to_d, view, projection, viewPos);
                 }
             } else {
                 if(all_grids[level]->grid[c][r].contains(N)){
                     //north
                     glm::vec3 move_to_n = move_to_walls;
                     move_to_n.x += r*m*wall_scale;
                     move_to_n.z -= c*m*wall_scale;
                     bool pd_walls = false;
                     renderCube(gsts->shader, wall_height, m+1, n, gsts->lrTexture, gsts->fbTexture, gsts->floorTex, gsts->ceilTex,
                                gsts->tex_dims, move_to_n, viewPos, view, wall_scale, pd_walls);

                     renderMirrors(0.6, gsts->mirror_shader, 1, 15, 20, gsts->mirrorTex, move_to_n, view, camera.Position, projection,true, N, level);
                     pd_walls = true;
                     gsts->shader.use();
                     renderCube(gsts->shader, wall_height, m+1, n, gsts->lrTexture, gsts->fbTexture, gsts->floorTex, gsts->ceilTex,
                                gsts->tex_dims, move_to_n, viewPos, view, wall_scale, pd_walls);

                 }

                 if(all_grids[level]->grid[c][r].contains(W)){
                     //west
                     glm::vec3 move_to_w = move_to_walls;
                     move_to_w.x += r*m*wall_scale;
                     move_to_w.z -= c*m*wall_scale;
                     bool pd_walls = false;
                     renderCube(gsts->shader, wall_height, n, m, gsts->lrTexture, gsts->fbTexture, gsts->lrTexture, gsts->ceilTex,
                                gsts->tex_dims, move_to_walls, viewPos, view, wall_scale, pd_walls);
                     pd_walls = true;
                     renderCube(gsts->shader, wall_height, n, m, gsts->lrTexture, gsts->fbTexture, gsts->floorTex, gsts->ceilTex,
                                gsts->tex_dims, move_to_walls, viewPos, view, wall_scale, pd_walls);
                     //renderMirrors(0.6, mirror_shader, 1, n, m, textureColorbuffer, move_to_w, view, camera.Position, projection, true, W);
                 }

                 if(all_grids[level]->grid[c][r].contains(E)){
                     //east
                     glm::vec3 move_to_e = move_to_walls;
                     move_to_e.x += r*m*wall_scale;
                     move_to_e.z -= (c+1)*m*wall_scale;
                     bool pd_walls = false;
                     renderCube(gsts->shader, wall_height, n, m, gsts->lrTexture, gsts->fbTexture, gsts->lrTexture, gsts->ceilTex,
                                gsts->tex_dims, move_to_e, viewPos, view, wall_scale, pd_walls);
                     pd_walls = true;
                     renderCube(gsts->shader, wall_height, n, m, gsts->lrTexture, gsts->fbTexture, gsts->lrTexture, gsts->ceilTex,
                                gsts->tex_dims, move_to_e, viewPos, view, wall_scale, pd_walls);
                     //renderMirrors(0.6, mirror_shader, 1, n, m, textureColorbuffer, move_to_e, view, camera.Position, projection, true, E);
                 }

                 if(all_grids[level]->grid[c][r].contains(S)){
                     //south
                     glm::vec3 move_to_s = move_to_walls;
                     move_to_s.x += (r+1)*m*wall_scale;
                     move_to_s.z -= c*m*wall_scale;
                     bool pd_walls = false;
                     renderCube(gsts->shader, wall_height, m+1, n, gsts->lrTexture, gsts->fbTexture, gsts->lrTexture, gsts->ceilTex,
                                gsts->tex_dims, move_to_s, viewPos, view, wall_scale, pd_walls);
                     renderMirrors(0.6, gsts->mirror_shader, 1, 15, 20, gsts->mirrorTex, move_to_s, view, camera.Position, projection,true, N, level);
                     gsts->shader.use();
                     pd_walls = true;
                     renderCube(gsts->shader, wall_height, m+1, n, gsts->lrTexture, gsts->fbTexture, gsts->lrTexture, gsts->ceilTex,
                                gsts->tex_dims, move_to_s, viewPos, view, wall_scale, pd_walls);
                     //renderMirrors(0.6, mirror_shader, 1, n, m, textureColorbuffer, move_to_s, view, camera.Position, projection, true, S);
                 }
             }
         }
     }


     if(!pass_door_walls && !pass_doors) {

         glActiveTexture(GL_TEXTURE0);
         glBindTexture(GL_TEXTURE_2D, gsts->ceilTex);
         //glActiveTexture(GL_TEXTURE1);
         //glBindTexture(GL_TEXTURE_2D, c_nMap);
         //gsts->shader.setVec2("texRes", gsts->tex_dims[c_dMap]);

         int map_sz = std::ceil(wall_scale*grid_map_size);
         for(int i = 0; i < map_sz; ++i){
             for(int j = 0; j < map_sz; ++j){
                 glm::mat4 model_cl = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
                 model_cl = glm::rotate(model_cl, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                 model_cl = glm::translate(model_cl, glm::vec3(0.0f,gsts->room_scale*0.5f,-gsts->room_scale*0.5));
                 model_cl = glm::translate(model_cl, glm::vec3(gsts->room_scale*static_cast<float>(i),gsts->room_scale*static_cast<float>(j),-static_cast<float>((gsts->room_scale*gsts->height)-1)));
                 model_cl = glm::translate(model_cl, glm::vec3(-level_move_to.z, level_move_to.x, -level_move_to.y));
                 model_cl = glm::scale(model_cl, glm::vec3(gsts->room_scale));
                 gsts->shader.setMat4("model", model_cl);
                 lastQuad = renderQuad();
             }
         }

         glActiveTexture(GL_TEXTURE0);
         glBindTexture(GL_TEXTURE_2D, gsts->floorTex);
         //glActiveTexture(GL_TEXTURE1);
         //glBindTexture(GL_TEXTURE_2D, f_nMap);
         //shader.setVec2("texRes", tex_dims[f_dMap]);

         for(int i = 0; i < map_sz; ++i){
             for(int j = 0; j < map_sz; ++j){
                 glm::mat4 model_fl = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
                 model_fl = glm::rotate(model_fl, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                 model_fl = glm::rotate(model_fl, glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                 model_fl = glm::translate(model_fl, glm::vec3(0.0f,gsts->room_scale*0.5f,-gsts->room_scale*0.5));
                 model_fl = glm::translate(model_fl, glm::vec3(-gsts->room_scale*static_cast<float>(i),gsts->room_scale*static_cast<float>(j),0.0f));
                 model_fl = glm::translate(model_fl, glm::vec3(level_move_to.z, level_move_to.x, level_move_to.y));
                 model_fl = glm::scale(model_fl, glm::vec3(gsts->room_scale));
                 gsts->shader.setMat4("model", model_fl);
                 lastQuad = renderQuad();
             }
         }
     }

}

unsigned int renderRoom(RoomSettings& rsts, int height, int m, int n, glm::vec3 move_to, bool& pass_door_walls, bool pass_doors){

    return _renderRoom(rsts.shader, height, m, n, rsts.fbTexture, rsts.lrTexture,
                       rsts.floorTex, rsts.ceilTex, move_to,
                       rsts.sz, rsts.mirror_shader, rsts.num_mirrors, rsts.mirrorTex, rsts.innerDoorTex,
                       rsts.num_lights, rsts.lightCubeVAO, rsts.lightTex,
                       rsts.lightShader, rsts.tex_dims,
                       rsts.view, rsts.projection, rsts.viewPos, rsts.room_scale, pass_door_walls, pass_doors);
}

//unsigned int mirrorVAO()
unsigned int _renderRoom(Shader& shader, int height, int m, int n, unsigned int lr_dMap, unsigned int fb_dMap,
                         unsigned int f_dMap, unsigned int c_dMap, glm::vec3 move_to,
                         float sz, Shader mirror_shader, int num_mirrors, unsigned int textureColorbuffer,
                         unsigned int innerDoorTex,
                         int num_lights, unsigned int lightCubeVAO, unsigned int lightTex,
                         Shader& lightShader, TexDims& tex_dims,
                         glm::mat4 view, glm::mat4 projection, glm::vec3 viewPos, float scale, bool& pass_door_walls, bool pass_doors){


    initShader(shader,2,view, true, -1, all_grids[0]->grid);
    unsigned int lastQuad = 0;
    float light_scale = 0.1f;

    if(pass_door_walls){
        glStencilFunc(GL_NOTEQUAL, 1, 0xFF);
        glStencilMask(0x00);
        lastQuad = renderCube(shader, height, m, n, lr_dMap, fb_dMap, f_dMap, c_dMap, tex_dims, move_to, viewPos, view, scale, pass_door_walls);
        glBindVertexArray(0);
        glStencilMask(0xFF);
        glStencilFunc(GL_ALWAYS, 0, 0xFF);
        glEnable(GL_DEPTH_TEST);

    } else {

        if(pass_doors){
            //glClear(GL_STENCIL_BUFFER_BIT);
            renderDoors(10, mirror_shader, shader, 9, n, m, height, textureColorbuffer, innerDoorTex,
                        tex_dims, false,
                        move_to, view, projection, viewPos);
        } else {
            lastQuad = renderCube(shader, height, m, n, lr_dMap, fb_dMap, f_dMap, c_dMap, tex_dims, move_to, viewPos, view, scale, pass_door_walls);

            // draw the lamp object
            renderLamps(shader, lightCubeVAO, num_lights, light_scale, move_to, m, n,
                        view, projection);

            renderMirrors(0.6, mirror_shader, num_mirrors, n, m, textureColorbuffer, move_to, view, camera.Position, projection, false, N, -1);

        }

    }
    return lastQuad != 0 ? 1 : 0;
}

unsigned int cubeVAO = 0;
unsigned int cubeVBO;
unsigned int renderCubeLayer(Shader& shader, int layer, unsigned int fbTexture, unsigned int lrTexture,
                             int m, int n, TexDims& tex_dims, glm::vec3 move_to, glm::vec3 viewPos, glm::mat4 view, float scale, bool& pass_door_walls)
{
    unsigned int lastQuad = 0;
    glm::vec3 lightPos(10.5f, 1.0f, 0.3f);        /////////////////
    glm::mat4 model = glm::mat4(1.0f);
    shader.setVec3("viewPos", viewPos);
    shader.setMat4("view", view);
    float quad_sz = 1.5;
    //model = glm::rotate(model, glm::radians((float)glfwGetTime() * -10.0f), glm::normalize(glm::vec3(1.0, 0.0, 1.0))); // rotate the quad to show normal mapping from multiple directions

    if(pass_door_walls){
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, lrTexture);
        //glActiveTexture(GL_TEXTURE1);
        shader.setVec2("texRes", tex_dims[lrTexture]);


        glm::mat4 model_f = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
        glm::mat4 model_b = glm::rotate(model, glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        model_b = glm::translate(model_b, glm::vec3(0.0f,0.0f,-scale*static_cast<float>(n)));

        model_f = glm::translate(model_f, glm::vec3(-move_to.z,move_to.y, move_to.x));
        model_b = glm::translate(model_b, glm::vec3(move_to.z,move_to.y, -move_to.x));

        for(int i = 0; i < m; ++i){
            glm::mat4 model_ff = glm::translate(model_f, glm::vec3(1.0f*i*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_ff, glm::vec3(scale)));
            lastQuad = renderQuad();

            glm::mat4 model_bb = glm::translate(model_b, glm::vec3(-1.0f*i*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_bb, glm::vec3(scale)));
            lastQuad = renderQuad();
        }
    } else {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fbTexture);
        //glActiveTexture(GL_TEXTURE1);
        shader.setVec2("texRes", tex_dims[fbTexture]);

        glm::mat4 model_l = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
        model_l = glm::rotate(model_l, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 model_r = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
        model_r = glm::rotate(model_r, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        model_r = glm::translate(model_r, glm::vec3(0.0f,0.0f,-scale*static_cast<float>(m)));
        model_l = glm::translate(model_l, glm::vec3(-scale*0.5f,scale*0.0f,-scale*0.5f));
        model_r = glm::translate(model_r, glm::vec3(scale*0.5f,scale*0.0f,scale*0.5f));

        model_l = glm::translate(model_l, glm::vec3(-move_to.x,move_to.y,-move_to.z));
        model_r = glm::translate(model_r, glm::vec3(move_to.x,move_to.y,move_to.z));

        for(int j = 0; j < n; ++j){

            glm::mat4 model_ll = glm::translate(model_l, glm::vec3(-1.0f*j*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_ll, glm::vec3(scale)));
            lastQuad = renderQuad();

            glm::mat4 model_rr = glm::translate(model_r, glm::vec3(1.0f*j*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_rr, glm::vec3(scale)));
            lastQuad = renderQuad();
        }
        //pass_door_walls = true;
    }

    return lastQuad;
}

unsigned int renderCubeDirection(Shader& shader, int layer, unsigned int fbTexture, unsigned int lrTexture,
                             int len_f_to_b, int len_l_to_r, TexDims& tex_dims, glm::vec3 move_to, glm::vec3 viewPos, glm::mat4 view, float scale, Direction drc)
{
    unsigned int lastQuad = 0;
    glm::vec3 lightPos(10.5f, 1.0f, 0.3f);        /////////////////
    glm::mat4 model = glm::mat4(1.0f);
    shader.setVec3("viewPos", viewPos);
    shader.setMat4("view", view);
    float quad_sz = 1.5;
    //model = glm::rotate(model, glm::radians((float)glfwGetTime() * -10.0f), glm::normalize(glm::vec3(1.0, 0.0, 1.0))); // rotate the quad to show normal mapping from multiple directions

    if(drc == N || drc == S){
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, lrTexture);
        shader.setVec2("texRes", tex_dims[lrTexture]);

        glm::mat4 model_l = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
        glm::mat4 model_r = glm::rotate(model, glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        model_r = glm::translate(model_r, glm::vec3(0.0f,0.0f,-scale*static_cast<float>(len_l_to_r)));

        model_l = glm::translate(model_l, glm::vec3(-move_to.z,move_to.y, move_to.x));
        model_r = glm::translate(model_r, glm::vec3(move_to.z,move_to.y, -move_to.x));

        if(drc == N){

        } else if(drc == S){

        }

        for(int i = 0; i < len_f_to_b; ++i){
            glm::mat4 model_ll = glm::translate(model_l, glm::vec3(1.0f*i*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_ll, glm::vec3(scale)));
            lastQuad = renderQuad();

            glm::mat4 model_rr = glm::translate(model_r, glm::vec3(-1.0f*i*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_rr, glm::vec3(scale)));
            lastQuad = renderQuad();
        }
    } else {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, lrTexture);
        shader.setVec2("texRes", tex_dims[lrTexture]);
        if(drc == W){
            glm::mat4 model_r = glm::rotate(model, glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
            model_r = glm::translate(model_r, glm::vec3(0.0f,0.0f,-scale*static_cast<float>(len_l_to_r)));
            model_r = glm::translate(model_r, glm::vec3(move_to.z,move_to.y, -move_to.x));
            for(int i = 0; i < len_f_to_b; ++i){
                glm::mat4 model_rr = glm::translate(model_r, glm::vec3(-1.0f*i*scale,scale*static_cast<float>(layer),0.0f));
                shader.setMat4("model", glm::scale(model_rr, glm::vec3(scale)));
                lastQuad = renderQuad();
            }
        } else if (drc == E){
            glm::mat4 model_l = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
            model_l = glm::translate(model_l, glm::vec3(-move_to.z,move_to.y, move_to.x));
        }
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fbTexture);
        //glActiveTexture(GL_TEXTURE1);
        shader.setVec2("texRes", tex_dims[fbTexture]);

        glm::mat4 model_f = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
        model_f = glm::rotate(model_f, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 model_b = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
        model_b = glm::rotate(model_b, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        model_b = glm::translate(model_b, glm::vec3(0.0f,0.0f,-scale*static_cast<float>(len_f_to_b)));
        model_f = glm::translate(model_f, glm::vec3(-scale*0.5f,scale*0.0f,-scale*0.5f));
        model_b = glm::translate(model_b, glm::vec3(scale*0.5f,scale*0.0f,scale*0.5f));

        model_f = glm::translate(model_f, glm::vec3(-move_to.x,move_to.y,-move_to.z));
        model_b = glm::translate(model_b, glm::vec3(move_to.x,move_to.y,move_to.z));

        for(int j = 0; j < len_l_to_r; ++j){

            glm::mat4 model_ff = glm::translate(model_f, glm::vec3(-1.0f*j*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_ff, glm::vec3(scale)));
            lastQuad = renderQuad();

            glm::mat4 model_bb = glm::translate(model_b, glm::vec3(1.0f*j*scale,scale*static_cast<float>(layer),0.0f));
            shader.setMat4("model", glm::scale(model_bb, glm::vec3(scale)));
            lastQuad = renderQuad();
        }
        //pass_door_walls = true;
    }

    return lastQuad;
}

unsigned int renderCube(Shader& shader, int height, int m, int n, unsigned int lr_dMap, unsigned int fb_dMap,
                        unsigned int f_dMap, unsigned int c_dMap,
                        TexDims& tex_dims, glm::vec3 move_to, glm::vec3 viewPos, glm::mat4 view,
                        float scale, bool& pass_door_walls)
{

    //for(int i = 0; i < 4; ++i){
    //    shader.setVec3("pointLights["+std::to_string(i)+"].position", 3*i,  0.6f, 0.3f);
    //}
    //float il = m / num_lights;

    unsigned int lastQuad = 0;
    float quad_sz = 1.5;
    glm::vec3 lightPos(10.5f, 1.0f, 0.3f);        /////////////////
    glm::mat4 model = glm::mat4(1.0f);
    //model = glm::rotate(model, glm::radians((float)glfwGetTime() * -10.0f), glm::normalize(glm::vec3(1.0, 0.0, 1.0))); // rotate the quad to show normal mapping from multiple directions
    shader.setMat4("model", model);
    shader.setMat4("view", view);
    shader.setVec3("viewPos", viewPos);
    shader.setVec3("lightPos", lightPos);


    //glStencilFunc(GL_NOTEQUAL, 1, 0xFF);
    //glStencilMask(0x00);

    for(int k = 0; k < height; ++k){
        //glStencilFunc(GL_NOTEQUAL, 1, 0xFF);
        //glStencilMask(0x00);
        renderCubeLayer(shader, k, lr_dMap, fb_dMap, m, n, tex_dims, move_to, viewPos, view, scale, pass_door_walls);
        //glBindVertexArray(0);
        //glStencilMask(0xFF);
        //glStencilFunc(GL_ALWAYS, 0, 0xFF);
        //glEnable(GL_DEPTH_TEST);
        //renderCubeLayer(shader, k, dMap, nMap, m, n, move_to, scale, pass_door_walls);
    }

    if(!pass_door_walls){

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, c_dMap);
        //glActiveTexture(GL_TEXTURE1);
        //glBindTexture(GL_TEXTURE_2D, c_nMap);
        shader.setVec2("texRes", tex_dims[c_dMap]);

        for(int i = 0; i < m; ++i){
            for(int j = 0; j < n; ++j){
                glm::mat4 model_cl = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
                model_cl = glm::rotate(model_cl, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                model_cl = glm::translate(model_cl, glm::vec3(0.0f,scale*0.5f,-scale*0.5));
                model_cl = glm::translate(model_cl, glm::vec3(scale*static_cast<float>(i),scale*static_cast<float>(j),-scale*static_cast<float>(height-1)));
                model_cl = glm::translate(model_cl, glm::vec3(-move_to.z, move_to.x, -move_to.y));
                model_cl = glm::scale(model_cl, glm::vec3(scale));
                shader.setMat4("model", model_cl);
                lastQuad = renderQuad();
            }
        }

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, f_dMap);
        //glActiveTexture(GL_TEXTURE1);
        //glBindTexture(GL_TEXTURE_2D, f_nMap);
        shader.setVec2("texRes", tex_dims[f_dMap]);

        for(int i = 0; i < m; ++i){
            for(int j = 0; j < n; ++j){
                glm::mat4 model_fl = glm::translate(model, glm::vec3(0.0f,0.0f,0.0f));
                model_fl = glm::rotate(model_fl, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                model_fl = glm::rotate(model_fl, glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                model_fl = glm::translate(model_fl, glm::vec3(0.0f,scale*0.5f,-scale*0.5));
                model_fl = glm::translate(model_fl, glm::vec3(-scale*static_cast<float>(i),scale*static_cast<float>(j),0.0f));
                model_fl = glm::translate(model_fl, glm::vec3(move_to.z, move_to.x, move_to.y));
                model_fl = glm::scale(model_fl, glm::vec3(scale));
                shader.setMat4("model", model_fl);
                lastQuad = renderQuad();
            }
        }
    }
    //glBindVertexArray(0);
    //glStencilMask(0xFF);
    //glStencilFunc(GL_ALWAYS, 0, 0xFF);
    //glEnable(GL_DEPTH_TEST);
    return lastQuad;
}


// renders a 1x1 quad in NDC with manually calculated tangent vectors
// ------------------------------------------------------------------
void _renderQuad(Shape shape)
{

    if(quadVAO == 0){
        glGenVertexArrays(1, &quadVAO);
        glGenBuffers(1, &quadVBO);


        // positions
        glm::vec3 pos1(-0.5f,  0.5f, 0.0f);
        glm::vec3 pos2(-0.5f, -0.5f, 0.0f);
        glm::vec3 pos3( 0.5f, -0.5f, 0.0f);
        glm::vec3 pos4( 0.5f,  0.5f, 0.0f);

        glm::vec3 pos5(-0.5f, 0.0f,  0.5f);
        glm::vec3 pos6(-0.5f, 0.0f, -0.5f);
        glm::vec3 pos7( 0.5f, 0.0f, -0.5f);
        glm::vec3 pos8( 0.5f, 0.0f,  0.5f);

        glm::vec3 pos9( -0.5f, 0.0f,  0.5f);
        glm::vec3 pos10(-0.5f, 0.0f, -0.5f);
        glm::vec3 pos11( 0.5f, 0.0f, -0.5f);
        glm::vec3 pos12( 0.5f, 0.0f,  0.5f);


        // texture coordinates
        glm::vec2 uv1(0.0f, 1.0f);
        glm::vec2 uv2(0.0f, 0.0f);
        glm::vec2 uv3(1.0f, 0.0f);
        glm::vec2 uv4(1.0f, 1.0f);

        // normal vector
        glm::vec3 nm(0.0f, 0.0f, 1.0f);

        glm::vec3 nm2(0.0f, 1.0f, 0.0f);

        glm::vec3 nm3(0.0f, 1.0f, 0.0f);

        // calculate tangent/bitangent vectors of both triangles
        glm::vec3 tangent1, bitangent1;
        glm::vec3 tangent2, bitangent2;

        glm::vec3 tangent3, bitangent3;
        glm::vec3 tangent4, bitangent4;

        glm::vec3 tangent5, bitangent5;
        glm::vec3 tangent6, bitangent6;

        // triangle 1
        // ----------
        glm::vec3 edge1 = (pos2 - pos1);
        glm::vec3 edge2 = (pos3 - pos1);
        glm::vec2 deltaUV1 = (uv2 - uv1);
        glm::vec2 deltaUV2 = (uv3 - uv1);

        float f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

        tangent1.x = f * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
        tangent1.y = f * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
        tangent1.z = f * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);

        bitangent1.x = f * (-deltaUV2.x * edge1.x + deltaUV1.x * edge2.x);
        bitangent1.y = f * (-deltaUV2.x * edge1.y + deltaUV1.x * edge2.y);
        bitangent1.z = f * (-deltaUV2.x * edge1.z + deltaUV1.x * edge2.z);

        // triangle 2
        // ----------
        edge1 = pos3 - pos1;
        edge2 = pos4 - pos1;
        deltaUV1 = uv3 - uv1;
        deltaUV2 = uv4 - uv1;

        f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

        tangent2.x = f * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
        tangent2.y = f * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
        tangent2.z = f * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);


        bitangent2.x = f * (-deltaUV2.x * edge1.x + deltaUV1.x * edge2.x);
        bitangent2.y = f * (-deltaUV2.x * edge1.y + deltaUV1.x * edge2.y);
        bitangent2.z = f * (-deltaUV2.x * edge1.z + deltaUV1.x * edge2.z);


        float quadVertices[] = {
                    // positions            // normal         // texcoords  // tangent                          // bitangent
                    pos1.x, pos1.y, pos1.z, nm.x, nm.y, nm.z, uv1.x, uv1.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,
                    pos2.x, pos2.y, pos2.z, nm.x, nm.y, nm.z, uv2.x, uv2.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,
                    pos3.x, pos3.y, pos3.z, nm.x, nm.y, nm.z, uv3.x, uv3.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,

                    pos1.x, pos1.y, pos1.z, nm.x, nm.y, nm.z, uv1.x, uv1.y, tangent2.x, tangent2.y, tangent2.z, bitangent2.x, bitangent2.y, bitangent2.z,
                    pos3.x, pos3.y, pos3.z, nm.x, nm.y, nm.z, uv3.x, uv3.y, tangent2.x, tangent2.y, tangent2.z, bitangent2.x, bitangent2.y, bitangent2.z,
                    pos4.x, pos4.y, pos4.z, nm.x, nm.y, nm.z, uv4.x, uv4.y, tangent2.x, tangent2.y, tangent2.z, bitangent2.x, bitangent2.y, bitangent2.z,
        };
        float triangleVertices[] = {
                    // positions            // normal         // texcoords  // tangent                          // bitangent
                    pos1.x, pos1.y, pos1.z, nm.x, nm.y, nm.z, uv1.x, uv1.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,
                    pos2.x, pos2.y, pos2.z, nm.x, nm.y, nm.z, uv2.x, uv2.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,
                    pos3.x, pos3.y, pos3.z, nm.x, nm.y, nm.z, uv3.x, uv3.y, tangent1.x, tangent1.y, tangent1.z, bitangent1.x, bitangent1.y, bitangent1.z,
        };

        // configure plane VAO
        if(shape == Quad){
            glBindVertexArray(quadVAO);
            glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
            glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(6 * sizeof(float)));
            glEnableVertexAttribArray(3);
            glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(8 * sizeof(float)));
            glEnableVertexAttribArray(4);
            glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(11 * sizeof(float)));
        } else {
            glBindVertexArray(quadVAO);
            glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
            glBufferData(GL_ARRAY_BUFFER, sizeof(triangleVertices), &triangleVertices, GL_STATIC_DRAW);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(6 * sizeof(float)));
            glEnableVertexAttribArray(3);
            glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(8 * sizeof(float)));
            glEnableVertexAttribArray(4);
            glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(11 * sizeof(float)));
        }
    }
}

unsigned int renderQuad()
{
    //  unsigned int quadVAO = 0;
    // unsigned int quadVBO;
    _renderQuad(Quad);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
    return quadVAO;
}

void renderWdws(float sz, Shader shader, int num_doors, int n, int m, int h,
                unsigned int mirrorTex, unsigned int mirror_normalMap, glm::vec3 move_to,
                glm::mat4 view, glm::mat4 projection){

}

bool flag = false;
void renderDoors(float width, Shader& shader, Shader& normal_shader, int num_doors, int n, int m, int h,
                 unsigned int outer_brickTex,
                 unsigned int inner_brickTex,
                 TexDims& tex_dims, bool grid_doors,
                 glm::vec3 move_to, glm::mat4 view, glm::mat4 projection, glm::vec3 viewPos)
{

    int w = 0.3*(m/num_doors);
    glm::vec3 move_t = move_to;
    //renderCube(shader, 0.6, w, 0.3, mirrorTex, mirror_normalMap,
    //          mirror_normalMap, mirrorTex, move_to);
    float scale = 0.1f;
    float border_scale = 0.02;
    if(!grid_doors){
        move_to.y -= 0.5*(1-scale);
    }
    int m_w = 15; int n_w = 1;
    int wall_height = 14;
    if(grid_doors){
        scale = 0.155;
    }

    glStencilMask(0xFF);
    glStencilFunc(GL_ALWAYS, 1, 0xFF);
    //glDisable(GL_DEPTH_TEST);
    //glDepthFunc(GL_LESS);
    shader.use();
    shader.setVec3("viewPos", camera.Position);
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);

    shader.setVec3("dirLight.direction", -6.2f, -5.0f, -6.3f);
    shader.setVec3("dirLight.ambient", 0.05f, 0.05f, 0.05f);
    shader.setVec3("dirLight.diffuse", 5.4f, 5.4f, 5.4f);
    shader.setVec3("dirLight.specular", 5.5f, 5.5f, 5.5f);


    bool pass_door_walls = false;
    renderCube(shader, 15, width, 1,
               outer_brickTex, inner_brickTex,
               outer_brickTex, outer_brickTex, tex_dims,
               move_to, viewPos, view,scale-border_scale, pass_door_walls);
    pass_door_walls = true;
    renderCube(shader, 15, width, 1,
               outer_brickTex, inner_brickTex,
               outer_brickTex, outer_brickTex, tex_dims,
               move_to, viewPos, view, scale-border_scale, pass_door_walls);
    //glDepthFunc(GL_ALWAYS);
    glStencilFunc(GL_NOTEQUAL, 1, 0xFF);
    glStencilMask(0x00);
    //lglDisable(GL_DEPTH_TEST);
    //glDepthFunc(GL_LESS);

    normal_shader.use();
    normal_shader.setVec3("viewPos", camera.Position);
    normal_shader.setMat4("view", view);
    unsigned int lastQuad = 0;

    normal_shader.setVec3("dirLight.direction", -0.2f, -1.0f, -0.3f);
    normal_shader.setVec3("dirLight.ambient", 0.05f, 0.05f, 0.05f);
    normal_shader.setVec3("dirLight.diffuse", 0.4f, 0.4f, 0.4f);
    normal_shader.setVec3("dirLight.specular", 0.5f, 0.5f, 0.5f);

    normal_shader.setInt("nr_point_lights", 8);
    for(int i = 0; i < 8; ++i){
        shader.setVec3("pointLights["+std::to_string(i)+"].ambient", 0.05f, 0.05f, 0.05f);
        shader.setVec3("pointLights["+std::to_string(i)+"].diffuse", 0.8f, 0.8f, 0.8f);
        shader.setVec3("pointLights["+std::to_string(i)+"].specular", 1.0f, 1.0f, 1.0f);
        shader.setFloat("pointLights["+std::to_string(i)+"].constant", 1.0f);
        shader.setFloat("pointLights["+std::to_string(i)+"].linear", 0.09f);
        shader.setFloat("pointLights["+std::to_string(i)+"].quadratic", 0.032f);
    }
    // spotLight
    normal_shader.setVec3("spotLight.position", camera.Position);
    normal_shader.setVec3("spotLight.direction", camera.Front);
    normal_shader.setVec3("spotLight.ambient", 0.3f, 0.3f, 0.3f);
    normal_shader.setVec3("spotLight.diffuse", 1.0f, 1.0f, 1.0f);
    normal_shader.setVec3("spotLight.specular", 1.0f, 1.0f, 1.0f);
    normal_shader.setFloat("spotLight.constant", 1.0f);
    normal_shader.setFloat("spotLight.linear", 0.09f);
    normal_shader.setFloat("spotLight.quadratic", 0.032f);
    normal_shader.setFloat("spotLight.cutOff", glm::cos(glm::radians(12.5f)));
    normal_shader.setFloat("spotLight.outerCutOff", glm::cos(glm::radians(15.0f)));

    // material properties
    normal_shader.setFloat("shininess", 52.0f);


    //if(!flag){
    //    std::cout << "moveto z " << move_to.z << '\n';
    //   // std::cout << "w*scale " << scale*(move_to.z + width) << '\n';
    //    //std::cout << "w" << (scale-0.02)*(move_to.z + width) << '\n';
    //    std::cout << (move_to.z + width)*(border_scale/2) << '\n';
    //    flag = true;
    //}

    // move_to.z *= (1+border_scale);
    move_to.z += width*(border_scale/2);
    //move_to.z += (width*border_scale)/2;
    pass_door_walls = false;
    renderCube(normal_shader, 13, width, 1, outer_brickTex, outer_brickTex,
               outer_brickTex, outer_brickTex, tex_dims, move_to, viewPos, view, scale, pass_door_walls);
    pass_door_walls = true;
    renderCube(normal_shader, 13, width, 1, outer_brickTex, outer_brickTex,
               outer_brickTex, outer_brickTex, tex_dims, move_to, viewPos, view, scale, pass_door_walls);

    move_to.x -= 0.1;
    move_to.z -= width*(border_scale/4);
    pass_door_walls = false;
    renderCube(normal_shader, 13, width, 1, outer_brickTex, outer_brickTex,
               outer_brickTex, outer_brickTex, tex_dims, move_to, viewPos, view, scale, pass_door_walls);
    pass_door_walls = true;
    renderCube(normal_shader, 13, width, 1, outer_brickTex, outer_brickTex,
               outer_brickTex, outer_brickTex, tex_dims, move_to, viewPos, view, scale, pass_door_walls);
    // glBindVertexArray(0);
    // glStencilMask(0x00);
    //glStencilFunc(GL_ALWAYS, 0, 0x00);

    glBindVertexArray(0);
    glStencilMask(0xFF);
    glStencilFunc(GL_ALWAYS, 0, 0xFF);
    glEnable(GL_DEPTH_TEST);

    //move_t.z += w;
    //for(int i = 0; i < num_doors; ++i){
    //    move_t.z += (m/num_doors);
    //    renderCube(shader, 0.6, w, 0.3, mirrorTex, mirror_normalMap,
    //               mirror_normalMap, mirrorTex, move_t);
    //}
}

void renderMirrors(float sz, Shader& shader, int num_mirrors, int n, int m, unsigned int textureColorbuffer,
                   glm::vec3 move_to, glm::mat4 view, glm::vec3 viewPos, glm::mat4 projection,
                   bool is_grid, Direction drc, int level)
{

    int lvl = is_grid ? level : 0;
    initShader(shader,2,view, is_grid, m, all_grids[lvl]->grid);
    //shader.use();
    //shader.setVec3("dirLight.direction", -6.2f, -5.0f, -6.3f);
    //shader.setVec3("dirLight.ambient", 0.05f, 0.05f, 0.05f);
    //shader.setVec3("dirLight.diffuse", 5.4f, 5.4f, 5.4f);
    //shader.setVec3("dirLight.specular", 5.5f, 5.5f, 5.5f);
    //shader.setFloat("shininess", 100.0f);
    //shader.setMat4("view", view);
    //shader.setMat4("projection", projection);


    float scale = sz;
    float il = n/num_mirrors;
    float il2 = il/2;
    glm::mat4 model = glm::mat4(1.0f);
    if(!is_grid){ // n=length,m=width
        model = glm::translate(model, glm::vec3(-move_to.z,move_to.y,-move_to.x));
        model = glm::translate(model, glm::vec3(0.0f, 0.05f,static_cast<float>(m)+0.1));
        for(int k = 0; k < num_mirrors; ++k){
            float il_pos = (il*k)+il2;
            float c;
            _renderQuad(Quad);
            glBindVertexArray(quadVAO);
            //  glm::mat4 model_b = glm::translate(model, glm::vec3(-il_pos,  0.0f,0.0f));
            glm::mat4 model_b = glm::translate(model, glm::vec3(-il_pos,  0.0f,0.0f));
            model_b = glm::scale(model_b, glm::vec3(scale));
            //static_cast<float>(n-0.01)
            shader.setMat4("model", model_b);
            //glActiveTexture(GL_TEXTURE0);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
            glDrawArrays(GL_TRIANGLES, 0, 6);
            glBindVertexArray(0);
        }
    } else { //n=cell size, m= grid_size
        float wall_scale = 0.155;
        //for(int r = 0; r < m; ++r) {
         //   for(int c = 0; c < m; ++c){
                switch(drc){
                    case N : {
                     //   move_to.z -= (m_crd.first+0.25)*n*wall_scale;
                        move_to.z += (n/4)*wall_scale;
                        move_to.x += 3*wall_scale;
                        move_to.x += 0.1f*wall_scale;
                        _renderQuad(Quad);
                        glBindVertexArray(quadVAO);
                        glm::mat4 model_rc = glm::translate(model, glm::vec3(-move_to.z,move_to.y,-move_to.x));
                        shader.setMat4("model", model_rc);
                        glActiveTexture(GL_TEXTURE0);
                        glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
                        // glActiveTexture(GL_TEXTURE1);
                        glDrawArrays(GL_TRIANGLES, 0, 6);
                    }
                    case S : {
                     //   move_to.z -= (c+0.25)*n*wall_scale;
                        //move_to.x += n*wall_scale;
                        move_to.z += (n/4)*wall_scale;
                        move_to.x -= 3*wall_scale;
                        move_to.x -= 0.1f*wall_scale;
                        _renderQuad(Quad);
                        glBindVertexArray(quadVAO);
                        glm::mat4 model_rc = glm::translate(model, glm::vec3(-move_to.z,move_to.y,-move_to.x));
                        shader.setMat4("model", model_rc);
                        glActiveTexture(GL_TEXTURE0);
                        glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
                        // glActiveTexture(GL_TEXTURE1);
                        glDrawArrays(GL_TRIANGLES, 0, 6);
                    }
                }
         //   }
        //}
    }
}

void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD,
                                   deltaTime );
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD,
                                   deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT,
                                   deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT,
                                   deltaTime);

    auto Position = camera.Position;
    // std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
}


//void processInput(GLFWwindow *window, PathTo& path, Grid& grid)
//{
//    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
//        glfwSetWindowShouldClose(window, true);
//
//    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
//        camera.ProcessKeyboardRoom(FORWARD,
//                                   deltaTime, 15, 0.155f,
//                                   path.first, grid);
//    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
//        camera.ProcessKeyboardRoom(BACKWARD,
//                                   deltaTime, 15, 0.155f,
//                                   path.first, grid);
//    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
//        camera.ProcessKeyboardRoom(LEFT,
//                                   deltaTime, 15, 0.155,
//                                   path.first, grid);
//    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
//        camera.ProcessKeyboardRoom(RIGHT,
//                                   deltaTime, 15, 0.155,
//                                   path.first, grid);
//    auto Position = camera.Position;
//   // std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
//}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    bool l_down,r_down;
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if(GLFW_PRESS == action)
            l_down= true;
        else if(GLFW_RELEASE == action)
            l_down= false;
    }
    if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if(GLFW_PRESS == action)
            r_down= true;
        else if(GLFW_RELEASE == action)
            r_down= false;
    }

    if(r_down) {
        if(!mouse_disabled){
            glfwSetCursorPosCallback(window, NULL);
        } else {
            glfwSetCursorPosCallback(window, mouse_callback);
        }
    }
    if(l_down) {
        if(!mouse_disabled){
         /////////
        }
    }
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{

    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}

// utility function for loading a 2D texture from file
// ---------------------------------------------------
unsigned int loadTexture(char const * path)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    unsigned char *data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, format == GL_RGBA ? GL_CLAMP_TO_EDGE : GL_LINEAR); // for this tutorial: use GL_CLAMP_TO_EDGE to prevent semi-transparent borders. Due to interpolation it takes texels from next repeat
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, format == GL_RGBA ? GL_CLAMP_TO_EDGE : GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}

// loads a cubemap texture from 6 individual texture faces
// order:
// +X (right)
// -X (left)
// +Y (top)
// -Y (bottom)
// +Z (front)
// -Z (back)
// -------------------------------------------------------
unsigned int loadCubemap(vector<std::string> faces)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

    int width, height, nrChannels;
    for (unsigned int i = 0; i < faces.size(); i++)
    {
        unsigned char *data = stbi_load(faces[i].c_str(), &width, &height, &nrChannels, 0);
        if (data)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);
        }
        else
        {
            std::cout << "Cubemap texture failed to load at path: " << faces[i] << std::endl;
            stbi_image_free(data);
        }
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    return textureID;
}