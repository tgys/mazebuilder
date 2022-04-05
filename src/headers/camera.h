#ifndef CAMERA_H
#define CAMERA_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>
#include <set>

// Defines several possible options for camera movement. Used as abstraction to stay away from window-system specific input methods
enum Camera_Movement {
    FORWARD,
    BACKWARD,
    LEFT,
    RIGHT,
    UP,
    DOWN
};

enum Direction{
    N = 0,
    S = 1,
    E = 2,
    W = 3
};

// Default camera values
const float YAW         = -90.0f;
const float PITCH       =  0.0f;
const float SPEED       =  9.5f;
const float SENSITIVITY =  0.1f;
const float ZOOM        =  45.0f;

class cubeDims{
    int h; int n; int m;
public:
    glm::vec3 move_to;
    cubeDims(int h,int n,int m,glm::vec3 mv_to) :
    h(h), n(n), m(m), move_to(mv_to) {}
    int height(){ return h; }
    int N(){ return n; }
    int M(){ return m; }
};

// An abstract camera class that processes input and calculates the corresponding Euler Angles, Vectors and Matrices for use in OpenGL
class Camera
{
public:
    // camera Attributes
    glm::vec3 Position;
    glm::vec3 Front;
    glm::vec3 Up;
    glm::vec3 Right;
    glm::vec3 WorldUp;
    // euler Angles
    float Yaw;
    float Pitch;
    // camera options
    float MovementSpeed;
    float MouseSensitivity;
    float Zoom;

    std::vector<std::unique_ptr<cubeDims>> room_map;
    std::vector<Camera*> cam_mp;

    void AddRoom(int h, int n, int m, glm::vec3 move_to){
        room_map.push_back(std::make_unique<cubeDims>(h,n,m,move_to));
    }

    // constructor with vectors
    Camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
    {
        Position = position;
        WorldUp = up;
        Yaw = yaw;
        Pitch = pitch;
        updateCameraVectors();
    }
    // constructor with scalar values
    Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
    {
        Position = glm::vec3(posX, posY, posZ);
        WorldUp = glm::vec3(upX, upY, upZ);
        Yaw = yaw;
        Pitch = pitch;
        updateCameraVectors();
    }

    // returns the view matrix calculated using Euler Angles and the LookAt Matrix
    glm::mat4 GetViewMatrix()
    {
        return glm::lookAt(Position, Position + Front, Up);
    }

    glm::mat4 GetViewMatrixPos(glm::vec3 pos)
    {
        return glm::lookAt(pos, pos + Front, Up);
    }

    // processes input received from any keyboard-like input system. Accepts input parameter in the form of camera defined ENUM (to abstract it from windowing systems)
    void ProcessKeyboard(Camera_Movement direction, float deltaTime)
    {
        float velocity = MovementSpeed * deltaTime;
        if (direction == FORWARD)
            Position += Front * velocity;
        if (direction == BACKWARD)
            Position -= Front * velocity;
        if (direction == LEFT)
            Position -= Right * velocity;
        if (direction == RIGHT)
            Position += Right * velocity;
        if (direction == UP)
            Position += Up * velocity;
        if (direction == DOWN)
            Position -= Up * velocity;
    }

    bool check_in_box(glm::vec3& pos, std::unique_ptr<cubeDims>& cb){
        std::cout << "pos " << '\n';
        //std::cout << pos.x << ' ' << pos.y << ' ' << pos.z << '\n';

        glm::vec3 tpos = pos + cb->move_to;
        std::cout << "tpos " << '\n';
        std::cout << tpos.x << ' ' << tpos.y << ' ' << tpos.z << '\n';

        float c_hm = cb->move_to.y + cb->height();
        //float c_nm = cb->move_to.x + cb->N();
        float c_nm = (-cb->move_to.z) + cb->N();
        //float c_mm = cb->move_to.z + cb->M();
        float c_mm = cb->move_to.x + cb->M();
        std::cout << "c_hm " << cb->move_to.y + cb->height() << '\n';
        std::cout << "c_nm " << cb->move_to.x + cb->N() << '\n';
        std::cout << "c_mm " << cb->move_to.z + cb->M() << '\n';
        std::cout << "c_h " << cb->move_to.y << '\n';
        std::cout << "c_n " << cb->move_to.x << '\n';
        std::cout << "c_m " << cb->move_to.z << '\n';

        if(tpos.z < (cb->move_to.x + 0.01) && tpos.x > ((-cb->move_to.z) + 0.01) && tpos.x < ((-cb->move_to.z) + 1.5)){
            std::cout << "checkinbox true "  << '\n';
            return true;
        } else {
            if(//tpos.y > (cb->move_to.y + 0.01) && tpos.y < (c_hm - 0.01) &&
            (tpos.x > ((-cb->move_to.z) + 0.01) && tpos.x < (c_nm - 0.01) &&
               tpos.z > (cb->move_to.z + 0.01) && tpos.z < (c_mm - 0.01))
               || tpos.y > 2){
                std::cout << "checkinbox true notdoor "  << '\n';
                  return true;
            }
        }
        std::cout << "checkinbox false "  << '\n';
        return false;
    }

    typedef std::pair<int ,int> Coord;
    typedef std::set<Direction> Cell;
    typedef std::vector<std::vector<Cell>> Grid;

    bool check_in_room(glm::vec3& npos){
        int idx = 0;
        for(auto& cptr : room_map){

            if(check_in_box(npos, cptr)){
                std::cout << "checkinroom true "  << '\n';
                std::cout << "idx " << idx << '\n';
                return true;
            }
            idx++;
        }
        std::cout << "checkinroom false "  << '\n';
        //return false;
        return true;
    }

    bool check_in_grid(glm::vec3& npos, Grid& grid, int m, float wall_scale){
        int c_start = std::floor((npos.x/m)/wall_scale);
        int r_start = std::floor((npos.z/m)/wall_scale);

        int n_bound = r_start*m*wall_scale;
        int w_bound = c_start*m*wall_scale;
        int e_bound = (c_start+1)*m*wall_scale;
        int s_bound = (r_start+1)*m*wall_scale;

        bool in_grid = true;
        if(grid[c_start][r_start].contains(N)){
            if(npos.x > (n_bound-wall_scale)){
                in_grid = false;
            }
        }
        if(grid[c_start][r_start].contains(S)){
            if(npos.x < (s_bound+wall_scale)){
                in_grid = false;
            }
        }
        if(grid[c_start][r_start].contains(W)){
            if(npos.z < (w_bound+wall_scale)){
                in_grid = false;
            }
        }
        if(grid[c_start][r_start].contains(E)){
            if(npos.z > (e_bound-wall_scale)){
                in_grid = false;
            }
        }

        return in_grid;

        std::cout << "CSTART: " << c_start << ' ' << "RSTART: " <<  r_start << '\n';
    }


    void ProcessKeyboardRoom(Camera_Movement direction,
                             float deltaTime, int m, float wall_scale,
                             std::vector<std::pair<Coord,Direction>>& path,
                             Grid& grid)
    {
        float velocity = MovementSpeed * deltaTime;
        glm::vec3 new_position = Position;
        if (direction == FORWARD){
            new_position = Position + Front * velocity;
            //std::cout << "front "  << Front.x << ' ' << Front.y << ' ' << Front.z << '\n';
            if(check_in_grid(new_position, grid, m, wall_scale)){
                Position = new_position;
            }
            //if(check_in_room(new_position)){
            //    std::cout << "position changed"  << '\n';
            //    std::cout << "old pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //    Position = new_position;
            //    std::cout << "new pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //}
        }
        if (direction == BACKWARD){
            new_position = Position - Front * velocity;
            //std::cout << "front "  << Front.x << ' ' << Front.y << ' ' << Front.z << '\n';
            if(check_in_grid(new_position, grid, m, wall_scale)){
                Position = new_position;
            }
            //if(check_in_room(new_position)){
            //    std::cout << "position changed"  << '\n';
            //    std::cout << "old pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //    Position = new_position;
            //    std::cout << "new pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //}
        }
        if (direction == LEFT){
            new_position = Position - Right * velocity;
            //if(check_in_room(new_position)){
            if(check_in_grid(new_position, grid, m, wall_scale)){
                Position = new_position;
            }
            //    std::cout << "position changed"  << '\n';
            //    std::cout << "old pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //    Position = new_position;
            //    std::cout << "new pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //}
        }
        if (direction == RIGHT){
            new_position = Position + Right * velocity;
            //if(check_in_room(new_position)){
            if(check_in_grid(new_position, grid, m, wall_scale)){
                Position = new_position;
            }
            //    std::cout << "position changed"  << '\n';
            //    std::cout << "old pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //    Position = new_position;
            //    std::cout << "new pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //}
        }
        if (direction == UP){
            new_position = Position + Up * velocity;
            //if(check_in_room(new_position)){
            if(check_in_grid(new_position, grid, m, wall_scale)){
                Position = new_position;
            }
            //    std::cout << "position changed"  << '\n';
            //    std::cout << "old pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //    Position = new_position;
            //    std::cout << "new pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //}
        }
        if (direction == DOWN){
            new_position = Position - Up * velocity;
            //if(check_in_room(new_position)){
            if(check_in_grid(new_position, grid, m, wall_scale)){
                Position = new_position;
            }
            //    std::cout << "position changed"  << '\n';
            //    std::cout << "old pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //    Position = new_position;
            //    std::cout << "new pos "  << '\n';
            //    std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
            //}
        }
        std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << '\n';
    }

    // processes input received from a mouse input system. Expects the offset value in both the x and y direction.
    void ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch = true)
    {
        xoffset *= MouseSensitivity;
        yoffset *= MouseSensitivity;

        Yaw   += xoffset;
        //std::cout << Yaw << '\n';
        Pitch += yoffset;

        // make sure that when pitch is out of bounds, screen doesn't get flipped
        if (constrainPitch)
        {
            if (Pitch > 89.0f)
                Pitch = 89.0f;
            if (Pitch < -89.0f)
                Pitch = -89.0f;
        }

        // update Front, Right and Up Vectors using the updated Euler angles
        updateCameraVectors();
    }

    // processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
    void ProcessMouseScroll(float yoffset)
    {
        Zoom -= (float)yoffset;
        if (Zoom < 1.0f)
            Zoom = 1.0f;
        if (Zoom > 45.0f)
            Zoom = 45.0f; 
    }

private:
    // calculates the front vector from the Camera's (updated) Euler Angles
    void updateCameraVectors()
    {
        // calculate the new Front vector
        glm::vec3 front;
        front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
        front.y = sin(glm::radians(Pitch));
       // front.y = 0;
        front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
        Front = glm::normalize(front);
        // also re-calculate the Right and Up vector
        Right = glm::normalize(glm::cross(Front, WorldUp));  // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
        Up    = glm::normalize(glm::cross(Right, Front));
    }
};
#endif