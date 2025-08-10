#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "container.h"
#include "particle.h"

#include <array>
#include <iostream>
using namespace std;

int num_particles = 125; //use perfect cubes to make life easier
float particle_radius = 0.06;
float drawing_radius = 0.02;
glm::vec3 gravity = glm::vec3(0.0, -.000005, 0.0);
//Vector3D gravity = Vector3D(0.0, 0.0, 0.0);


// for mouse dragging
bool gDragging = false;
double gLastX = 0.0, gLastY = 0.0;
float gSensitivity = 0.25f;


unsigned int shaderProgram{};
glm::mat4 rotationMatrix = glm::mat4(1.0f);
Container container;

GLuint VBO, VAO;
std::vector<Particle*> particles;

float lastFrameStartTime = 0.0f;

constexpr auto vertexShaderSource = R"(
    #version 330 core
    
    layout (location = 0) in vec3 aPos;
    uniform mat4 rotationMatrix;

    void main()
    {
        gl_Position =  rotationMatrix * vec4(aPos.x, aPos.y, aPos.z, 1.0);
    }
)";

constexpr auto fragmentShaderSource = R"(
    #version 330 core

    out vec4 FragColor;

    void main()
    {
        FragColor = vec4(1.0f, 1.0f, 1.0f, 1.0f);
    }
)";

void framebufferSizeChanged(GLFWwindow *window, int width, int height);
void processInput(GLFWwindow *window, float deltaTime);
void render(GLFWwindow *window, float deltaTime);
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
void addOneParticle();
void removeLastParticle();


void initializeParticles() {
    //free memory associated with current particles if any 
    for (int i = 0; i < particles.size(); i++) {
        delete particles[i];
    }

    if (num_particles == 1) {
        Particle* p = new Particle(0.0, 0.0, 0.0, particle_radius, drawing_radius);
        particles.push_back(p);
        return;
    }

    int dim = cbrt(num_particles);

    if (dim % 2 == 1) {
        float min_x = -0 - ((dim - 1) * particle_radius);
        float max_x = 0 + (((dim - 1)) * particle_radius);
        float min_y = -0 - (((dim - 1)) * particle_radius);
        float max_y = 0 + (((dim - 1)) * particle_radius);
        float min_z = -0 - (((dim - 1)) * particle_radius);
        float max_z = 0 + (((dim - 1)) * particle_radius);
        for (float i = min_x; i <= max_x; i += 2 * particle_radius) {
            for (float j = min_y; j <= max_y; j += 2 * particle_radius) {
                for (float k = min_z; k <= max_z; k += 2 * particle_radius) {
                    Particle* p = new Particle(i, j, k, particle_radius, drawing_radius);
                    particles.push_back(p);
                }
            }
        }
    }
    else {
        float min_x = -particle_radius - (((dim / 2) - 1) * 2 * particle_radius);
        float max_x = particle_radius + (((dim / 2) - 1) * 2 * particle_radius);
        float min_y = -particle_radius - (((dim / 2) - 1) * 2 * particle_radius);
        float max_y = particle_radius + (((dim / 2) - 1) * 2 * particle_radius);
        float min_z = -particle_radius - (((dim / 2) - 1) * 2 * particle_radius);
        float max_z = particle_radius + (((dim / 2) - 1) * 2 * particle_radius);
        for (float i = min_x; i <= max_x; i += 2 * particle_radius) {
            for (float j = min_y; j <= max_y; j += 2 * particle_radius) {
                for (float k = min_z; k <= max_z; k += 2 * particle_radius) {
                    Particle* p = new Particle(i, j, k, particle_radius, drawing_radius);
                    particles.push_back(p);
                }
            }
        }
    }

}

int main()
{
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    const GLFWvidmode *VideoMode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    const auto ScreenHeight = VideoMode->height;

    GLFWwindow *window = glfwCreateWindow(ScreenHeight / 2, ScreenHeight / 2,
                                          "Waterbender", NULL, NULL);
    if (window == NULL)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);

    if (glewInit() != GLEW_OK)
    {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }

    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);

    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, nullptr, infoLog);
        std::cout << "Vertex Shader Compilation Failed:" << infoLog << std::endl;
    }

    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, nullptr, infoLog);
        std::cout << "Fragment Shader Compilation Failed:" << infoLog << std::endl;
    }

    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);

    if (!success)
    {
        glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
        std::cout << "Shader Program Linking Failed: %s" << infoLog << std::endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    initializeParticles();
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);


    glfwSetFramebufferSizeCallback(window, framebufferSizeChanged);


    while (!glfwWindowShouldClose(window))
    {
        float currentFrameStartTime = static_cast<float>(glfwGetTime());
        float deltaTime = currentFrameStartTime - lastFrameStartTime;
        lastFrameStartTime = currentFrameStartTime;

        processInput(window, deltaTime);
        render(window, deltaTime);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);

    glfwTerminate();
    return 0;
}
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            gDragging = true;
            glfwGetCursorPos(window, &gLastX, &gLastY);
            
        }
        else if (action == GLFW_RELEASE) {
            gDragging = false;
            
        }
    }
}

void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    if (!gDragging) return;
    double dx = xpos - gLastX;
    double dy = ypos - gLastY;
    gLastX = xpos; gLastY = ypos;
    float yawDeg = gSensitivity * static_cast<float>(dx);
    rotationMatrix = glm::rotate(rotationMatrix, glm::radians(yawDeg), glm::vec3(0.f, 1.f, 0.f));
    float pitchDeg = gSensitivity * static_cast<float>(dy);
    rotationMatrix = glm::rotate(rotationMatrix, glm::radians(pitchDeg), glm::vec3(1.f, 0.f, 0.f));
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    if (yoffset > 0.0) {
        addOneParticle();
    }
    else if (yoffset < 0.0) {
        removeLastParticle();
    }
}
void addOneParticle() {
    //generate particle in a random location where x, y, and z are in the range [-0.1, 0.1]
    Particle* p = new Particle((((float)rand() / (RAND_MAX)) / 5) - 0.1f, (((float)rand() / (RAND_MAX)) / 5) - 0.1f, (((float)rand() / (RAND_MAX)) / 5) - 0.1f, particle_radius, drawing_radius);
    particles.push_back(p);
    num_particles = static_cast<int>(particles.size());
}

void removeLastParticle() {
    if (particles.empty()) return;
    delete particles.back();
    particles.pop_back();
    num_particles = static_cast<int>(particles.size());
}
void framebufferSizeChanged(GLFWwindow *window, int width, int height)
{
    float currentFrameStartTime = static_cast<float>(glfwGetTime());
    float deltaTime = currentFrameStartTime - lastFrameStartTime;
    lastFrameStartTime = currentFrameStartTime;
    glViewport(0, 0, width, height);
    render(window, deltaTime);
}

void processInput(GLFWwindow *window, float deltaTime)
{
    const float rotationSpeedDegreesPerSecond = 90.0f;
    const float rotationAnglePerFrame = rotationSpeedDegreesPerSecond * deltaTime;

    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
    {
        rotationMatrix = glm::rotate(rotationMatrix,
                                     glm::radians(-rotationAnglePerFrame),
                                     glm::vec3(0.0f, 1.0f, 0.0f));
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
    {
        rotationMatrix = glm::rotate(rotationMatrix,
                                     glm::radians(rotationAnglePerFrame),
                                     glm::vec3(0.0f, 1.0f, 0.0f));
    }
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
    {
        rotationMatrix = glm::rotate(rotationMatrix,
                                    glm::radians(rotationAnglePerFrame),
                                    glm::vec3(1.0f, 0.0f, 0.0f));
    }
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
    {
        rotationMatrix = glm::rotate(rotationMatrix,
                                    glm::radians(-rotationAnglePerFrame),
                                    glm::vec3(1.0f, 0.0f, 0.0f));
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) 
    {
        rotationMatrix = glm::mat4(1.0f);
    }
}

void render(GLFWwindow *window, float deltaTime)
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glUseProgram(shaderProgram);

    auto rotation_location = glGetUniformLocation(shaderProgram, "rotationMatrix");
    glUniformMatrix4fv(rotation_location, 1, GL_FALSE, glm::value_ptr(rotationMatrix));

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER,
        container.vertices.size() * sizeof(float),
        container.vertices.data(),
        GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(container.vertices.size() / 3));

    glm::mat4 invRotation = glm::transpose(rotationMatrix);
    glm::vec4 grav_vector = invRotation * glm::vec4(gravity.x, gravity.y, gravity.z, 1.0);
    for (int i = 0; i < particles.size(); i++) {
        std::vector<float> curr;
        particles[i]->updatePosition(&grav_vector, &container, &invRotation, &particles, deltaTime);
        particles[i]->draw(10, &curr, &invRotation);
        glBufferData(GL_ARRAY_BUFFER,
            curr.size() * sizeof(float),
            curr.data(),
            GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(curr.size() / 3));
    }

    glfwSwapBuffers(window);
}

