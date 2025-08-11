#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <unordered_set>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>

#include "container.h"
#include "particle.h"

using std::vector;

<<<<<<< HEAD
int   num_particles = 125;
float particle_radius = 0.06f;
float drawing_radius = 0.02f;
glm::vec3 gravity = glm::vec3(0.0f, -0.001f, 0.0f);

int   k_neighbors = 4;
float max_edge_dist = 0.35f;
=======
int num_particles = 125; //use perfect cubes to make life easier
float particle_radius = 0.06;
float drawing_radius = 0.02;
glm::vec3 gravity = glm::vec3(0.0, -.000005, 0.0);
>>>>>>> 2579bf8221f6876b715f34a7a0db0678eec93bbe

const float POS_ALPHA = 0.12f;   // 0..1, smaller -> more lag
const float NORM_ALPHA = 0.20f;   // EMA for normals
const float TOPOLOGY_DT = 0.08f;   // seconds between topology rebuilds

bool   gDragging = false;
double gLastX = 0.0, gLastY = 0.0;
float  gSensitivity = 0.25f;

glm::mat4 rotationMatrix = glm::mat4(1.0f);
Container container;
vector<Particle*> particles;


vector<glm::vec3> fluidPos;   // lagged particle positions
vector<glm::vec3> prevNorm;   // EMA per-vertex normals (indexed by particle index)
vector<glm::ivec3> fluidTris; // triangle index list
float topoTimer = 0.0f;

GLuint progTris = 0;            // filled mesh triangles
GLint  uRot_tris = -1;

GLuint progLines = 0;            // box wireframe
GLint  uRot_lines = -1;

GLuint trisVAO = 0, trisVBO = 0; // interleaved pos(3)+normal(3)
size_t trisBytes = 0;
GLsizei trisVertexCount = 0;

GLuint containerVAO = 0, containerVBO = 0;

float lastFrameStartTime = 0.0f;

static const char* VS_tris = R"(
#version 330 core
layout (location=0) in vec3 aPos;
layout (location=1) in vec3 aNor;
uniform mat4 uRot;
out vec3 vN;
void main() {
    gl_Position = uRot * vec4(aPos, 1.0);
    vN = mat3(uRot) * aNor; // rotation only
}
)";

static const char* FS_tris = R"(
#version 330 core
in vec3 vN;
out vec4 FragColor;
void main() {
    vec3 N = normalize(vN);
    vec3 L = normalize(vec3(0.45, 0.75, 0.5));
    float lambert = max(dot(N, L), 0.0);
    vec3 base = vec3(0.2, 0.6, 1.0);
    vec3 col = base * (0.2 + 0.8 * lambert);
    FragColor = vec4(col, 0.50); // semi-transparent
}
)";

static const char* VS_lines = R"(
#version 330 core
layout (location=0) in vec3 aPos;
uniform mat4 uRot;
void main(){ gl_Position = uRot * vec4(aPos,1.0); }
)";

static const char* FS_lines = R"(
#version 330 core
out vec4 FragColor;
void main(){ FragColor = vec4(1.0); }
)";

static GLuint compile(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) { char log[1024]; glGetShaderInfoLog(s, 1024, nullptr, log); std::cerr << "Shader compile error:\n" << log << "\n"; }
    return s;
}
static GLuint link(GLuint vs, GLuint fs) {
    GLuint p = glCreateProgram();
    glAttachShader(p, vs);
    glAttachShader(p, fs);
    glLinkProgram(p);
    GLint ok = 0; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) { char log[1024]; glGetProgramInfoLog(p, 1024, nullptr, log); std::cerr << "Program link error:\n" << log << "\n"; }
    return p;
}

void framebufferSizeChanged(GLFWwindow*, int, int);
void processInput(GLFWwindow*, float);
void mouseButtonCallback(GLFWwindow*, int, int, int);
void cursorPosCallback(GLFWwindow*, double, double);
void scrollCallback(GLFWwindow*, double, double);
void addOneParticle();
void removeLastParticle();
void initializeParticles();

static inline float dist2(const glm::vec3& a, const glm::vec3& b) {
    glm::vec3 d = b - a; return glm::dot(d, d);
}

static inline bool finite3(const glm::vec3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

static void rebuildTopologyFromPositions(const vector<glm::vec3>& P, int k, float maxDist, vector<glm::ivec3>& outTris)
{
    outTris.clear();
    const int N = (int)P.size();
    if (N < 3) return;
    const float maxD2 = maxDist * maxDist;

    struct Cand { float d2; int j; };
    vector<vector<int>> nbr(N);
    for (int i = 0; i < N; ++i) {
        vector<Cand> cand; cand.reserve(N - 1);
        for (int j = 0; j < N; ++j) if (j != i) {
            float d2 = dist2(P[i], P[j]);
            if (d2 <= maxD2) cand.push_back({ d2, j });
        }
        int take = std::min(k, (int)cand.size());
        if (take > 0 && take < (int)cand.size()) {
            std::nth_element(cand.begin(), cand.begin() + take, cand.end(),
                [](const Cand& a, const Cand& b) { return a.d2 < b.d2; });
        }
        nbr[i].reserve(take);
        for (int t = 0; t < take; ++t) nbr[i].push_back(cand[t].j);
    }

    std::unordered_set<unsigned long long> triSet;
    triSet.reserve(N * 8);
    auto keyOf = [](int a, int b, int c)->unsigned long long {
        unsigned int ia = a, ib = b, ic = c;
        return (unsigned long long(ia) << 42) | (unsigned long long(ib) << 21) | unsigned long long(ic);
        };

    for (int i = 0; i < N; ++i) {
        const auto& L = nbr[i];
        const int M = (int)L.size();
        for (int u = 0; u < M; ++u) {
            int j = L[u];
            for (int v = u + 1; v < M; ++v) {
                int k2 = L[v];
                if (dist2(P[j], P[k2]) > maxD2) continue;

                int a = i, b = j, c = k2;
                if (a > b) std::swap(a, b);
                if (b > c) std::swap(b, c);
                if (a > b) std::swap(a, b);
                auto key = keyOf(a, b, c);
                if (triSet.insert(key).second) outTris.emplace_back(a, b, c);
            }
        }
    }
}

static void accumulateVertexNormals(const vector<glm::vec3>& P, const vector<glm::ivec3>& T, vector<glm::vec3>& outVn)
{
    const int N = (int)P.size();
    outVn.assign(N, glm::vec3(0.0f));
    for (auto t : T) {
        glm::vec3 pa = P[t.x], pb = P[t.y], pc = P[t.z];
        glm::vec3 n = glm::cross(pb - pa, pc - pa); // 2*area * face normal
        if (!finite3(n)) continue;
        outVn[t.x] += n; outVn[t.y] += n; outVn[t.z] += n;
    }
    for (int i = 0; i < N; ++i) {
        float len2 = glm::dot(outVn[i], outVn[i]);
        if (len2 > 1e-12f && finite3(outVn[i])) outVn[i] *= glm::inversesqrt(len2);
        else outVn[i] = glm::vec3(0, 0, 1);
    }
}

static void syncFluidArrays()
{
    const int N = (int)particles.size();
    if ((int)fluidPos.size() != N) {
        vector<glm::vec3> newPos; newPos.reserve(N);
        for (auto* p : particles)
            newPos.emplace_back(float(p->center.x), float(p->center.y), float(p->center.z));
        fluidPos.swap(newPos);
        prevNorm.assign(N, glm::vec3(0.0f));
        fluidTris.clear(); // force rebuild
        topoTimer = TOPOLOGY_DT;
    }
}

static void updateFluidMesh(float dt)
{
    syncFluidArrays();
    const int N = (int)particles.size();
    if (N < 3) { trisVertexCount = 0; return; }

    for (int i = 0; i < N; ++i) {
        glm::vec3 target(float(particles[i]->center.x),
            float(particles[i]->center.y),
            float(particles[i]->center.z));
        fluidPos[i] = (1.0f - POS_ALPHA) * fluidPos[i] + POS_ALPHA * target;
    }

    topoTimer += dt;
    if (topoTimer >= TOPOLOGY_DT || fluidTris.empty()) {
        rebuildTopologyFromPositions(fluidPos, k_neighbors, max_edge_dist, fluidTris);
        topoTimer = 0.0f;
    }
    if (fluidTris.empty()) { trisVertexCount = 0; return; }

    vector<glm::vec3> Vn; Vn.reserve(N);
    accumulateVertexNormals(fluidPos, fluidTris, Vn);

    if ((int)prevNorm.size() != N) prevNorm.assign(N, glm::vec3(0));
    glm::vec3 L = glm::normalize(glm::vec3(0.45f, 0.75f, 0.5f));
    for (int i = 0; i < N; ++i) {
        glm::vec3 cur = Vn[i];
        glm::vec3 ema = (prevNorm[i] == glm::vec3(0)) ? cur
            : glm::normalize((1.0f - NORM_ALPHA) * prevNorm[i] + NORM_ALPHA * cur);
        if (glm::dot(ema, L) < 0.0f) ema = -ema;
        prevNorm[i] = ema;
        Vn[i] = ema;
    }

    vector<float> interleaved;
    interleaved.reserve(fluidTris.size() * 3 * 6);
    auto pushVN = [&](int idx) {
        const glm::vec3& p = fluidPos[idx];
        const glm::vec3& n = Vn[idx];
        interleaved.push_back(p.x);
        interleaved.push_back(p.y);
        interleaved.push_back(p.z);
        interleaved.push_back(n.x);
        interleaved.push_back(n.y);
        interleaved.push_back(n.z);
        };
    for (auto t : fluidTris) { pushVN(t.x); pushVN(t.y); pushVN(t.z); }
    trisVertexCount = static_cast<GLsizei>(interleaved.size() / 6);

    glBindVertexArray(trisVAO);
    glBindBuffer(GL_ARRAY_BUFFER, trisVBO);
    size_t bytes = interleaved.size() * sizeof(float);
    if (bytes > trisBytes) { glBufferData(GL_ARRAY_BUFFER, bytes, nullptr, GL_STREAM_DRAW); trisBytes = bytes; }
    glBufferSubData(GL_ARRAY_BUFFER, 0, bytes, interleaved.data());
    glBindVertexArray(0);
}

void initializeParticles() {
    for (auto* p : particles) delete p;
    particles.clear();

    if (num_particles <= 0) return;
    if (num_particles == 1) {
        particles.push_back(new Particle(0.0f, 0.0f, 0.0f, particle_radius, drawing_radius));
    }
    else {
        int dim = static_cast<int>(std::round(std::cbrt(static_cast<double>(num_particles))));
        if (dim < 1) dim = 1;
        float step = 2.0f * particle_radius;
        float start = -step * (dim - 1) * 0.5f;
        for (int x = 0; x < dim; ++x)
            for (int y = 0; y < dim; ++y)
                for (int z = 0; z < dim; ++z) {
                    if ((int)particles.size() >= num_particles) break;
                    particles.push_back(new Particle(start + x * step, start + y * step, start + z * step,
                        particle_radius, drawing_radius));
                }
    }
    syncFluidArrays();
}

int main() {
    if (!glfwInit()) { std::cerr << "Failed to initialize GLFW\n"; return -1; }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    const GLFWvidmode* vm = glfwGetVideoMode(glfwGetPrimaryMonitor());
    int H = vm ? vm->height : 900;

    GLFWwindow* window = glfwCreateWindow(H / 2, H / 2, "Waterbender", nullptr, nullptr);
    if (!window) { std::cerr << "Failed to create window\n"; glfwTerminate(); return -1; }
    glfwMakeContextCurrent(window);

    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetFramebufferSizeCallback(window, framebufferSizeChanged);

    if (glewInit() != GLEW_OK) { std::cerr << "Failed to init GLEW\n"; glfwTerminate(); return -1; }

    {
        GLuint vs = compile(GL_VERTEX_SHADER, VS_tris);
        GLuint fs = compile(GL_FRAGMENT_SHADER, FS_tris);
        progTris = link(vs, fs);
        glDeleteShader(vs); glDeleteShader(fs);
        uRot_tris = glGetUniformLocation(progTris, "uRot");
    }

    {
        GLuint vs = compile(GL_VERTEX_SHADER, VS_lines);
        GLuint fs = compile(GL_FRAGMENT_SHADER, FS_lines);
        progLines = link(vs, fs);
        glDeleteShader(vs); glDeleteShader(fs);
        uRot_lines = glGetUniformLocation(progLines, "uRot");
    }

    // triangles buffer: interleaved pos(3)+normal(3)
    glGenVertexArrays(1, &trisVAO);
    glGenBuffers(1, &trisVBO);
    glBindVertexArray(trisVAO);
    glBindBuffer(GL_ARRAY_BUFFER, trisVBO);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    glEnableVertexAttribArray(0); // pos
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); // normal
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glBindVertexArray(0);

    // container VBO
    glGenVertexArrays(1, &containerVAO);
    glGenBuffers(1, &containerVBO);
    glBindVertexArray(containerVAO);
    glBindBuffer(GL_ARRAY_BUFFER, containerVBO);
    glBufferData(GL_ARRAY_BUFFER, container.vertices.size() * sizeof(float),
        container.vertices.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);

    initializeParticles();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    float last_fps_update = 0.0f;
    int num_frames = 0;
    float fps = 0.0f;
    while (!glfwWindowShouldClose(window)) {
        float now = static_cast<float>(glfwGetTime());
        float dt = now - lastFrameStartTime;
        lastFrameStartTime = now;

        num_frames++;
        if (now - last_fps_update >= 1.0) {
            fps = (float)num_frames / (now - last_fps_update);
            last_fps_update = now;
            num_frames = 0;
            string title = "Waterbender- " + to_string((int)fps) + " FPS";
            glfwSetWindowTitle(window, title.c_str());
        }

        processInput(window, dt);

        // physics via Particle::updatePosition (particles move; we render only fluid)
        glm::mat4 invRotation = glm::transpose(rotationMatrix);
        glm::vec4 grav_vector = invRotation * glm::vec4(gravity, 1.0f);
        for (auto* p : particles) {
            p->updatePosition(&grav_vector, &container, &invRotation, &particles, dt);
        }

        // update the lagged fluid mesh and upload
        updateFluidMesh(dt);

        // render
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // fluid
        glUseProgram(progTris);
        glUniformMatrix4fv(uRot_tris, 1, GL_FALSE, glm::value_ptr(rotationMatrix));
        glBindVertexArray(trisVAO);
        glDisable(GL_CULL_FACE); // show both sides for translucent sheet
        glDepthMask(GL_FALSE);   // avoid popping with partial transparency
        glDrawArrays(GL_TRIANGLES, 0, trisVertexCount);
        glDepthMask(GL_TRUE);
        glBindVertexArray(0);

        // box
        glUseProgram(progLines);
        glUniformMatrix4fv(uRot_lines, 1, GL_FALSE, glm::value_ptr(rotationMatrix));
        glBindVertexArray(containerVAO);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(container.vertices.size() / 3));
        glBindVertexArray(0);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    for (auto* p : particles) delete p;
    glDeleteVertexArrays(1, &trisVAO);
    glDeleteBuffers(1, &trisVBO);
    glDeleteVertexArrays(1, &containerVAO);
    glDeleteBuffers(1, &containerVBO);
    glDeleteProgram(progTris);
    glDeleteProgram(progLines);

    glfwTerminate();
    return 0;
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) { gDragging = true;  glfwGetCursorPos(window, &gLastX, &gLastY); }
        if (action == GLFW_RELEASE) { gDragging = false; }
    }
}
void cursorPosCallback(GLFWwindow*, double xpos, double ypos) {
    if (!gDragging) return;
    double dx = xpos - gLastX, dy = ypos - gLastY; gLastX = xpos; gLastY = ypos;
    float yawDeg = gSensitivity * static_cast<float>(dx);
    float pitchDeg = gSensitivity * static_cast<float>(dy);
    rotationMatrix = glm::rotate(rotationMatrix, glm::radians(yawDeg), glm::vec3(0, 1, 0));
    rotationMatrix = glm::rotate(rotationMatrix, glm::radians(pitchDeg), glm::vec3(1, 0, 0));
}
void scrollCallback(GLFWwindow*, double, double yoffset) {
    if (yoffset > 0.0) addOneParticle();
    else if (yoffset < 0.0) removeLastParticle();
}
void addOneParticle() {
    auto rnd = []() { return static_cast<float>(rand()) / static_cast<float>(RAND_MAX); };
    particles.push_back(new Particle(rnd() * 0.2f - 0.1f,
        rnd() * 0.2f - 0.1f,
        rnd() * 0.2f - 0.1f,
        particle_radius,
        drawing_radius));
    num_particles = static_cast<int>(particles.size());
    syncFluidArrays();
}
void removeLastParticle() {
    if (particles.empty()) return;
    delete particles.back(); particles.pop_back();
    num_particles = static_cast<int>(particles.size());
    syncFluidArrays();
}
void framebufferSizeChanged(GLFWwindow*, int width, int height) { glViewport(0, 0, width, height); }
void processInput(GLFWwindow* window, float dt) {
    const float speed = 90.0f * dt;
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(-speed), glm::vec3(0, 1, 0));
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(speed), glm::vec3(0, 1, 0));
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(speed), glm::vec3(1, 0, 0));
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(-speed), glm::vec3(1, 0, 0));
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) rotationMatrix = glm::mat4(1.0f);
}
<<<<<<< HEAD
=======

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
    //cout << grav_vector.x << " " << grav_vector.y << " " << grav_vector.z << endl;
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

>>>>>>> 2579bf8221f6876b715f34a7a0db0678eec93bbe
