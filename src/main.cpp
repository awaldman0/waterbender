// main.cpp - Fluid mesh (lagged kNN triangles) + box + MODERATE bloom with -30% intensity
// Adds toggle: press 'P' to switch between shaded mesh and raw particle points.
// GL 3.3 core (no immediate mode). Requires: GLEW, GLFW, GLM, your container.h, particle.h.

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
#include <string>

#include "container.h"
#include "particle.h"

using std::vector;
using std::string;

// ---------------- Parameters ----------------
int   num_particles = 512;
float particle_radius = 0.04f;
float drawing_radius = 0.02f;
glm::vec3 gravity = glm::vec3(0.0f, -.000005f, 0.0f);

int   k_neighbors = 4;
float max_edge_dist = 0.35f;

const float POS_ALPHA = 0.12f;   // position EMA
const float NORM_ALPHA = 0.20f;   // normal EMA
const float TOPOLOGY_DT = 0.08f;   // topology rebuild cadence (s)

// Bloom controls (moderate; 30% reduction applied)
float BRIGHT_THRESHOLD = 1.10f; // slight raise to preserve detail
float BLOOM_STRENGTH = 0.63f; // 30% lower than 0.9
float COMPOSITE_EXPOSURE = 0.84f; // 30% lower than 1.2
int   BLUR_PASSES = 6;

// ---------------- Input ----------------
bool   gDragging = false;
double gLastX = 0.0, gLastY = 0.0;
float  gSensitivity = 0.25f;
bool   gShowParticles = false; // toggle target

// ---------------- Scene ----------------
glm::mat4 rotationMatrix = glm::mat4(1.0f);
Container container;
vector<Particle*> particles;

// ---------------- Fluid (lagged) state ----------------
static inline bool finite3(const glm::vec3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}
vector<glm::vec3> fluidPos;   // lagged positions
vector<glm::vec3> prevNorm;   // EMA normals
vector<glm::ivec3> fluidTris; // topology
float topoTimer = 0.0f;

// ---------------- GL objects (scene) ----------------
GLuint progTris = 0;            // shaded translucent mesh
GLint  uRot_tris = -1;

GLuint progLines = 0;            // wireframe box
GLint  uRot_lines = -1;

GLuint progPoints = 0;            // particle points
GLint  uRot_points = -1;

GLuint trisVAO = 0, trisVBO = 0; // interleaved pos(3)+normal(3)
size_t trisBytes = 0;
GLsizei trisVertexCount = 0;

GLuint pointsVAO = 0, pointsVBO = 0; // particle positions
size_t pointsBytes = 0;
GLsizei pointsCount = 0;

GLuint containerVAO = 0, containerVBO = 0;

// ---------------- HDR + Bloom ----------------
GLuint hdrFBO = 0, colorTex = 0, depthRBO = 0;
GLuint pingFBO[2] = { 0,0 }, pingTex[2] = { 0,0 };
GLuint quadVAO = 0, quadVBO = 0;

GLuint progExtract = 0, progBlur = 0, progComposite = 0;
GLint  uExtract_threshold = -1, uBlur_horizontal = -1, uComp_bloomStrength = -1, uComp_exposure = -1;

int screenW = 960, screenH = 960;
int bloomW = 480, bloomH = 480;     // half-res blur

// ---------------- Time ----------------
float lastFrameStartTime = 0.0f;

// ---------------- Shaders ----------------
static const char* VS_tris = R"(
#version 330 core
layout (location=0) in vec3 aPos;
layout (location=1) in vec3 aNor;
uniform mat4 uRot;
out vec3 vN;
void main() {
    gl_Position = uRot * vec4(aPos, 1.0);
    vN = mat3(uRot) * aNor; // rotate normal only
}
)";

// Moderate emissive water to keep geometry readable (alpha 0.6)
static const char* FS_tris = R"(
#version 330 core
in vec3 vN;
out vec4 FragColor;
void main() {
    vec3 N = normalize(vN);
    vec3 L = normalize(vec3(0.45, 0.75, 0.5));
    float lambert = max(dot(N, L), 0.0);

    vec3 base = vec3(0.2, 0.6, 1.0);

    vec3 shaded   = base * (0.25 + 1.2 * lambert);
    vec3 emissive = base * 0.84;         // 30% lower than 1.2

    vec3 col = shaded + emissive;
    FragColor = vec4(col, 0.60);
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

// Particle points (GL_POINTS) - core profile compliant
static const char* VS_points = R"(
#version 330 core
layout (location=0) in vec3 aPos;
uniform mat4 uRot;
void main(){
    gl_Position = uRot * vec4(aPos, 1.0);
    gl_PointSize = 5.0; // size in pixels
}
)";
static const char* FS_points = R"(
#version 330 core
out vec4 FragColor;
void main(){
    // simple round-ish point: discard corners for a circle look
    vec2 c = (gl_PointCoord - vec2(0.5)) * 2.0;
    if (dot(c,c) > 1.0) discard;
    FragColor = vec4(0.9, 0.95, 1.0, 1.0);
}
)";

// Fullscreen quad
static const char* VS_quad = R"(
#version 330 core
layout(location=0) in vec2 aPos;
layout(location=1) in vec2 aUV;
out vec2 vUV;
void main(){ vUV=aUV; gl_Position=vec4(aPos,0.0,1.0); }
)";

// Bright extract with moderate gain
static const char* FS_extract = R"(
#version 330 core
in vec2 vUV;
out vec4 FragColor;
uniform sampler2D uScene;
uniform float uThreshold; // ~1.1
void main(){
    vec3 c = texture(uScene, vUV).rgb;
    vec3 br = max(c - vec3(uThreshold), vec3(0.0)) * 1.2;
    FragColor = vec4(br, 1.0);
}
)";

// Gaussian blur (separable)
static const char* FS_blur = R"(
#version 330 core
in vec2 vUV;
out vec4 FragColor;
uniform sampler2D uImage;
uniform bool uHorizontal;
void main(){
    vec2 texel = 1.0 / vec2(textureSize(uImage, 0));
    float w[5] = float[5](0.227027, 0.1945946, 0.1216216, 0.054054, 0.016216);
    vec3 sum = texture(uImage, vUV).rgb * w[0];
    for(int i=1;i<5;++i){
        vec2 off = (uHorizontal ? vec2(texel.x*i,0) : vec2(0,texel.y*i));
        sum += texture(uImage, vUV + off).rgb * w[i];
        sum += texture(uImage, vUV - off).rgb * w[i];
    }
    FragColor = vec4(sum, 1.0);
}
)";

// Composite + ACES tonemap + exposure
static const char* FS_composite = R"(
#version 330 core
in vec2 vUV;
out vec4 FragColor;
uniform sampler2D uScene;
uniform sampler2D uBloom;
uniform float uBloomStrength; // ~0.63
uniform float uExposure;      // ~0.84

vec3 ACESFilm(vec3 x){
    const float a=2.51, b=0.03, c=2.43, d=0.59, e=0.14;
    return clamp((x*(a*x + b)) / (x*(c*x + d) + e), 0.0, 1.0);
}

void main(){
    vec3 scene = texture(uScene, vUV).rgb;
    vec3 bloom = texture(uBloom, vUV).rgb;

    vec3 hdr = scene + uBloomStrength * bloom;
    hdr *= uExposure;

    vec3 ldr = ACESFilm(hdr);
    FragColor = vec4(ldr, 1.0);
}
)";

// ---------------- GL utilities ----------------
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
template<typename F>
static GLuint makeProgram(const char* vsrc, const char* fsrc, F&& afterLink) {
    GLuint vs = compile(GL_VERTEX_SHADER, vsrc);
    GLuint fs = compile(GL_FRAGMENT_SHADER, fsrc);
    GLuint p = link(vs, fs);
    glDeleteShader(vs); glDeleteShader(fs);
    afterLink(p);
    return p;
}

// ---------------- Decls ----------------
void framebufferSizeChanged(GLFWwindow*, int, int);
void processInput(GLFWwindow*, float);
void mouseButtonCallback(GLFWwindow*, int, int, int);
void cursorPosCallback(GLFWwindow*, double, double);
void scrollCallback(GLFWwindow*, double, double);
void addOneParticle();
void removeLastParticle();
void initializeParticles();

// ---------------- Helpers ----------------
static inline float dist2(const glm::vec3& a, const glm::vec3& b) {
    glm::vec3 d = b - a; return glm::dot(d, d);
}

// kNN triangles (lagged P)
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

// area-weighted vertex normals
static void accumulateVertexNormals(const vector<glm::vec3>& P, const vector<glm::ivec3>& T, vector<glm::vec3>& outVn)
{
    const int N = (int)P.size();
    outVn.assign(N, glm::vec3(0.0f));
    for (auto t : T) {
        glm::vec3 pa = P[t.x], pb = P[t.y], pc = P[t.z];
        glm::vec3 n = glm::cross(pb - pa, pc - pa); // 2*area * unit normal
        if (!finite3(n)) continue;
        outVn[t.x] += n; outVn[t.y] += n; outVn[t.z] += n;
    }
    for (int i = 0; i < N; ++i) {
        float len2 = glm::dot(outVn[i], outVn[i]);
        if (len2 > 1e-12f && finite3(outVn[i])) outVn[i] *= glm::inversesqrt(len2);
        else outVn[i] = glm::vec3(0, 0, 1);
    }
}

// sync lagged arrays to particle count
static void syncFluidArrays()
{
    const int N = (int)particles.size();
    if ((int)fluidPos.size() != N) {
        vector<glm::vec3> newPos; newPos.reserve(N);
        for (auto* p : particles)
            newPos.emplace_back(float(p->center.x), float(p->center.y), float(p->center.z));
        fluidPos.swap(newPos);
        prevNorm.assign(N, glm::vec3(0.0f));
        fluidTris.clear();
        topoTimer = TOPOLOGY_DT; // force immediate rebuild
    }
}

// update lagged mesh + upload VBO
static void updateFluidMesh(float dt)
{
    syncFluidArrays();
    const int N = (int)particles.size();
    if (N < 3) { trisVertexCount = 0; return; }

    // position EMA
    for (int i = 0; i < N; ++i) {
        glm::vec3 target(float(particles[i]->center.x),
            float(particles[i]->center.y),
            float(particles[i]->center.z));
        fluidPos[i] = (1.0f - POS_ALPHA) * fluidPos[i] + POS_ALPHA * target;
    }

    // topology throttle
    topoTimer += dt;
    if (topoTimer >= TOPOLOGY_DT || fluidTris.empty()) {
        rebuildTopologyFromPositions(fluidPos, k_neighbors, max_edge_dist, fluidTris);
        topoTimer = 0.0f;
    }
    if (fluidTris.empty()) { trisVertexCount = 0; return; }

    // normals
    vector<glm::vec3> Vn; Vn.reserve(N);
    accumulateVertexNormals(fluidPos, fluidTris, Vn);

    // EMA normals + hemisphere
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

    // interleaved buffer
    vector<float> interleaved;
    interleaved.reserve(fluidTris.size() * 3 * 6);
    auto pushVN = [&](int idx) {
        const glm::vec3& p = fluidPos[idx];
        const glm::vec3& n = Vn[idx];
        interleaved.push_back(p.x); interleaved.push_back(p.y); interleaved.push_back(p.z);
        interleaved.push_back(n.x); interleaved.push_back(n.y); interleaved.push_back(n.z);
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

// upload particle positions for point rendering
static void updatePointsVBO()
{
    vector<glm::vec3> P;
    P.reserve(particles.size());
    for (auto* p : particles) P.emplace_back(float(p->center.x), float(p->center.y), float(p->center.z));
    pointsCount = static_cast<GLsizei>(P.size());
    glBindVertexArray(pointsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
    size_t bytes = P.size() * sizeof(glm::vec3);
    if (bytes > pointsBytes) { glBufferData(GL_ARRAY_BUFFER, bytes, nullptr, GL_STREAM_DRAW); pointsBytes = bytes; }
    glBufferSubData(GL_ARRAY_BUFFER, 0, bytes, P.data());
    glBindVertexArray(0);
}

// ---------------- HDR/Bloom helpers ----------------
void makeHDRFBO(int w, int h) {
    if (hdrFBO) { glDeleteFramebuffers(1, &hdrFBO); glDeleteTextures(1, &colorTex); glDeleteRenderbuffers(1, &depthRBO); }
    glGenFramebuffers(1, &hdrFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, hdrFBO);

    glGenTextures(1, &colorTex);
    glBindTexture(GL_TEXTURE_2D, colorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorTex, 0);

    glGenRenderbuffers(1, &depthRBO);
    glBindRenderbuffer(GL_RENDERBUFFER, depthRBO);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, w, h);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, depthRBO);

    GLenum bufs[1] = { GL_COLOR_ATTACHMENT0 };
    glDrawBuffers(1, bufs);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) std::cerr << "HDR FBO incomplete\n";
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void makePingPong(int w, int h) {
    if (pingFBO[0]) { glDeleteFramebuffers(2, pingFBO); glDeleteTextures(2, pingTex); pingFBO[0] = pingFBO[1] = 0; }
    glGenFramebuffers(2, pingFBO);
    glGenTextures(2, pingTex);
    for (int i = 0; i < 2; ++i) {
        glBindFramebuffer(GL_FRAMEBUFFER, pingFBO[i]);
        glBindTexture(GL_TEXTURE_2D, pingTex[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pingTex[i], 0);
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) std::cerr << "Ping FBO incomplete\n";
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void makeQuad() {
    if (quadVAO) return;
    float verts[] = {
        // pos   // uv
        -1, -1,  0,0,
         1, -1,  1,0,
        -1,  1,  0,1,
         1,  1,  1,1
    };
    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &quadVBO);
    glBindVertexArray(quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);
}

// ---------------- Particles init ----------------
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

// ---------------- Callbacks ----------------
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
void framebufferSizeChanged(GLFWwindow*, int w, int h) {
    screenW = w; screenH = h;
    bloomW = std::max(1, w / 2);
    bloomH = std::max(1, h / 2);
    glViewport(0, 0, w, h);
    makeHDRFBO(w, h);
    makePingPong(bloomW, bloomH);
}

// ---------------- Input per-frame ----------------
void processInput(GLFWwindow* window, float dt)
{
    const float speed = 90.0f * dt;
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(-speed), glm::vec3(0, 1, 0));
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(speed), glm::vec3(0, 1, 0));
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(speed), glm::vec3(1, 0, 0));
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) rotationMatrix = glm::rotate(rotationMatrix, glm::radians(-speed), glm::vec3(1, 0, 0));

    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
      rotationMatrix = glm::mat4(1.0f);
      container.resetSize();
    }

    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) container.setWidth(container.width + 0.0005f);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) container.setWidth(container.width - 0.0005f);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) container.setHeight(container.height + 0.0005f);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) container.setHeight(container.height - 0.0005f);
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) container.setLength(container.length + 0.0005f);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) container.setLength(container.length - 0.0005f);

    // 'P' edge toggle for particle view
    static bool pPrev = false;
    bool pNow = glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS;
    if (pNow && !pPrev) gShowParticles = !gShowParticles;
    pPrev = pNow;
}

// ---------------- Main ----------------
int main() {
    if (!glfwInit()) { std::cerr << "Failed to initialize GLFW\n"; return -1; }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWmonitor* mon = glfwGetPrimaryMonitor();
    const GLFWvidmode* vm = glfwGetVideoMode(mon);
    int H = vm ? vm->height : 900;
    screenW = screenH = H / 2;
    bloomW = screenW / 2; bloomH = screenH / 2;

    GLFWwindow* window = glfwCreateWindow(screenW, screenH, "Waterbender", nullptr, nullptr);
    if (!window) { std::cerr << "Failed to create window\n"; glfwTerminate(); return -1; }
    glfwMakeContextCurrent(window);

    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetFramebufferSizeCallback(window, framebufferSizeChanged);

    if (glewInit() != GLEW_OK) { std::cerr << "Failed to init GLEW\n"; glfwTerminate(); return -1; }

    glEnable(GL_PROGRAM_POINT_SIZE);

    // Programs
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
    {
        GLuint vs = compile(GL_VERTEX_SHADER, VS_points);
        GLuint fs = compile(GL_FRAGMENT_SHADER, FS_points);
        progPoints = link(vs, fs);
        glDeleteShader(vs); glDeleteShader(fs);
        uRot_points = glGetUniformLocation(progPoints, "uRot");
    }

    // Post-process programs
    progExtract = makeProgram(VS_quad, FS_extract, [](GLuint) {});
    progBlur = makeProgram(VS_quad, FS_blur, [](GLuint) {});
    progComposite = makeProgram(VS_quad, FS_composite, [](GLuint) {});

    // Sampler bindings and uniforms
    glUseProgram(progExtract);
    glUniform1i(glGetUniformLocation(progExtract, "uScene"), 0);
    uExtract_threshold = glGetUniformLocation(progExtract, "uThreshold");

    glUseProgram(progBlur);
    glUniform1i(glGetUniformLocation(progBlur, "uImage"), 0);
    uBlur_horizontal = glGetUniformLocation(progBlur, "uHorizontal");

    glUseProgram(progComposite);
    glUniform1i(glGetUniformLocation(progComposite, "uScene"), 0);
    glUniform1i(glGetUniformLocation(progComposite, "uBloom"), 1);
    uComp_bloomStrength = glGetUniformLocation(progComposite, "uBloomStrength");
    uComp_exposure = glGetUniformLocation(progComposite, "uExposure");

    // Mesh VAO/VBO
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

    // Points VAO/VBO
    glGenVertexArrays(1, &pointsVAO);
    glGenBuffers(1, &pointsVBO);
    glBindVertexArray(pointsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);

    // Box VAO/VBO
    glGenVertexArrays(1, &containerVAO);
    glGenBuffers(1, &containerVBO);
    glBindVertexArray(containerVAO);
    glBindBuffer(GL_ARRAY_BUFFER, containerVBO);
    glBufferData(GL_ARRAY_BUFFER, container.vertices.size() * sizeof(float),
        container.vertices.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);

    // HDR + Bloom targets
    makeHDRFBO(screenW, screenH);
    makePingPong(bloomW, bloomH);
    makeQuad();

    initializeParticles();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // default

    float last_fps_update = 0.0f; int num_frames = 0; float fps = 0.0f;

    while (!glfwWindowShouldClose(window)) {
        float now = static_cast<float>(glfwGetTime());
        float dt = now - lastFrameStartTime;
        lastFrameStartTime = now;

        num_frames++;
        if (now - last_fps_update >= 1.0f) {
            fps = (float)num_frames / (now - last_fps_update);
            last_fps_update = now;
            num_frames = 0;
            string title = "Waterbender " + std::to_string((int)fps) + " FPS";
            glfwSetWindowTitle(window, title.c_str());
        }

        processInput(window, dt);

        // Physics
        glm::mat4 invRotation = glm::transpose(rotationMatrix);
        glm::vec4 grav_vector = invRotation * glm::vec4(gravity, 1.0f);

        for (auto* p : particles) {
            p->updatePosition(&grav_vector, &container, &invRotation, &particles, dt);
        }

        // Update visuals
        updateFluidMesh(dt);
        updatePointsVBO();

        // ---------- Pass 1: render scene into HDR FBO ----------
        glBindFramebuffer(GL_FRAMEBUFFER, hdrFBO);
        glViewport(0, 0, screenW, screenH);
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (gShowParticles) {
            // draw particles as points (opaque so they read clearly)
            glUseProgram(progPoints);
            glUniformMatrix4fv(uRot_points, 1, GL_FALSE, glm::value_ptr(rotationMatrix));
            glBindVertexArray(pointsVAO);
            glDisable(GL_CULL_FACE);
            glDepthMask(GL_TRUE);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glDrawArrays(GL_POINTS, 0, pointsCount);
            glBindVertexArray(0);
        }
        else {
            // draw shaded mesh (alpha to preserve geometry)
            glUseProgram(progTris);
            glUniformMatrix4fv(uRot_tris, 1, GL_FALSE, glm::value_ptr(rotationMatrix));
            glBindVertexArray(trisVAO);
            glDisable(GL_CULL_FACE);
            glDepthMask(GL_FALSE);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glDrawArrays(GL_TRIANGLES, 0, trisVertexCount);
            glDepthMask(GL_TRUE);
            glBindVertexArray(0);
        }

        // box
        glUseProgram(progLines);
        glUniformMatrix4fv(uRot_lines, 1, GL_FALSE, glm::value_ptr(rotationMatrix));
        glBindVertexArray(containerVAO);
        glBindBuffer(GL_ARRAY_BUFFER, containerVBO);
        glBufferData(GL_ARRAY_BUFFER, container.vertices.size() * sizeof(float),
          container.vertices.data(), GL_STATIC_DRAW);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(container.vertices.size() / 3));
        glBindVertexArray(0);

        // ---------- Pass 2: bright extract (half res) ----------
        glDisable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_FRAMEBUFFER, pingFBO[0]);
        glViewport(0, 0, bloomW, bloomH);
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(progExtract);
        glUniform1f(uExtract_threshold, BRIGHT_THRESHOLD);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, colorTex);
        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

        // ---------- Pass 3: Gaussian blur ping-pong (half res) ----------
        bool horizontal = true;
        for (int i = 0; i < BLUR_PASSES; ++i) {
            glBindFramebuffer(GL_FRAMEBUFFER, pingFBO[horizontal ? 1 : 0]);
            glViewport(0, 0, bloomW, bloomH);
            glUseProgram(progBlur);
            glUniform1i(uBlur_horizontal, horizontal ? 1 : 0);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, pingTex[horizontal ? 0 : 1]);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
            horizontal = !horizontal;
        }
        GLuint bloomTex = pingTex[horizontal ? 0 : 1];

        // ---------- Pass 4: composite to default framebuffer ----------
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, screenW, screenH);
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(progComposite);
        glUniform1f(uComp_bloomStrength, BLOOM_STRENGTH);
        glUniform1f(uComp_exposure, COMPOSITE_EXPOSURE);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, colorTex);   // scene HDR
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, bloomTex);   // blurred brights
        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glEnable(GL_DEPTH_TEST); // restore

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    for (auto* p : particles) delete p;

    glDeleteVertexArrays(1, &trisVAO);
    glDeleteBuffers(1, &trisVBO);
    glDeleteVertexArrays(1, &pointsVAO);
    glDeleteBuffers(1, &pointsVBO);
    glDeleteVertexArrays(1, &containerVAO);
    glDeleteBuffers(1, &containerVBO);

    if (hdrFBO) { glDeleteFramebuffers(1, &hdrFBO); }
    if (colorTex) { glDeleteTextures(1, &colorTex); }
    if (depthRBO) { glDeleteRenderbuffers(1, &depthRBO); }
    if (pingFBO[0]) glDeleteFramebuffers(2, pingFBO);
    if (pingTex[0]) glDeleteTextures(2, pingTex);
    if (quadVAO) { glDeleteVertexArrays(1, &quadVAO); glDeleteBuffers(1, &quadVBO); }

    glDeleteProgram(progTris);
    glDeleteProgram(progLines);
    glDeleteProgram(progPoints);
    glDeleteProgram(progExtract);
    glDeleteProgram(progBlur);
    glDeleteProgram(progComposite);

    glfwTerminate();
    return 0;
}
