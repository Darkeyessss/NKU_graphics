#pragma once
#ifndef W_CUBE_SPLINE_H
#define W_CUBE_SPLINE_H

#include <glm/glm.hpp>
#include <vector>


class spline {
public:
    spline() = delete;
    explicit spline(float h);
    ~spline();

    float Value(float distance);

    glm::vec2 Grad(glm::vec2 radius);

private:
    float CalValue(float distance);

    glm::vec2 CalGrad(glm::vec2 radius);

    float mH;
    float mH2;
    float sigma;
    glm::uvec2 mBufferSize;
    std::vector<std::vector<glm::vec2>> mGradBuffer;
    std::vector<float> mValueBuffer;
};
#endif 