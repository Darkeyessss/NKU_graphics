#include "L2d.h"
#include <glm/ext/scalar_constants.hpp>
#include <iostream>
using namespace std;

spline::spline(float h) {
    mH = h;
    mH2 = h * h;
    sigma = 40.0 / (7.0 * glm::pi<float>() * mH2);

    mBufferSize = glm::uvec2(128, 128);
    mGradBuffer = vector<vector<glm::vec2>>(mBufferSize.x, vector<glm::vec2>(mBufferSize.y));
    mValueBuffer = vector<float>(mBufferSize.x);

    for (int i = 0; i < mBufferSize.x; i++) {
        for (int j = 0; j < mBufferSize.y; j++) {
            float x = ((float)i + 0.5f) * mH / mBufferSize.x;
            float y = ((float)j + 0.5f) * mH / mBufferSize.y;
            glm::vec2 radius(x, y);
            mGradBuffer[i][j] = CalGrad(radius);
        }
    }

    for (int i = 0; i < mBufferSize.x; i++) {
        float distance = ((float)i + 0.5f) * mH / mBufferSize.x;
        mValueBuffer[i] = CalValue(distance);
    }
}

spline::~spline() {}

float spline::Value(float distance)
{
    float res = 0;
    int i = (abs(distance) * mBufferSize.x / mH);
    if (i >= mBufferSize.x) {
        return res;
    }
    res = mValueBuffer[i];
    return res;
}

glm::vec2 spline::Grad(glm::vec2 radius) {
    glm::vec2 res(0.0f, 0.0f);

    int i = (abs(radius.x) * mBufferSize.x / mH);
    int j = (abs(radius.y) * mBufferSize.x / mH);

    if (i >= mBufferSize.x || j >= mBufferSize.y) {
        return res;
    }

    res = mGradBuffer[i][j];

    if (radius.x < 0) {
        res.x = -res.x;
    }
    if (radius.y < 0) {
        res.y = -res.y;
    }

    return res;
}

glm::vec2 spline::CalGrad(glm::vec2 radius) {
    glm::vec2 res(0.0f, 0.0f);
    float r = glm::length(radius);
    float distance = r;
    if (distance < 1e-5) {
        return res;
    }

    float q = distance / mH;
    glm::vec2 qGrad = radius / (mH * distance);

    if (q < 0.5f) {
        res = 6.0f * (3.0f * q * q - 2.0f * q) * sigma * qGrad;
        return res;
    }
    else if (q >= 0.5 && q < 1.0f) {
        res = -6.0f * powf(1.0f - q, 2) * sigma * qGrad;
        return res;
    }
    return res;
}

float spline::CalValue(float dist) {
    float absDist = abs(dist);
    float q = absDist / mH;
    float q2 = q * q;
    float q3 = q * q2;
    float value = 0.0f;

    if (q < 0.5f) {
        value = 6.0f * (q3 - q2) + 1.0f;
        value *= sigma;
    }
    else if (q >= 0.5f && q < 1.0f) {
        value = 1.0f - q;
        value = pow(value, 3) * 2.0f;
        value *= sigma;
    }

    return value;
}
