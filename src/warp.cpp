/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include<iostream>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = sqrtf(sample.x());
    float theta = 2 * M_PI * sample.y();
    return Point2f(r*cos(theta), r*sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return  (p.norm() <= 1 && p.norm() >= 0) ? 1.f / M_PI : 0.0f;
}

Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    float theta = acos(1 - (1 - cosThetaMax)*sample.x());
    float phi = 2 * M_PI * sample.y();
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {

    return (v.norm() - 1.f <= Epsilon && v.z() >= cosThetaMax) ? 1.f / (2 * M_PI * (1 - cosThetaMax)) : 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float theta=acos(1-2*sample.x());
    float phi = 2 * M_PI * sample.y();
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return (v.norm() - 1.f <= Epsilon) ? 1.f / (4 * M_PI) : 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float theta = acos(1 - sample.x());
    float phi = 2 * M_PI * sample.y();
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return (v.norm() - 1.f <= Epsilon && v.z()>=0) ? 1.f / (2 * M_PI) : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float theta = acos(sqrt(1-sample.x()));
    float phi = 2 * M_PI * sample.y();
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return (v.norm() - 1.f <= Epsilon && v.z() >= 0) ? v.z() / M_PI : 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float tan2theta = -alpha*alpha * log(1- sample.x());
    float costheta = sqrt(1 / (1 + tan2theta));
    float sintheta = sqrt(tan2theta/(1+tan2theta));

    float phi= 2 * M_PI * sample.y();
    float x = sintheta * cos(phi);
    float y = sintheta * sin(phi);
    float z = costheta;
    return Vector3f(x, y, z);

}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {

    if (m.norm() - 1.f <= Epsilon && m.z() >0)
    {
        float cos = m.z();
        float sin = sqrtf(1 - cos * cos);
        float tan = sin / cos;
        float alpha2 = pow(alpha, 2);

        float res = exp(-pow(tan, 2) / alpha2)  / (M_PI * alpha2 * pow(cos, 3));   
        return res;
    }
    return 0.0f;

}

Vector3f Warp::squareToUniformTriangle(const Point2f &sample) {
    float su1 = sqrtf(sample.x());
    float u = 1.f - su1, v = sample.y() * su1;
    return Vector3f(u,v,1.f-u-v);
}

NORI_NAMESPACE_END
