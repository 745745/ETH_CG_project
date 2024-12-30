#pragma once

#if !defined(__NORI_MEDIUM_H)
#define __NORI_MEDIUM_H

#include <nori/object.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class Medium;

struct MediumInteraction {
    Point3f p;          // transmit or scatter point
    float t;            // transmit distance
    Vector3f wo;        // outgoing direction
    Medium* medium=nullptr;
    bool scatter=false;
};

class Phase
{
public:
    float g;
    Phase(float g) {  this->g = g;};
    Phase() { g = 0; };
    float eval(const Vector3f& wi, const Vector3f& wo)
    {  
        float cos = wi.dot(wo);
        float denom = 1 + g * g + 2 * g * cos;
        return INV_FOURPI * (1 - g * g) / std::pow(denom, 1.5f);    
    }

    Vector3f sample(const Vector3f& wi, const Point2f& sample)
    {
        float cosTheta;
        if (abs(g) < 1e-3)
        {
            cosTheta = 1 - 2 * sample.x();
        }
        else
        {
            float t = (1 - g * g) / (1 - g + 2 * g * sample.x());
            cosTheta = (1 + g * g - t * t) / (2 * g);
        }
        
        float sinTheta = std::sqrt(std::max(0.f, 1 - cosTheta * cosTheta));
        float phi = 2 * M_PI * sample.y(); 

        Vector3f wo(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
        Frame frame(-wi); 
        wo= frame.toWorld(wo);
        return wo;
    }
};


class Medium : public NoriObject {
public:
    virtual ~Medium() = default;

    virtual Color3f transmittance(float dist) const = 0;

    virtual Color3f sample(const Ray3f& ray, Point2f& sample, MediumInteraction& mi,float& maxDist) const = 0;

    virtual EClassType getClassType() const override {
        return EMedium;
    }

    int type; //type=0: int medium  =1: ext medium
    Phase* phase;
    Color3f m_sigma_s;
    Color3f m_sigma_a;
    Color3f sigma_t;
};

NORI_NAMESPACE_END

#endif