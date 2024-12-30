#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>

NORI_NAMESPACE_BEGIN

class AverageVisibility : public Integrator {
public:
    float length;

    AverageVisibility(const PropertyList& props) {
        length = (float)props.getFloat("length");
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(1.0f);

        Vector3f r = Warp::sampleUniformHemisphere(sampler, its.shFrame.n);
        Ray3f shadow(its.p, r);
        shadow.mint = Epsilon;
        shadow.maxt = length;

        Intersection its2;
        if (!scene->rayIntersect(shadow, its2))
            return Color3f(1.0f);
        else
            return Color3f(0.0f);
    }

    std::string toString() const {
        return "NormalIntegrator[]";
    }
};

NORI_REGISTER_CLASS(AverageVisibility, "av");
NORI_NAMESPACE_END