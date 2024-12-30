#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/bsdf.h>
#include<nori/warp.h>

NORI_NAMESPACE_BEGIN

class MaterialIntegrator: public Integrator {
public:

    MaterialIntegrator(const PropertyList& props) {
        
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        BSDFQueryRecord BRec(its.toLocal(-ray.d), Vector3f(0.f, 0.f, 1.f), ESolidAngle, its.uv);

        const BSDF* bsdf = its.mesh->getBSDF();

        return bsdf->eval(BRec);
    }   

    std::string toString() const {
        return "NormalIntegrator[]";
    }
};

NORI_REGISTER_CLASS(MaterialIntegrator, "material");
NORI_NAMESPACE_END