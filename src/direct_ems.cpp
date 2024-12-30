#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
NORI_NAMESPACE_BEGIN


class Direct_EMS :public Integrator
{
public:
    Direct_EMS(const PropertyList& props) {}
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {     
        
        Intersection its;
        if (!scene->rayIntersect(ray, its))
        {
            if (scene->getEnvLights() != nullptr)
            {
                EmitterQueryRecord lRec;
                lRec.wi = ray.d;

                return scene->getEnvLights()->eval(lRec);
            }
            else return 0.f;
        }

        Color3f Lr(0.f); 
        Color3f Le(0.f); 
        auto bsdf = its.mesh->getBSDF();
        auto emitter = scene->getRandomEmitter(sampler->next1D());

        if (its.mesh->isEmitter())
        {
            auto area_emitter = its.mesh->getEmitter();
            EmitterQueryRecord lRec(ray.o, its.p, its.shFrame.n);
            Le = area_emitter->eval(lRec);
        }

        EmitterQueryRecord lRec(its.p);
        Color3f Li = emitter->sample(lRec, sampler->next2D());
        if (scene->rayIntersect(lRec.shadowRay))
            Li = 0.f;

        BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lRec.wi), ESolidAngle, its.uv);
        float cos = its.shFrame.cosTheta(its.shFrame.toLocal(lRec.wi));


        Lr = scene->getLights().size() * bsdf->eval(bRec) * Li * std::max(cos, 0.f);
        
        return Lr + Le;
    }
    


    std::string toString() const {
        return "Direct_EMS";
    }
};


NORI_REGISTER_CLASS(Direct_EMS, "direct_ems");
NORI_NAMESPACE_END