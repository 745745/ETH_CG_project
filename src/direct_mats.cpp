#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
NORI_NAMESPACE_BEGIN


class Direct_MATS :public Integrator
{
public:
    Direct_MATS(const PropertyList& props) {}
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {

        Intersection its;
        if (!scene->rayIntersect(ray, its))
        {
            return Color3f(0.f);
        }

        Color3f Lr(0.f);
        Color3f Le(0.f);
        auto bsdf = its.mesh->getBSDF();

        if (its.mesh->isEmitter())
        {
            auto area_emitter = its.mesh->getEmitter();
            EmitterQueryRecord lRec(ray.o, its.p, its.shFrame.n);
            Le = area_emitter->eval(lRec);
        }

        BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d),its.uv);
        Color3f brdf = bsdf->sample(bRec, sampler->next2D());
        Intersection its2;
        Ray3f reflect_ray(its.p, its.shFrame.toWorld(bRec.wo));
        if (scene->rayIntersect(reflect_ray, its2))
        {
            if (its2.mesh->isEmitter())
            {
                EmitterQueryRecord lRec = EmitterQueryRecord(reflect_ray.o, its2.p, its2.shFrame.n);
                Lr = brdf *its2.mesh->getEmitter()->eval(lRec);
            }
        }

        return Le + Lr;

    }



    std::string toString() const {
        return "Direct_EMS";
    }
};


NORI_REGISTER_CLASS(Direct_MATS, "direct_mats");
NORI_NAMESPACE_END