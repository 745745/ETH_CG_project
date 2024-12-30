#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
NORI_NAMESPACE_BEGIN


class Path_MATS :public Integrator
{
public:
    Path_MATS(const PropertyList& props) {}
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {                
        Intersection its;
        int bounces = 0;
        Color3f L(0.f);
        Color3f t(1.f);
        Ray3f contin_ray = ray;
        while (true)
        {
            if (!scene->rayIntersect(contin_ray, its))
            {
                return L;
            }

            if (its.mesh->isEmitter())
            {
                auto area_emitter = its.mesh->getEmitter();
                EmitterQueryRecord lRec(contin_ray.o, its.p, its.shFrame.n);
                L += t*area_emitter->eval(lRec);
            }
           
            if (bounces >= 3)
            {
                float q = std::min(t.maxCoeff(), 0.99f);
                if (sampler->next1D() > q)
                    break;
                else t /= q;
            }

            BSDFQueryRecord bRec(its.shFrame.toLocal(-contin_ray.d), its.uv);
            Color3f brdf;
            if (its.mesh->getBSDF() != nullptr)
                brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            else
                brdf = 0;
            t *= brdf;   
            bounces++;
            contin_ray = Ray3f(its.p, its.toWorld(bRec.wo),Epsilon,INFINITY);
        }
        

        return L;
    }



    std::string toString() const {
        return "Path_MATS";
    }
};


NORI_REGISTER_CLASS(Path_MATS, "path_mats");
NORI_NAMESPACE_END