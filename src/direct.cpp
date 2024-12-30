#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
NORI_NAMESPACE_BEGIN


class Direct :public Integrator
{
public:
    Direct(const PropertyList& props) {}
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Color3f res(0.f);
        Vector2f sample;

        for (auto light : scene->getLights())
        {
            EmitterQueryRecord lrec;
            lrec.ref = its.p;
            Color3f res_p = light->sample(lrec, sample);
            if (scene->rayIntersect(lrec.shadowRay))
                continue;

            float cosine = its.shFrame.cosTheta(its.shFrame.toLocal(lrec.wi));
            auto bsdf = its.mesh->getBSDF();
            BSDFQueryRecord brec(its.shFrame.toLocal(lrec.wi), its.shFrame.toLocal(-ray.d), ESolidAngle,its.uv);
            brec.uv = its.uv;
            Color3f t = bsdf->eval(brec);
            res += res_p * t * cosine;
        }
        return res;
    }

    std::string toString() const {
        return "Direct[]";
    }
};


NORI_REGISTER_CLASS(Direct, "direct");
NORI_NAMESPACE_END