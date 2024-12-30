#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
NORI_NAMESPACE_BEGIN


class Direct_MIS :public Integrator
{
public:
    Direct_MIS(const PropertyList& props) {}
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

        Color3f Lr_BRDF(0.f);
        Color3f Lr_EM(0.f);
        Color3f Le(0.f);
        float pdf_BRDF = 0.f;
        float pdf_EM = 0.f;
        auto bsdf = its.mesh->getBSDF();
        int n_lights = scene->getLights().size();
        
        EmitterQueryRecord lRec(its.p);
        Color3f brdf;
        if (its.mesh->isEmitter())
        {
            auto area_emitter = its.mesh->getEmitter();
            EmitterQueryRecord l(ray.o, its.p, its.shFrame.n);
            Le = area_emitter->eval(l);
        }
        

        //BRDF      
        BSDFQueryRecord bRec(its.toLocal(-ray.d), its.uv);

        brdf = bsdf->sample(bRec, sampler->next2D());
        pdf_BRDF = bsdf->pdf(bRec);
        Ray3f reflect_ray;
        reflect_ray=Ray3f(its.p, its.toWorld(bRec.wo), Epsilon, INFINITY);

        Intersection its2;
        if (scene->rayIntersect(reflect_ray, its2))
        {
            if (its2.mesh->isEmitter())
            {
                lRec = EmitterQueryRecord(its.p, its2.p, its2.shFrame.n);
                pdf_EM = its2.mesh->getEmitter()->pdf(lRec)/ n_lights;

                float w_BRDF = (pdf_BRDF + pdf_EM) > 0 ? pdf_BRDF / (pdf_BRDF + pdf_EM) : 0.f;
                //delta brdf
                if (bRec.measure == EDiscrete)
                {
                    Lr_BRDF = brdf * its2.mesh->getEmitter()->eval(lRec);
                }
                else Lr_BRDF = w_BRDF * brdf * its2.mesh->getEmitter()->eval(lRec);
            }
        }
        else
        {
            if (scene->getEnvLights() != nullptr)
            {
                EmitterQueryRecord lRec;
                lRec.wi = its.toWorld(bRec.wo);
                pdf_EM = scene->getEnvLights()->pdf(lRec);
                float w_BRDF = (pdf_BRDF + pdf_EM) > 0 ? pdf_BRDF / (pdf_BRDF + pdf_EM) : 0.f;
                Lr_BRDF = w_BRDF * brdf*scene->getEnvLights()->eval(lRec);
            }
        }

        //Emitter
        auto emitter = scene->getRandomEmitter(sampler->next1D());

        Color3f Li = emitter->sample(lRec, sampler->next2D());
        pdf_EM = emitter->pdf(lRec)/ n_lights;
        if (!scene->rayIntersect(lRec.shadowRay))
        {
            BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lRec.wi), ESolidAngle, its.uv);
            float cos = its.shFrame.cosTheta(its.toLocal(lRec.wi));                 
            pdf_BRDF = bsdf->pdf(bRec);
            brdf = bsdf->eval(bRec);
            float w_EM = (pdf_BRDF + pdf_EM) > 0 ? pdf_EM / (pdf_BRDF + pdf_EM) : 0.f;
            Lr_EM = w_EM * n_lights * brdf * Li * std::max(cos,0.f);
        }
 
        return Le +  Lr_BRDF +  Lr_EM;
    }


   
    std::string toString() const {
        return "Direct_MIS";
    }
};


NORI_REGISTER_CLASS(Direct_MIS, "direct_mis");
NORI_NAMESPACE_END