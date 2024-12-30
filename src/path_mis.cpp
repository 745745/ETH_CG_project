#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
NORI_NAMESPACE_BEGIN


class Path_MIS :public Integrator
{
public:
    Path_MIS(const PropertyList& props) {}
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        int bounces = 0;
        Color3f L(0.f);
        Color3f t(1.f);
        Ray3f contin_ray = ray;
        int n_lights = scene->getLights().size();
        float pdf_em = 0, pdf_mat = 0;
        float w_mat = 1, w_em = 1;

        if (!scene->rayIntersect(contin_ray, its))
        {
            if (scene->getEnvLights() != nullptr)
            {
                EmitterQueryRecord lRec;
                lRec.wi = ray.d;
                return scene->getEnvLights()->eval(lRec);
            }
            else return L;
        }

        if (its.mesh->isEmitter())
        {
            EmitterQueryRecord lRec(contin_ray.o, its.p, its.shFrame.n);
            L += t * its.mesh->getEmitter()->eval(lRec);
        }

        while (true)
        {

            if (bounces >= 3)
            {
                float q = std::min(t.maxCoeff(), 0.99f);
                if (sampler->next1D() > q)
                    break;
                else t /= q;
            }

            //emitter
            auto bsdf = its.mesh->getBSDF();
            auto emitter = scene->getRandomEmitter(sampler->next1D());

            EmitterQueryRecord lRec(its.p);
            Color3f Li = emitter->sample(lRec, sampler->next2D());
            pdf_em = emitter->pdf(lRec) / n_lights;
            if (!scene->rayIntersect(lRec.shadowRay))
            {
                BSDFQueryRecord bRec(its.toLocal(-contin_ray.d), its.toLocal(lRec.wi), ESolidAngle, its.uv);

                float cos = its.shFrame.cosTheta(its.toLocal(lRec.wi));
                pdf_mat = bsdf->pdf(bRec);

                Color3f brdf = bsdf->eval(bRec);
                w_em = (pdf_mat + pdf_em) > 0 ? pdf_em / (pdf_mat + pdf_em) : 0.f;
                L += w_em * t * n_lights * brdf * Li * std::max(cos, 0.f);
            }


            BSDFQueryRecord bRec(its.toLocal(-contin_ray.d), its.uv);
            Color3f brdf = bsdf->sample(bRec, sampler->next2D());
            if (brdf == 0.f)
                break;
            pdf_mat = bsdf->pdf(bRec);
            t *= brdf;
            bounces++;
            contin_ray = Ray3f(its.p, its.toWorld(bRec.wo), Epsilon, INFINITY);

            if (!scene->rayIntersect(contin_ray, its))
            {
                if (scene->getEnvLights() != nullptr)
                {
                    EmitterQueryRecord lRec;
                    lRec.wi = its.toWorld(bRec.wo);

                    L+=t* scene->getEnvLights()->eval(lRec);
                    return L;
                }
                else return L;
            }

            if (its.mesh->isEmitter())
            {
                auto area_emitter = its.mesh->getEmitter();
                EmitterQueryRecord lRec(contin_ray.o, its.p, its.shFrame.n);
                pdf_em = area_emitter->pdf(lRec) / n_lights;
                w_mat = (pdf_mat + pdf_em) > 0 ? pdf_mat / (pdf_mat + pdf_em) : 0.f;
                if (bRec.measure == EDiscrete)
                {
                    w_mat = 1.f;
                }
                L += w_mat * t * its.mesh->getEmitter()->eval(lRec);
            }
        }
        return L;
    }



    std::string toString() const {
        return "Path_MIS";
    }
};


NORI_REGISTER_CLASS(Path_MIS, "path_mis");
NORI_NAMESPACE_END

/*
#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
NORI_NAMESPACE_BEGIN


class Path_MIS :public Integrator
{
public:
    Path_MIS(const PropertyList& props) {}
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        int bounces = 0;
        Color3f L(0.f);
        Color3f t(1.f);
        Ray3f contin_ray = ray;
        int n_lights = scene->getLights().size();
        float pdf_em=0, pdf_mat=0;
        float w_mat = 1, w_em = 1;
        bool specular = false;

        if (!scene->rayIntersect(contin_ray, its))
        {
            if (scene->getEnvLights() != nullptr)
            {
                EmitterQueryRecord lRec;
                lRec.wi = ray.d;
                return scene->getEnvLights()->eval(lRec);
            }
            else return L;
        }

        while (true)
        {

            if (bounces >= 3)
            {
                float q = std::min(t.maxCoeff(), 0.99f);
                if (sampler->next1D() > q)
                    break;
                else t /= q;
            }
            
            if( bounces==0 || specular==true)
            {
                if (its.mesh->isEmitter())
                {
                    EmitterQueryRecord lRec(contin_ray.o, its.p, its.shFrame.n);
                    L += t * its.mesh->getEmitter()->eval(lRec);
                }
            }


            //emitter
            auto bsdf = its.mesh->getBSDF();
            auto emitter = scene->getRandomEmitter(sampler->next1D());
            bool BSSRDF;
            if (bsdf != nullptr)
                BSSRDF = bsdf->isBSSRDF();
            else BSSRDF = false;

            EmitterQueryRecord lRec(its.p);
            Color3f Li = emitter->sample(lRec, sampler->next2D());
            pdf_em = emitter->pdf(lRec) / n_lights;
            if (!scene->rayIntersect(lRec.shadowRay))
            {
                BSDFQueryRecord bRec(its.toLocal(-contin_ray.d), its.toLocal(lRec.wi), ESolidAngle, its.uv);
                if (BSSRDF)
                {
                    Frame f(its.shFrame.n);
                    float channel = sampler->next1D();
                    int channel_idx;
                    if (channel < 1.f / 3.f)
                        channel_idx = 0;
                    else if (channel < 2.f / 3.f)
                        channel_idx = 1;
                    else channel_idx = 2;

                    BSSRDFQuery bssrdfRec(its.shFrame.n,its.p,its.p,channel_idx);                    
                    bRec.bssrdfRec = bssrdfRec;                    
                }

                float cos = its.shFrame.cosTheta(its.toLocal(lRec.wi));
                pdf_mat = bsdf->pdf(bRec);
                
                Color3f brdf = bsdf->eval(bRec);
                w_em = (pdf_mat + pdf_em) > 0 ? pdf_em / (pdf_mat + pdf_em) : 0.f;
                L += w_em * t * n_lights * brdf * Li * std::max(cos, 0.f);
            }    

            
            BSDFQueryRecord bRec(its.toLocal(-contin_ray.d), its.uv);
            if (BSSRDF)
            {
                Frame f(its.shFrame.n);
                float channel = sampler->next1D();
                int channel_idx;
                if (channel < 1.f / 3.f)
                    channel_idx = 0;
                else if (channel < 2.f / 3.f)
                    channel_idx = 1;
                else channel_idx = 2;
                Point2f rnd = sampler->next2D();

                BSSRDFQuery bssrdfRec(f, rnd, channel_idx, its.mesh, scene,its.p);
                bRec.bssrdfRec = bssrdfRec;
            }
            Color3f brdf = bsdf->sample(bRec, sampler->next2D());
            if (brdf == 0.f)
                break;
            if (bRec.measure == EDiscrete)
                specular = true;
            else
                specular = false;
            pdf_mat = bsdf->pdf(bRec);
            t *= brdf;
            bounces++;
            if (BSSRDF)
            {
                Frame f(bRec.bssrdfRec.s_normal);
                contin_ray = Ray3f(bRec.bssrdfRec.scatterPoint,
                    f.toWorld(bRec.wo));
            }
            else contin_ray = Ray3f(its.p, its.toWorld(bRec.wo));

            if (!scene->rayIntersect(contin_ray, its))
            {
                if (scene->getEnvLights() != nullptr)
                {
                    EmitterQueryRecord lRec;
                    lRec.wi = its.toWorld(bRec.wo);
                    return scene->getEnvLights()->eval(lRec);
                }
                else return L;
            }
            
            if (its.mesh->isEmitter())
            {
                auto area_emitter = its.mesh->getEmitter();
                EmitterQueryRecord lRec(contin_ray.o, its.p, its.shFrame.n);
                pdf_em = area_emitter->pdf(lRec)/n_lights;
                w_mat = (pdf_mat + pdf_em) > 0 ? pdf_mat / (pdf_mat + pdf_em) : 0.f;
                if (bRec.measure == EDiscrete)
                {
                    w_mat = 1.f;
                }
                L += w_mat * t * its.mesh->getEmitter()->eval(lRec);
            }            
        }
        return L;
    }



    std::string toString() const {
        return "Path_MIS";
    }
};


NORI_REGISTER_CLASS(Path_MIS, "path_mis");
NORI_NAMESPACE_END
*/