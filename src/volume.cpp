#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/warp.h>
#include<nori/bsdf.h>
#include<nori/BSSRDF.h>
NORI_NAMESPACE_BEGIN


class Volume :public Integrator
{
public:
    Volume(const PropertyList& props) {}
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f L(0.f);
        Color3f q(1.f);
        MediumInteraction mi;
        Ray3f r = ray;
        r.d = r.d.normalized();
        int bounces = 0;
        bool specular = false;
        Intersection its;
        bool intersect = scene->rayIntersect(r, its, true);
        if (!intersect)
        {
            if (scene->getEnvLights() != nullptr)
            {
                EmitterQueryRecord lRec;
                lRec.wi = ray.d;
                return Tr(scene, r) * scene->getEnvLights()->eval(lRec);
            }
            else return L;
        }

        while (true)
        {
            intersect = scene->rayIntersect(r, its, false);
            if (!intersect)
            {
                if (scene->getEnvLights() != nullptr)
                {
                    EmitterQueryRecord lRec;
                    lRec.wi = ray.d;
                    return q*Tr(scene, r) * scene->getEnvLights()->eval(lRec);
                }
                else return L;
            }
            //ray is in medium, sample next point
            if (r.medium != nullptr)
            {
                q *= r.medium->sample(r, sampler->next2D(), mi, its.t);
            }
            if (q == 0.f)
                break;
            // if scattering happens
            if (mi.scatter)
            {
                L += q * sampleLight(scene, sampler, r.d, mi);
                Vector3f wo = mi.medium->phase->sample(-r.d, sampler->next2D());      
                r = Ray3f(mi.p, wo);
                r.medium = mi.medium;
                continue;
            }

            // no scattering, ray passes through medium or is not in medium
            if (intersect)
            {
                auto bsdf = its.mesh->getBSDF();

                if (its.mesh->isEmitter() && (bounces ==0 || specular))
                {
                    EmitterQueryRecord lRec(r.o, its.p, its.shFrame.n);
                    L += q * its.mesh->getEmitter()->eval(lRec);
                }

                // mesh doesn't have bsdf, which means ray passes out of medium
                if (bsdf == nullptr)
                {                                                   
                    L += q * sampleLight(scene, sampler, r.d, its, its.mesh->get_int_medium());
                    //Vector3f wo = its.mesh->get_int_medium()->phase->sample(-r.d, sampler->next2D());
                    r = Ray3f(its.p, r.d);
                    r.medium = new_medium(r.d, its);                                                            
                    continue;
                }
                else  // ray intersects with a bsdf surface
                {
                    L += q * sampleLight(scene, sampler, r.d, its, r.medium);
                    BSDFQueryRecord bRec(its.toLocal(-r.d), its.uv);
                    Color3f brdf = bsdf->sample(bRec, sampler->next2D());
                    if (brdf == 0.f)
                        break;
                    q *= brdf;
                    bounces++;
                    r = Ray3f(its.p, its.toWorld(bRec.wo));
                    r.medium = new_medium(r.d, its);
                    if (bRec.measure = EDiscrete)
                        specular = true;
                    else
                        specular = false;
                    float cos = its.shFrame.cosTheta(bRec.wo);
                    // if has bssrdf
                    if ( cos<0 && its.mesh->hasBSSRDF())
                    {
                        Frame f(its.shFrame.n);
                        float channel = sampler->next1D();
                        int channel_idx;
                        if (channel < 1.f / 3.f)
                            channel_idx = 0;
                        else if (channel < 2.f / 3.f)
                            channel_idx = 1;
                        else channel_idx = 2;
                        BSSRDFQuery bssrdfRec(f, sampler->next2D(), channel_idx, its.mesh, scene, its.p);
                        bRec.bssrdfRec = bssrdfRec;
                        Color3f t = its.mesh->getBSSRDF()->sample(bRec, sampler->next2D());
                        if (t == 0)
                            break;
                        q *= t;
                        Intersection its2;
                        its2.mesh = bRec.bssrdfRec.mesh;
                        its2.p = bRec.bssrdfRec.scatterPoint;
                        its2.shFrame = Frame(bRec.bssrdfRec.s_normal);
                        L += q * sampleLight(scene, sampler, its2.shFrame.toWorld(-bRec.wo), its2);

                        bRec= BSDFQueryRecord(-bRec.wo, its2.uv);
                        t = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                        if (t == 0.f)
                            break;
                        q *= t;
                        r = Ray3f(bRec.bssrdfRec.scatterPoint, its2.shFrame.toWorld(bRec.wo));
                        r.medium = new_medium(r.d, its2);
                    }

                }
               
            }

            if (bounces >= 3)
            {
                float t = std::min(q.maxCoeff(), 0.99f);
                if (sampler->next1D() > t)
                    break;
                else q /= t;
            }
        }
        return L;
    }


    // the point is a intersection point
    Color3f sampleLight(const Scene* scene, Sampler* sampler, const Vector3f wi, const Intersection& its,Medium* medium=nullptr) const
    {
        if (scene->getLights().size() == 0)
        {
            EmitterQueryRecord lRec;
            lRec.wi = wi;
            Color3f Le = scene->getEnvLights()->eval(lRec);
            Ray3f r(its.p, wi);
            Color3f tr = Tr(scene, r);
            return Le * tr;
        }
        const Emitter* light = scene->getRandomEmitter(sampler->next1D());
        return EstimateDirect(scene,sampler,wi,light,its, medium);
        //return 0.f; 
    }

    // the point is a scatter point
    Color3f sampleLight(const Scene* scene, Sampler* sampler, const Vector3f wi, const MediumInteraction& mi) const
    {
        if (scene->getLights().size() == 0)
        {
            EmitterQueryRecord lRec;
            lRec.wi = wi;
            Color3f Le = scene->getEnvLights()->eval(lRec);
            Ray3f r(mi.p, wi);
            Color3f tr = Tr(scene, r);
            return Le * tr;
        }
        const Emitter* light = scene->getRandomEmitter(sampler->next1D());
        Intersection its;
        its.p = mi.p;
        return EstimateDirect(scene, sampler, wi, light,its, mi.medium);
        //return 0.f;
    }

    // estimate L between a point and a emitter. If the point is in medium, pass medium
    Color3f EstimateDirect(const Scene* scene, Sampler* sampler,const Vector3f wi, const Emitter* light, const Intersection& its, Medium* medium = nullptr) const
    {        
        int n_lights = scene->getLights().size();
        Color3f L(0.f);
        EmitterQueryRecord lRec(its.p);
        Color3f Li = light->sample(lRec,sampler->next2D());
        float pdf_light = 0.f;
        float pdf_bsdf = 0.f;
        pdf_light = light->pdf(lRec)/n_lights;
        if (pdf_light > 0.f)
        {
            Color3f f;
            // point is on a mesh
            if (its.mesh!= nullptr && its.mesh->getBSDF()!=nullptr)
            {                
                BSDFQueryRecord bRec(its.toLocal( - wi), its.toLocal(lRec.wi), ESolidAngle, its.uv);
                float cos = its.shFrame.cosTheta(its.toLocal(lRec.wi));
                f = its.mesh->getBSDF()->eval(bRec) * std::max(cos, 0.f);                                   
                pdf_bsdf = its.mesh->getBSDF()->pdf(bRec);
            }
            else
            {
                // point is a scattering point
                float p = medium->phase->eval(-wi,lRec.wi);
                f = p;
                pdf_bsdf = p;
            }

            if ( !(f == 0.f) )
            {
                Li *= Tr(scene, its.p, lRec.p,medium);
                if (!(Li == 0.f))
                {
                    if (light->isDelta)
                        L += f * Li *n_lights;
                    else
                    {
                        float w = powerHeuristic(pdf_light, pdf_bsdf,2);
                        L += f * Li * w * n_lights;
                    }
                }
            }

        }
        
        // sample bsdf
        if (!light->isDelta)
        {
            Vector3f wo;
            Color3f f(0.f);
            bool isDiscrete = false;
            BSDFQueryRecord bRec(its.toLocal(-wi), its.uv);
            if (its.mesh != nullptr && its.mesh->getBSDF()!=nullptr)
            {                
                f = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                pdf_bsdf = its.mesh->getBSDF()->pdf(bRec);
                wo = its.toWorld(bRec.wo);
                isDiscrete = (bRec.measure == EDiscrete);                
            }
            else
            {
                wo = medium->phase->sample(wi, sampler->next2D());
                float p = medium->phase->eval(-wi, wo);
                pdf_bsdf = p;
                f = p;
            }


            if (pdf_bsdf > 0.f)
            {
                Color3f Le(0.f);
                Color3f tr(1.f);
                float w = 1;
                Ray3f reflectRay(its.p, wo);
                reflectRay.medium = medium;
                Intersection its2;
                //ignore medium here, only find the intersection on light
                bool intersect = scene->rayIntersect(reflectRay, its2);
                if (intersect)
                {
                    if (its2.mesh->isEmitter())
                    {
                        EmitterQueryRecord lRec(its.p, its2.p, its2.shFrame.n);
                        pdf_light = its2.mesh->getEmitter()->pdf(lRec) / n_lights;
                        Le = its2.mesh->getEmitter()->eval(lRec);
                        tr = Tr(scene, its.p, its2.p, reflectRay.medium);
                    }
                }
                else
                {
                    if (scene->getEnvLights() != nullptr)
                    {
                        EmitterQueryRecord lRec;
                        lRec.wi = reflectRay.d;
                        Le = scene->getEnvLights()->eval(lRec);
                        tr = Tr(scene, reflectRay);
                        L += f * Le * tr;
                        return L;
                    }
                }

                if (!isDiscrete)
                {                    
                    if (pdf_light == 0.f)
                        return L;
                    w = powerHeuristic(pdf_bsdf, pdf_light, 2);
                }

                if (!(Le == 0.f))
                    L += f * Le * tr * w;
            }
        }
           
        return L;
    }

    // Tr between two points. Return 0.f if occluded
    Color3f Tr(const Scene* scene,const Point3f& point1, const Point3f& point2,Medium* start_medium=nullptr) const
    {
        Ray3f ray(point1, (point2 - point1).normalized());
        ray.maxt = (point2 - point1).norm();
        ray.medium = start_medium;
        Color3f Tr(1.f);
        while (true)
        {
            Intersection its;
            if (scene->rayIntersect(ray, its, false))
            {
                // medium, compute tr               
                if (ray.medium != nullptr)
                {
                    Tr *= ray.medium->transmittance(its.t);
                } 
                
                ray.o = its.p;
                ray.maxt -= its.t;
                if (ray.maxt < Epsilon)
                    break;
                if (!its.mesh->isMedium()) // not medium, then ray is occluded
                {
                    return Color3f(0.f);
                }
                ray.medium = new_medium(ray.d, its);
            }
            else
                break;
        }
        return Tr;
    }

    // ray may not be occluded, then compute tr between ray and envlight
    Color3f Tr(const Scene* scene, Ray3f& ray) const
    {
        Color3f Tr(1.f);
        while (true)
        {
            Intersection its;
            if (scene->rayIntersect(ray, its, false))
            {
                // medium, compute tr               
                if (ray.medium != nullptr)
                {
                    Tr *= ray.medium->transmittance(its.t);
                }
                if (!its.mesh->isMedium()) // not medium, then ray is occluded
                    return Color3f(0.f);


                ray.o = its.p;
                ray.medium = new_medium(ray.d, its);
            }
            else
                break;
        }
        return Tr;
    }

    float powerHeuristic(float pdf1, float pdf2,float power=1) const
    {
        float f = pow(pdf1,power);
        float g = pow(pdf2, power);
        return f / (f + g);
    }

    Medium* new_medium(Vector3f d, Intersection& its) const
    {
        if (d.dot(its.shFrame.n) < 0)
            return its.mesh->get_int_medium();
        else
        {
            return its.mesh->get_ext_medium();
        }
    }

    std::string toString() const {
        return "Volume";
    }
};


NORI_REGISTER_CLASS(Volume, "volume");
NORI_NAMESPACE_END
