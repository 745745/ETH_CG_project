#include <nori/medium.h>


NORI_NAMESPACE_BEGIN

class HomogeneousMedium : public Medium {
public:
    HomogeneousMedium(const PropertyList& props) {
        m_sigma_a = props.getColor("sigmaS", Color3f(0.0014f, 0.0001f, 0.00002f));
        m_sigma_s = props.getColor("sigmaA", Color3f(0.002f, 0.003f, 0.004f));
        g = props.getFloat("g", 0.f);
        sigma_t = m_sigma_a + m_sigma_s;
        phase =new Phase(g);
    }

    Color3f transmittance(float dist) const override {
        return Color3f(std::exp(-sigma_t.r() * dist), std::exp(-sigma_t.g() * dist), std::exp(-sigma_t.b() * dist));
    }

    Color3f sample(const Ray3f& ray, Point2f& sample, MediumInteraction& mi, float& maxDist) const override {
        int idx;
        if (sample[0] < 1.f / 3.f)
            idx = 0;
        else if (sample[0] < 2.f / 3.f)
            idx = 1;
        else idx = 2;

        float t = -std::log(1 - sample[1]) / sigma_t[idx];
        t = t * ray.d.norm();
        mi.medium = (Medium*)this;
        // scatter
        if (t < maxDist) {
            mi.t = t;
            mi.p = ray.o + t * ray.d;
            
            mi.scatter = true;
            Color3f tr = transmittance(t);
            Color3f density = sigma_t * tr;
            float pdf=0;
            for (int i = 0; i < 3; i++)
            {
                pdf += density[i];
            }
            pdf /= 3.f;
            return tr*m_sigma_s/pdf;
        }
        // transmission
        else
        {
            mi.scatter = false;
            mi.p= ray.o + maxDist * ray.d;
            mi.t = maxDist;
            Color3f tr = transmittance(maxDist);
            float pdf = 0;
            for (int i = 0; i < 3; i++)
            {
                pdf += tr[i];
            }
            pdf /= 3.f;
            return tr / pdf;
        }
    }

    std::string toString() const override {
        return tfm::format(
            "HomogeneousMedium[\n"
            "  sigmaS = %s,\n"
            "  sigmaA = %s\n"
            "]",
            m_sigma_s.toString(),
            m_sigma_a.toString()
        );
    }

private:
    float g;
};

NORI_REGISTER_CLASS(HomogeneousMedium, "homogeneous");
NORI_NAMESPACE_END