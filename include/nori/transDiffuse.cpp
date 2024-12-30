#include <nori/bsdf.h>
#include <nori/frame.h>
#include <Eigen/Dense>
#include <nori/warp.h>
#include<nori/BSSRDF.h>
#include<nori/transDiffuse.h>


NORI_NAMESPACE_BEGIN

TransDiffuse::TransDiffuse(const PropertyList& propList)
{
    /* Interior IOR (default: BK7 borosilicate optical glass) */
    m_intIOR = propList.getFloat("intIOR", 1.3046f);

    /* Exterior IOR (default: air) */
    m_extIOR = propList.getFloat("extIOR", 1.000277f);

    /* Albedo of the diffuse base material (a.k.a "kd") */
    m_albedo = propList.getColor("albedo", Color3f(0.9f));

}

Color3f TransDiffuse::eval(const BSDFQueryRecord& bRec) const
{
    if (bRec.wi.dot(bRec.wo) < 0)
        return 0.f;

    float cosThetaI = Frame::cosTheta(bRec.wi);


    if (Frame::cosTheta(bRec.wo) > 0) {
        return m_albedo * INV_PI;
    }
    
}

float TransDiffuse::pdf(const BSDFQueryRecord& bRec) const
{
    if (bRec.wi.dot(bRec.wo) < 0)
        return 0.f;
    else
        return INV_PI * abs(Frame::cosTheta(bRec.wo));
}

Color3f TransDiffuse::sample(BSDFQueryRecord& bRec, const Point2f& sample) const
{
    Normal3f n(0, 0, 1);
    float cos = Frame::cosTheta(bRec.wi);

    float cosThetaT = fresnel(cos, m_extIOR, m_intIOR);
    float eta;
    if (cos < 0)
    {
        eta = m_intIOR / m_extIOR;
        n = Normal3f(0, 0, -1);
    }
    else
    {
        eta = m_extIOR / m_intIOR;
        n = Normal3f(0, 0, 1);
    }

    if (sample.x() < cosThetaT) {
        bRec.wo = Warp::squareToCosineHemisphere(sample);
        bRec.measure = EMeasure::ESolidAngle;
        bRec.eta = 1.0f;
        float pdf = Warp::squareToCosineHemispherePdf(bRec.wo);
        return m_albedo  * INV_PI/pdf;
    }
    else {
        
        Point2f Sample = Point2f((sample.x() - cosThetaT) / (1.0f - cosThetaT), sample.y());
        float sin_I2 = 1.0f - cos * cos;
        float sin_T2 = sin_I2 * eta * eta;
        Vector3f wi = bRec.wi;
        bRec.measure = EMeasure::EDiscrete;
        if (sin_T2 >= 1)
        {
            bRec.eta = 1;
            bRec.wo = Vector3f(-wi.x(), -wi.y(), wi.z());
            return Color3f(1.f);
        }
        else
        {
            float cos_T = sqrt(1 - sin_T2);
            bRec.wo = (eta * (-wi) + (eta * wi.dot(n) - cos_T) * n).normalized();
            bRec.eta = eta;
            return  (Color3f(1.0f) - m_albedo)  *eta*eta;           
        }

    }
}

std::string TransDiffuse::toString() const
{
	return std::string();
}

NORI_REGISTER_CLASS(TransDiffuse, "td");
NORI_NAMESPACE_END