#pragma once
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <Eigen/Dense>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class TransDiffuse : public BSDF {
public:
    TransDiffuse(const PropertyList& propList);


    /// Evaluate the BRDF for the given pair of directions
    virtual Color3f eval(const BSDFQueryRecord& bRec) const override;

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord& bRec) const override;

    /// Sample the BRDF
    virtual Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const override;

    virtual std::string toString() const override;
private:
    float m_intIOR, m_extIOR;
    Color3f m_albedo;
};
NORI_NAMESPACE_END


