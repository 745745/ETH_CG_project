#pragma once
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <Eigen/Dense>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class BSSRDF : public BSDF {
public:
	
	float eta;
	float g;
	float A;

	// v and R_m are used for sampling
	Color3f R_m;

	float pdf_channel = 1.f / 3.f;

	float C1, C1_inv;
	float C2;
	float Cphi, Ce;
	// max intersection in point sampling
	int max_intersection = 20;

	Color3f m_sigma_a, sigma_tr, sigma_s,sigma_t;
	Color3f s_prime, t_prime, alpha_prime, rho;
	Color3f D;
	Color3f z_r, z_v, z_e;
	Color3f albedo;

	float m_intIOR, m_extIOR;
	BSSRDF(const PropertyList& propList);
	Color3f BSSRDF::compute_Rd(float r) const;
	virtual Color3f eval(const BSDFQueryRecord&) const override;
	virtual float pdf(const BSDFQueryRecord& bRec) const override;
	float sample_point(BSDFQueryRecord& bRec, const Point2f& sample) const;
	virtual Color3f sample(BSDFQueryRecord& bRec, const Point2f& sample) const override;
	virtual std::string toString() const;
	virtual float disk_pdf(float r,float R) const;
	float mis_weight(BSDFQueryRecord& bRec, float& pointpdf) const;
};






NORI_NAMESPACE_END