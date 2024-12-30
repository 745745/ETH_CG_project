#include <nori/bsdf.h>
#include <nori/frame.h>
#include <Eigen/Dense>
#include <nori/warp.h>
#include<nori/BSSRDF.h>

NORI_NAMESPACE_BEGIN

float mis(float x, float y, float z)
{
	float r2 = y / x;
	float r3 = z / x;

	return 1.f / (1 + r2 * r2 + r3 * r3);
}

// from paper "A Practical Model for Subsurface Light Transport"
BSSRDF::BSSRDF(const PropertyList& propList) {

	/* Interior IOR (default: glass) */
	m_intIOR = propList.getFloat("intIOR", 1.5046f);

	/* Exterior IOR (default: air) */
	m_extIOR = propList.getFloat("extIOR", 1.000277f);

	eta = m_intIOR / m_extIOR;

	/* sigma_a (default: skin 1)*/
	m_sigma_a = propList.getColor("sigma_a", Color3f(0.0021, 0.0041, 0.0071));

	/* s_prime (default: skin 1)*/
	s_prime = propList.getColor("s_prime", Color3f(2.19, 2.62, 3.0));

	albedo = propList.getColor("albedo", Color3f(0.83 ,0.79, 0.75));

	sigma_s = albedo / (1 - albedo) * m_sigma_a;


	t_prime = s_prime + m_sigma_a;
	sigma_t = m_sigma_a + sigma_s;
	alpha_prime = s_prime / t_prime;

	if (eta < 1.f)
	{
		C1 = 0.919317 + eta * (-3.4793 + eta * (6.75335 + eta * (-7.80989 + eta * (4.98554 - eta * 1.36881))));
		C2 = (0.828421f
			- 2.62051f * eta
			+ 3.36231f * eta * eta
			- 1.95284f * eta * eta * eta
			+ 0.236494f * pow(eta, 4)
			+ 0.145787f * pow(eta, 5));
		float eta_inv = 1.f / eta;
		C1_inv = C1 = -9.23372 + eta_inv * (22.2272 + eta_inv * (-20.9292 + eta_inv * (10.2291 + eta_inv * (-2.54396 + eta_inv * 0.254913))));
	}
	else
	{
		C1 = -9.23372 + eta * (22.2272 + eta * (-20.9292 + eta * (10.2291 + eta * (-2.54396 + eta * 0.254913))));
		C2 = (-1641.1f
			+ 135.926f / pow(eta, 3)
			- 656.175f / pow(eta, 2)
			+ 1376.53f / eta
			+ 1213.67f * eta
			- 568.556f * eta * eta
			+ 164.798f * pow(eta, 3)
			- 27.0181f * pow(eta, 4)
			+ 1.91826f * pow(eta, 5));
		float eta_inv = 1.f / eta;
		C1_inv = C1 = 0.919317 + eta_inv * (-3.4793 + eta_inv * (6.75335 + eta_inv * (-7.80989 + eta_inv * (4.98554 - eta_inv * 1.36881))));
	}

	A = (1 + C2) / (1 - C1);
	D = (2.f * m_sigma_a + s_prime) / (3.f * t_prime * t_prime);

	sigma_tr = sqrt(m_sigma_a / D);
	z_e = -2 * A * D;
	z_r = 1.f / t_prime;

	z_v = -z_r + 2 * z_e;
	Cphi = (1.f - C1) * 0.25f;
	Ce = (1.f - C2) * 0.5f;
	R_m = -log(0.01f) / sigma_tr;
}



Color3f BSSRDF::compute_Rd(float r) const
{
	Color3f d_r = sqrt(r * r + z_r * z_r);
	Color3f d_v = sqrt(r * r + z_v * z_v);
	Color3f real_contribution = (Ce * z_r * (sigma_tr * d_r + 1.f) / (d_r * d_r) + Cphi/D) * exp(-sigma_tr * d_r) / d_r;
	Color3f virtual_contribution = (Ce * z_v * (sigma_tr * d_v + 1.f) / (d_v * d_v) + Cphi/D) * exp(-sigma_tr * d_v) / d_v;
	Color3f R_d = alpha_prime * alpha_prime * (real_contribution - virtual_contribution) * INV_FOURPI;
	return R_d;
}

/*  In "A Practical Model for Subsurface Light Transport", S = Sd+ S1.
 *  But according to "A better dipole", it points out that add S1 to S will not improve performance and itself is inaccurate. So S1 will be discarded here
 *  Besides, Sd's computation will change a little.
 *  To eval, bRec has to have wi, wo, scatter point, incident point
*/
Color3f BSSRDF::eval(const BSDFQueryRecord& bRec) const
{
	float r = (bRec.bssrdfRec.scatterPoint - bRec.bssrdfRec.incidentPoint).norm();
	float cos_i = abs(bRec.bssrdfRec.frame.n.dot(bRec.wi));
	float cos_o = abs(bRec.bssrdfRec.s_normal.dot(bRec.wo));
	float Ft_i = 1 - fresnel(cos_i, m_extIOR, m_intIOR);
	float Ft_o = 1 - fresnel(cos_o, m_extIOR, m_intIOR);
	//float Ft_i = 1.f;
	//float Ft_o = 1.f;
	// S_d

	if (r > R_m[bRec.bssrdfRec.channel_idx])
		return 0.f;
	Color3f R_d = compute_Rd(r);

	Color3f S_d = Ft_i * R_d * Ft_o / (1.f - C1);

	return S_d;

}

float BSSRDF::pdf(const BSDFQueryRecord& bRec) const {
	float r = (bRec.bssrdfRec.scatterPoint - bRec.bssrdfRec.incidentPoint).norm();
	if (r > R_m[bRec.bssrdfRec.channel_idx])
		return 0.f;


	float pdf_axis= bRec.bssrdfRec.pdf_axis;;
	Vector3f n = bRec.bssrdfRec.frame.n;
	Vector3f z = bRec.bssrdfRec.s_normal;
	float cos = abs(n.dot(z));


	float pdf_disk = disk_pdf(r, R_m[bRec.bssrdfRec.channel_idx]);

	return pdf_disk * pdf_axis * pdf_channel * cos;
}

float BSSRDF::sample_point(BSDFQueryRecord& bRec, const Point2f& sample) const
{
	float tr = sigma_tr[bRec.bssrdfRec.channel_idx];
	float r = -log(1.f - sample.x()) / tr;
	// reject sample
	if (r > R_m[bRec.bssrdfRec.channel_idx])
	{
		return 0.f;
	}

	float phi = 2.f * M_PI * sample.y();

	Vector3f s = bRec.bssrdfRec.frame.s;
	Vector3f t = bRec.bssrdfRec.frame.t;
	Vector3f n = bRec.bssrdfRec.frame.n;

	Frame bFrame;
	float axis_rand = bRec.bssrdfRec.rnd[0];
	if (axis_rand < 0.5f)
	{
		bFrame = Frame(s, t, n);
		bRec.bssrdfRec.pdf_axis = 0.5;
		bRec.bssrdfRec.choose_axis = 0;
	}
	else if (axis_rand < 0.75f)
	{
		bFrame = Frame(t, n, s);
		bRec.bssrdfRec.pdf_axis = 0.25;
		bRec.bssrdfRec.choose_axis = 1;
	}
	else
	{
		bFrame = Frame(n, s, t);
		bRec.bssrdfRec.pdf_axis = 0.25;
		bRec.bssrdfRec.choose_axis = 2;
	}


	Vector3f disk_point = Vector3f(r * cos(phi), r * sin(phi), 0);

	float pdf_disk = disk_pdf(r,R_m[bRec.bssrdfRec.channel_idx]);

	float h = sqrt(R_m[bRec.bssrdfRec.channel_idx] * R_m[bRec.bssrdfRec.channel_idx] - r * r);
	Vector3f hn = h * bFrame.n;
	Vector3f entry_pt = bRec.bssrdfRec.incidentPoint + bFrame.toWorld(disk_point) + hn;
	auto ray_dir = -n;
	float tmax = 2.f * h;
	std::vector<Intersection> intersects;
	Ray3f ray(entry_pt, ray_dir, Epsilon, tmax);
	// shoot prob ray to find possible points
	for (int i = 0; i < max_intersection; i++)
	{
		
		Intersection its;
		if (!bRec.bssrdfRec.scene->rayIntersect(ray, its))
		{
			break;
		}

		// make sure the scatter point is on the same mesh			
		if (its.mesh == bRec.bssrdfRec.mesh)
		{
			intersects.push_back(its);
		}
		ray.o = its.p;
		ray.maxt -= its.t;
		if (tmax < Epsilon)
			break;
	}

	if (intersects.size() == 0)
		return 0.f;

	float point_rnd = bRec.bssrdfRec.rnd[1];
	int idx = floor(point_rnd * intersects.size());
	bRec.bssrdfRec.scatterPoint = intersects[idx].p;
	bRec.bssrdfRec.s_normal = intersects[idx].shFrame.n;
	bRec.bssrdfRec.frame = intersects[idx].shFrame;
	float pdf = bRec.bssrdfRec.pdf_axis * pdf_disk * abs(bRec.bssrdfRec.s_normal.dot(bRec.bssrdfRec.frame.n))*pdf_channel;
	pdf /= intersects.size();
	return pdf;
}


/*
	using method from paper "BSSRDF Importance Sampling" and PBRT 15.4
*/
Color3f BSSRDF::sample(BSDFQueryRecord& bRec, const Point2f& sample) const {
	bRec.measure = EBSSRDF;

	float pdf_point = sample_point(bRec, sample);
	if (pdf_point == 0.f)
		return 0.f;

	// uniformly sample the direction?
	bRec.wo = Warp::squareToCosineHemisphere(sample);
	//bRec.wo = Warp::squareToUniformHemisphere(sample);
	Color3f bssrdf = eval(bRec);

	return  bssrdf / pdf_point;
}


std::string BSSRDF::toString() const {
	return tfm::format(
		"BSSRDF[\n"
		"  intIOR = %f,\n"
		"  extIOR = %f\n"
		"  g = %f\n"
		"]",
		m_intIOR, m_extIOR);
}

float BSSRDF::disk_pdf(float r,float R) const
{
	float pdf = 0.f;
	if (r > R)
		return 0.f;

	for (int i = 0; i < 3; i++)
	{
		float tr = sigma_tr[i];
		pdf += tr * exp(-tr * r);
	}
	pdf = pdf * INV_TWOPI / r;

	return pdf;
}

float BSSRDF::mis_weight(BSDFQueryRecord& bRec, float& pointpdf) const
{
	Vector3f d = -bRec.bssrdfRec.scatterPoint + bRec.bssrdfRec.incidentPoint;
	Frame f=bRec.bssrdfRec.frame;
	float du = abs(d.dot(f.s));
	float dv = abs(d.dot(f.t));

	float dot_un = abs(bRec.bssrdfRec.s_normal.dot(f.s));
	float dot_vn = abs(bRec.bssrdfRec.s_normal.dot(f.t));

	float pdf_u = disk_pdf(du,R_m[1]) * dot_un;
	float pdf_v = disk_pdf(dv, R_m[1])*dot_vn;

	if (pdf_u == 0.f || pdf_v == 0.f)
		return 1.f;

	switch (bRec.bssrdfRec.choose_axis)
	{
		case 0:
			return mis(pointpdf, 0.25f * pdf_u, 0.25f * pdf_v);
		case 1:
			return mis(pointpdf, 0.25f * pdf_u, 0.5f * pdf_v);
		case 2:
			return mis(pointpdf, 0.5f * pdf_u, 0.25f * pdf_v);
	}
}





NORI_REGISTER_CLASS(BSSRDF, "bssrdf");
NORI_NAMESPACE_END