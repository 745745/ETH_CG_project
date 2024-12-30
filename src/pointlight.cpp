#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter
{
public:
	Vector3f position;
	Color3f power;
	EMeasure EDiscrete;

	PointLight(const PropertyList& props)
	{
		position = props.getPoint3("position");
		power = props.getColor("power");
		isDelta = true;
	}
	
	Color3f sample(EmitterQueryRecord& lRec, const Point2f& sampler) const
	{
		lRec.measure = EDiscrete;
		lRec.p = position;
		lRec.wi = (position - lRec.ref).normalized();
		lRec.pdf = 1.f;
		lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, (position - lRec.ref).norm()- Epsilon);
		float distance = (position - lRec.ref).squaredNorm();
		return power / (4.f * M_PI * distance);
	};

	Color3f eval(const EmitterQueryRecord& lRec) const
	{
		return power / (4.f * M_PI * (position - lRec.ref).squaredNorm());
	};

	float pdf(const EmitterQueryRecord& lRec) const
	{
		return 1.f;
	};


	Color3f samplePhoton(Ray3f& ray, const Point2f& sample1, const Point2f& sample2) const
	{
		Vector3f direction = Warp::squareToUniformSphere(sample1);
		ray = Ray3f(position, direction);
		return power / (4.f * M_PI);
	};

	std::string toString() const {
		return "PointLight[]";
	}
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END
