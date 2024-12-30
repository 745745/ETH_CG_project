/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Prévost

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Sphere : public Shape {
public:
    Sphere(const PropertyList & propList) {
        m_position = propList.getPoint3("center", Point3f());
        m_radius = propList.getFloat("radius", 1.f);

        m_bbox.expandBy(m_position - Vector3f(m_radius));
        m_bbox.expandBy(m_position + Vector3f(m_radius));

        if (propList.has("int_medium"))
            int_medium_name = propList.getString("int_medium");
        if (propList.has("ext_medium"))
            ext_medium_name = propList.getString("ext_medium");
    }

    virtual BoundingBox3f getBoundingBox(uint32_t index) const override { return m_bbox; }

    virtual Point3f getCentroid(uint32_t index) const override { return m_position; }

    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {

        Vector3f oc = m_position - ray.o;
        float d = oc.dot(ray.d);
        if (d < 0.f)
            return false;
        // To avoid expensive norm computation for d<0. Normalize it after rejection.
        d = d / ray.d.norm();

        float projection_square = oc.squaredNorm() - d * d;
        if (projection_square > pow(m_radius, 2))
            return false;

        float half_chord =sqrtf(m_radius* m_radius - projection_square);
        float t1 = (d - half_chord) / ray.d.norm();
        float t2 = (d + half_chord) / ray.d.norm();

        if (t1 >= ray.mint && t1 <= ray.maxt)
        {
            t = t1 ;
            return true;
        }
            
        else if (t2 >= ray.mint  && t2 <= ray.maxt )
        {
            t = t2 ;
            return true;
        }
        return false;
    }

    virtual void setHitInformation(uint32_t index, const Ray3f& ray, Intersection& its) const override {

            its.p = ray.o + its.t * ray.d;
            Vector3f normals = (its.p - m_position).normalized();
            Frame x(normals);
            its.geoFrame = its.shFrame = x;

            Point2f uv = sphericalCoordinates(normals);
            its.uv = Point2f(uv.y() / (2 * M_PI), uv.x() / M_PI);
            return; 
    }

    virtual void sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const override {
        Vector3f q = Warp::squareToUniformSphere(sample);
        sRec.p = m_position + m_radius * q;
        sRec.n = q;
        sRec.pdf = std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }
    virtual float pdfSurface(const ShapeQueryRecord & sRec) const override {
        return std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }


    virtual std::string toString() const override {
        return tfm::format(
                "Sphere[\n"
                "  center = %s,\n"
                "  radius = %f,\n"
                "  bsdf = %s,\n"
                "  emitter = %s\n"
                "]",
                m_position.toString(),
                m_radius,
                m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
                m_emitter ? indent(m_emitter->toString()) : std::string("null"));
    }

protected:
    Point3f m_position;
    float m_radius;
};

NORI_REGISTER_CLASS(Sphere, "sphere");
NORI_NAMESPACE_END
