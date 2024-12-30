/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob, Romain Prévost

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

#if !defined(__NORI_BSDF_H)
#define __NORI_BSDF_H

#include <nori/object.h>
#include <nori/scene.h>
#include <nori/frame.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN


struct BSSRDFQuery {
    // below has to be set before use
    // local frame of incident point
    Frame frame;
    Point2f rnd;

    // random idx for channel selection(need setup before using)
    int channel_idx;

    // mesh of incident point located in
    const Shape* mesh;
    const Scene* scene;

    // below should be set in sample
    // the normal of scatter point
    Normal3f s_normal;
    Point3f scatterPoint;
    Point3f incidentPoint;
    float pdf_axis;
    int choose_axis;

    // for sampling
    BSSRDFQuery(const Frame& frame,const Point2f& rnd,const int& channel_idx,
        const Shape* mesh, const Scene* scene, const Point3f& incidentPoint) :frame(frame), rnd(rnd), channel_idx(channel_idx),
        mesh(mesh), scene(scene), incidentPoint(incidentPoint){};

    // for pdf&eval
    BSSRDFQuery(const Normal3f& s_normal,const Point3f& scatter,const Point3f& incident,const int& channel_idx): s_normal(s_normal),scatterPoint(scatter),incidentPoint(incident),frame(Frame(s_normal)) , channel_idx(channel_idx) {};

    BSSRDFQuery() {};
};

/**
 * \brief Convenience data structure used to pass multiple
 * parameters to the evaluation and sampling routines in \ref BSDF
 */
struct BSDFQueryRecord {

    BSSRDFQuery bssrdfRec;

    /// Incident direction (in the local frame)
    Vector3f wi;

    /// Outgoing direction (in the local frame)
    Vector3f wo;

    /// Relative refractive index in the sampled direction
    float eta;

    /// Measure associated with the sample
    EMeasure measure;

    /// UV associated with the point
    Point2f uv;

    /// Create a new record for sampling the BSDF
    BSDFQueryRecord(const Vector3f &wi, const Point2f &uv)
        : wi(wi), measure(EUnknownMeasure), uv(uv) { }

    /// Create a new record for querying the BSDF
    BSDFQueryRecord(const Vector3f &wi,
            const Vector3f &wo, EMeasure measure, const Point2f& uv)
        : wi(wi), wo(wo), measure(measure), uv(uv) { }

};

/**
 * \brief Superclass of all bidirectional scattering distribution functions
 */
class BSDF : public NoriObject {
public:
    /**
     * \brief Sample the BSDF and return the importance weight (i.e. the
     * value of the BSDF * cos(theta_o) divided by the probability density
     * of the sample with respect to solid angles).
     *
     * \param bRec    A BSDF query record
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The BSDF value divided by the probability density of the sample
     *         sample. The returned value also includes the cosine
     *         foreshortening factor associated with the outgoing direction,
     *         when this is appropriate. A zero value means that sampling
     *         failed.
     */
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const = 0;

    /**
     * \brief Evaluate the BSDF for a pair of directions and measure
     * specified in \code bRec
     *
     * \param bRec
     *     A record with detailed information on the BSDF query
     * \return
     *     The BSDF value, evaluated for each color channel
     */
    virtual Color3f eval(const BSDFQueryRecord &bRec) const = 0;

    /**
     * \brief Compute the probability of sampling \c bRec.wo
     * (conditioned on \c bRec.wi).
     *
     * This method provides access to the probability density that
     * is realized by the \ref sample() method.
     *
     * \param bRec
     *     A record with detailed information on the BSDF query
     *
     * \return
     *     A probability/density value expressed with respect
     *     to the specified measure
     */

    virtual float pdf(const BSDFQueryRecord &bRec) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    virtual EClassType getClassType() const override { return EBSDF; }

    /**
     * \brief Return whether or not this BRDF is diffuse. This
     * is primarily used by photon mapping to decide whether
     * or not to store photons on a surface
     */
    virtual bool isDiffuse() const { return false; }


};



NORI_NAMESPACE_END

#endif /* __NORI_BSDF_H */
