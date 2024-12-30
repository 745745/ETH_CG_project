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
#include<nori/BSSRDF.h>
//#include <nori/warp.h>
//#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Shape::~Shape() {
    delete m_bsdf;
    //delete m_emitter; // scene is responsible for deleting the emitter
}

void Shape::activate() {

    if (!m_bsdf && m_bssrdf != nullptr)
    {
        /*
        Color3f albedo = m_bssrdf->albedo;
        PropertyList l;
        l.setColor("albedo", albedo);
        m_bsdf = static_cast<BSDF*>(
            NoriObjectFactory::createInstance("diffuse", l));
        m_bsdf->activate();
        */
        float m_ext = m_bssrdf->m_extIOR;
        float m_int = m_bssrdf->m_intIOR;
        PropertyList l;
        l.setFloat("intIOR", m_int);
        l.setFloat("extIOR", m_ext);
        l.setColor("albedo", m_bssrdf->albedo);
        m_bsdf = static_cast<BSDF*>(
            NoriObjectFactory::createInstance("dielectric", l));
        m_bsdf->activate();
    }
    if (!m_bsdf && int_medium_name.length()==0 && int_medium_name.length() == 0) {
        // If no material and medium was assigned, instantiate a diffuse BRDF 
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
        m_bsdf->activate();
    }
    

        
}


void Shape::setParent(NoriObject* parent)
{
    Scene* scene = static_cast<Scene*>(parent);
    if(int_medium_name.size()!=0)
        int_medium = scene->lookMedium(int_medium_name);
    if (ext_medium_name.size() != 0)
        ext_medium = scene->lookMedium(ext_medium_name);
}

void Shape::addChild(NoriObject *obj) {
    std::string t = obj->getIdName();
    switch (obj->getClassType()) {
        
    case EBSDF:
    {
        
        if (t == "bssrdf")
        {
            m_bssrdf = (BSSRDF*)(obj);
            break;
        }
        if (m_bsdf)
            throw NoriException(
                "Shape: tried to register multiple BSDF instances!");
        m_bsdf = static_cast<BSDF*>(obj);
        break;

    }
        case EEmitter:
            if (m_emitter)
                throw NoriException(
                    "Shape: tried to register multiple Emitter instances!");
            m_emitter = static_cast<Emitter *>(obj);
            m_emitter->setShape(static_cast<Shape*>(this));
            break;

        default:
            throw NoriException("Shape::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END
